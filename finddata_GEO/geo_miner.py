#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
geo_miner.py v2 — GEO 数据挖掘：按瘤种+建库类型关键词爬取高价值 GSE 数据集
（v2 新增：样本级标题核查 + 技术类型过滤 + 精准 GBM 判定词）

v1 教训（2026-08-12 实测）:
  - 摘要关键词分类不可靠：103 个"A 类"里只有 5 个高置信 GBM
  - 原因: ①glioma 关键词匹配到髓母/DIPG/室管膜瘤 ②SuperSeries 混入甲基化/ChIP ③细胞系处理实验混入
  - v2 改为: 用 esummary 的 samples 字段（样本标题）按占比判定 + gdstype 技术类型过滤

功能:
  1. 多关键词组合检索 + RNA-seq 策略 + 人类 + gse[Filter]
  2. 样本级核查分类:
     A  = 高置信瘤种组织队列（表达谱, 瘤种词样本≥70%）
     A? = 疑似（30-70%）
     B  = 非目标瘤种/液体活检/正常组织
     E  = 表达谱但细胞系/处理实验主导
     C  = 单细胞   D = 信息不足
  3. 输出: output/candidates_all_*.csv + bulk_tissue_top_*.csv + report_*.md
  4. 本地缓存（cache/），重跑秒出、限速友好

用法:
  python3 geo_miner.py --cancer gbm --min-samples 10
  python3 geo_miner.py --cancer luad --min-samples 50
  python3 geo_miner.py --gse GSE107559   # 单查一个 GSE（快速核查模式）

依赖: Python3 标准库，无需第三方包
"""

import os
import re
import csv
import json
import time
import argparse
import urllib.request
import urllib.parse
from collections import Counter
from datetime import datetime

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "output")
CACHE_DIR = os.path.join(BASE_DIR, "cache")
CACHE_FILE = os.path.join(CACHE_DIR, "geo_cache.json")
EMAIL = "3099928014@qq.com"
SLEEP = 0.4  # 限速间隔（秒）

# ─────────────────────────────────────────────
# 癌种 → 关键词组
# ─────────────────────────────────────────────
CANCER_KEYWORDS = {
    "gbm": ['glioblastoma', 'glioma', '"high-grade glioma"', 'IDH-wildtype'],
    "luad": ['"lung adenocarcinoma"', '"non-small cell lung cancer"', 'NSCLC'],
    "brca": ['"breast cancer"', 'breast carcinoma'],
    "crc":  ['"colorectal cancer"', '"colon cancer"', '"rectal cancer"'],
    "panc": ['"pancreatic cancer"', '"pancreatic ductal adenocarcinoma"', 'PDAC'],
}

# 目标瘤种的样本级判定词（小写正则片段）
CANCER_SAMPLE_WORDS = {
    "gbm": [r"glioblastoma", r"\bgbm\b", r"who grade ?iv", r"glioblastoma multiforme"],
    "luad": [r"lung adenocarcinoma", r"luad", r"nsclc", r"adenocarcinoma of the lung"],
    "brca": [r"breast cancer", r"breast carcinoma", r"\bbrca\b"],
    "crc":  [r"colorectal", r"colon cancer", r"rectal cancer", r"\bcrc\b"],
    "panc": [r"pancreatic cancer", r"pancreatic ductal", r"pdac", r"pancreas"],
}

# 泛瘤种词（目标瘤种词命中率低时兜底，如 "glioma" 泛称）
CANCER_FALLBACK_WORDS = {
    "gbm": [r"glioma", r"astrocytoma", r"oligodendroglioma", r"glial tumor", r"brain tumor", r"brain tumour"],
}

# 排除特征（样本级，按占比判定）
EXCLUDE_SAMPLE_FEATURES = {
    "液体活检/血液": [r"platelet", r"\bblood\b", r"serum", r"plasma", r"liquid biops", r"cfrna", r"\bpbmc\b"],
    "其他癌种": [r"medulloblastoma", r"ependymoma", r"\bdipg\b", r"cscc", r"squamous", r"melanoma",
               r"carcinoma", r"sarcoma", r"lymphoma", r"leukemia", r"breast", r"lung cancer", r"colon",
               r"lgg", r"pilocytic", r"low-grade", r"ganglioglioma"],
    "正常/非肿瘤组织": [r"\bnormal\b", r"astrocyte", r"neuron", r"ipsc", r"ips cell", r"embryonic",
                     r"hippocampus", r"cortex", r"white matter"],
    "细胞系/处理实验": [r"u87", r"u251", r"ln229", r"a172", r"t98g", r"gsc", r"neurosphere", r"cell line",
                     r"knockdown", r"knockout", r"crispr", r"sirna", r"shrna", r"sg[a-z]+", r"\bdox\b",
                     r"vehicle", r"dmso", r"treated", r"treatment", r"inhibitor", r"overexpression",
                     r"transfected", r"\brep\d", r"sh[0-9]", r"tmz", r"temozolomide", r"radiation"],
    "非人类模型": [r"mouse", r"\bmice\b", r"rat\b", r"zebrafish", r"\bsh\b"],
}

# 单细胞标记
SC_WORDS = [r"single-cell", r"single cell", r"snrna", r"scrna", r"single-nucleus", r"10x", r"spatial"]


# ─────────────────────────────────────────────
# NCBI API（带缓存 + 重试）
# ─────────────────────────────────────────────
def _fetch(url, retries=3, backoff=(10, 30, 60)):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": f"geo_miner/1.0 (mailto:{EMAIL})"})
            with urllib.request.urlopen(req, timeout=45) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as e:
            if attempt == retries - 1:
                raise
            wait = backoff[attempt] if attempt < len(backoff) else 60
            print(f"    [RETRY] 请求失败({e})，等待 {wait}s 后重试 ({attempt+2}/{retries})...")
            time.sleep(wait)


def _load_cache():
    if os.path.exists(CACHE_FILE):
        try:
            with open(CACHE_FILE) as f:
                return json.load(f)
        except Exception:
            pass
    return {"esearch": {}, "esummary": {}}


def _save_cache(cache):
    os.makedirs(CACHE_DIR, exist_ok=True)
    with open(CACHE_FILE, "w") as f:
        json.dump(cache, f, ensure_ascii=False)


def esearch(db, term, retmax=600, sort="relevance"):
    cache = _load_cache()
    key = f"{db}|{term}|{retmax}|{sort}"
    if key in cache["esearch"]:
        return cache["esearch"][key]
    url = (f"{EUTILS}/esearch.fcgi?db={db}&term={urllib.parse.quote(term)}"
           f"&retmax={retmax}&sort={sort}&retmode=json")
    try:
        raw = _fetch(url)
        d = json.loads(raw)
        res = d.get("esearchresult", {})
        if "count" not in res or "idlist" not in res:
            raise RuntimeError(f"esearch 异常响应: {raw[:150]}")
        ids = res.get("idlist", [])
        time.sleep(SLEEP)
        cache["esearch"][key] = ids
        _save_cache(cache)
        return ids
    except Exception as e:
        print(f"  [ERR] esearch 失败: {e}")
        return []


def esummary(db, ids):
    if not ids:
        return []
    cache = _load_cache()
    out, todo = [], []
    for uid in ids:
        if uid in cache["esummary"]:
            out.append(cache["esummary"][uid])
        else:
            todo.append(uid)
    for i in range(0, len(todo), 150):
        batch = todo[i:i + 150]
        url = f"{EUTILS}/esummary.fcgi?db={db}&id={','.join(batch)}&retmode=json"
        try:
            d = json.loads(_fetch(url))
            for uid in batch:
                if uid in d.get("result", {}):
                    info = d["result"][uid]
                    cache["esummary"][uid] = info
                    out.append(info)
            time.sleep(SLEEP)
        except Exception as e:
            print(f"  [ERR] esummary 失败: {e}")
    _save_cache(cache)
    return out


def get_gse_info(gse):
    """按 GSE 号查（用于单查模式）"""
    ids = esearch("gds", f"{gse}[Accession]", retmax=5)
    if not ids:
        return None
    sums = esummary("gds", ids)
    return sums[0] if sums else None


# ─────────────────────────────────────────────
# v2 分类器：样本级 + 技术类型
# ─────────────────────────────────────────────
def classify_v2(s, cancer="gbm"):
    """返回 (类别, 理由, gbm_pct, 排除特征, 样本示例)"""
    gdstype = (s.get("gdstype") or "").lower()
    samples = s.get("samples") or []
    n = len(samples) or s.get("n_samples", 0)
    sample_titles = [(x.get("title") or "") for x in samples]
    st = " ".join(sample_titles).lower()

    # ── 技术类型判定 ──
    is_expr = "expression profiling" in gdstype
    is_seq = "high throughput sequencing" in gdstype
    is_methyl = "methylation" in gdstype
    is_chip = ("binding" in gdstype) or ("occupancy" in gdstype) or ("chip" in gdstype)
    tech_note = []
    if not is_expr:
        tech_note.append("非表达谱" if not is_seq else "非Expression类型")
    if is_methyl:
        tech_note.append("甲基化")
    if is_chip:
        tech_note.append("ChIP/结合")

    # ── 样本级瘤种判定 ──
    gbm_pat = re.compile("|".join(CANCER_SAMPLE_WORDS.get(cancer, [])))
    gbm_n = sum(1 for t in sample_titles if gbm_pat.search(t.lower()))
    fall_pat = re.compile("|".join(CANCER_FALLBACK_WORDS.get(cancer, [])))
    fall_n = sum(1 for t in sample_titles if fall_pat.search(t.lower()))
    gbm_pct = gbm_n / n * 100 if n else 0

    # ── 样本级排除特征（占比）──
    feat_hits = {}
    for fname, pats in EXCLUDE_SAMPLE_FEATURES.items():
        pat = re.compile("|".join(pats))
        c = sum(1 for t in sample_titles if pat.search(t.lower()))
        if c and c / n >= 0.05:
            feat_hits[fname] = c

    # ── 单细胞优先 ──
    sc_pat = re.compile("|".join(SC_WORDS))
    sc_n = sum(1 for t in sample_titles if sc_pat.search(t.lower()))
    if sc_n / n >= 0.5 if n else False:
        return "C", "单细胞/空间数据", gbm_pct, feat_hits, (sample_titles[0] if sample_titles else "")[:45]

    # ── 判定 ──
    # 表达谱缺失 → 甲基化/ChIP 等
    if not is_expr and (is_methyl or is_chip or not is_seq):
        return "E2", f"非表达谱技术({';'.join(tech_note)})", gbm_pct, feat_hits, (sample_titles[0] if sample_titles else "")[:45]

    # 高置信目标瘤种：瘤种词样本占比 ≥70% 且无非目标主导、且非细胞系/处理实验主导
    cell_feat = feat_hits.get("细胞系/处理实验", 0)
    cell_ratio = cell_feat / n if n else 0
    other_feat = feat_hits.get("其他癌种", 0)
    if gbm_pct >= 70 and not (other_feat and other_feat / n > 0.3) and cell_ratio <= 0.3:
        return "A", f"瘤种词样本{gbm_pct:.0f}%", gbm_pct, feat_hits, (sample_titles[0] if sample_titles else "")[:45]

    # 瘤种词高但细胞系/处理实验主导 → 机制类（E）
    if gbm_pct >= 70 and cell_ratio > 0.3:
        return "E", f"瘤种词{gbm_pct:.0f}%但细胞系/处理主导({cell_feat})", gbm_pct, feat_hits, (sample_titles[0] if sample_titles else "")[:45]

    # 疑似
    if gbm_pct >= 30:
        return "A?", f"瘤种词样本{gbm_pct:.0f}%", gbm_pct, feat_hits, (sample_titles[0] if sample_titles else "")[:45]

    # 泛瘤种词（如 glioma 泛称）占比高但无明确目标瘤种词 → 疑似
    if fall_n / n >= 0.5 if n else False:
        return "A?", f"泛瘤种词样本{fall_n/n*100:.0f}%(非精确瘤种)", gbm_pct, feat_hits, (sample_titles[0] if sample_titles else "")[:45]

    # 表达谱但非目标瘤种
    if feat_hits:
        main_feat = max(feat_hits, key=lambda k: feat_hits[k])
        return "B", f"非目标瘤种/类型({main_feat}={feat_hits[main_feat]})", gbm_pct, feat_hits, (sample_titles[0] if sample_titles else "")[:45]

    # 表达谱、样本标题无特征词（可能是编号式样本）→ 用摘要兜底
    summ = (s.get("summary") or "").lower()
    if gbm_pat.search(summ):
        return "A?", "摘要含瘤种词(样本标题无特征)", gbm_pct, feat_hits, (sample_titles[0] if sample_titles else "")[:45]
    return "D", "样本标题/摘要无特征词", gbm_pct, feat_hits, (sample_titles[0] if sample_titles else "")[:45]


# ─────────────────────────────────────────────
# 单查模式
# ─────────────────────────────────────────────
def single_check(gse, cancer):
    print(f"查询 {gse} ...")
    s = get_gse_info(gse)
    if not s:
        print("未找到")
        return
    cat, reason, pct, feats, ex = classify_v2(s, cancer)
    print(f"  GSE: {s.get('accession')} | 样本: {s.get('n_samples')} | 类型: {s.get('gdstype')}")
    print(f"  标题: {(s.get('title') or '')[:80]}")
    print(f"  判定: {cat} | {reason} | 瘤种词样本占比: {pct:.0f}%")
    print(f"  排除特征: {feats if feats else '-'}")
    print(f"  样本示例: {ex}")


# ─────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description="GEO 数据挖掘 v2（样本级核查）")
    ap.add_argument("--cancer", default="gbm", help="癌种键名: gbm/luad/brca/crc/panc")
    ap.add_argument("--min-samples", type=int, default=10, help="候选最小样本数")
    ap.add_argument("--max-results", type=int, default=600, help="esearch 最大返回数")
    ap.add_argument("--gse", default=None, help="单查模式: 核查一个 GSE")
    args = ap.parse_args()

    if args.gse:
        single_check(args.gse.upper(), args.cancer)
        return

    kws = CANCER_KEYWORDS.get(args.cancer, [args.cancer])
    print(f"检索癌种: {args.cancer} | 关键词: {kws}")

    # 1. 多关键词检索合并
    all_ids, seen = [], set()
    for kw in kws:
        term = (f'{kw}[All Fields] AND "RNA-seq"[Strategy] '
                f'AND "Homo sapiens"[Organism] AND gse[Filter]')
        ids = esearch("gds", term, retmax=args.max_results)
        new = [i for i in ids if i not in seen]
        seen.update(new)
        all_ids.extend(new)
        print(f"  [{kw}] 命中 {len(ids)} (新增 {len(new)})")
        time.sleep(SLEEP)

    print(f"\n合并去重共 {len(all_ids)} 个 GSE，拉取详情...")
    sums = esummary("gds", all_ids)
    sums.sort(key=lambda s: -s.get("n_samples", 0))

    # 2. v2 分类
    rows = []
    for s in sums:
        cat, reason, pct, feats, ex = classify_v2(s, args.cancer)
        n = s.get("n_samples", 0)
        if n < args.min_samples and cat in ("A", "A?"):
            cat = "A?(小样本)"
        rows.append({
            "GSE": s.get("accession", "?"),
            "UID": s.get("uid", ""),
            "n_samples": n,
            "title": (s.get("title") or "").replace("\n", " ")[:80],
            "date": (s.get("pdat") or "")[:10],
            "gdstype": (s.get("gdstype") or "")[:60],
            "tumor_pct": round(pct, 1),
            "summary": (s.get("summary") or "").replace("\n", " ")[:150],
            "category": cat,
            "reason": reason,
            "exclude_features": "; ".join(f"{k}={v}" for k, v in feats.items()) or "-",
            "sample_example": ex,
        })

    # 3. 输出
    os.makedirs(OUT_DIR, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d")
    fields = list(rows[0].keys()) if rows else ["GSE"]
    f_all = os.path.join(OUT_DIR, f"candidates_all_{ts}.csv")
    with open(f_all, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # A 类精选（高置信，按样本量）
    a_rows = [r for r in rows if r["category"] == "A"]
    a_rows.sort(key=lambda r: -r["n_samples"])
    f_a = os.path.join(OUT_DIR, f"bulk_tissue_top_{ts}.csv")
    with open(f_a, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in a_rows:
            w.writerow(r)

    # 4. 报告
    cc = Counter(r["category"] for r in rows)
    rep = [f"# GEO 挖掘报告 v2：{args.cancer}（{ts}）", ""]
    rep.append(f"- 总命中 GSE: {len(rows)}")
    rep.append(f"- 分类: {dict(cc)}")
    rep.append(f"- A = 高置信{args.cancer}组织队列（瘤种词样本≥70%）")
    rep.append("")
    rep.append(f"## A 类精选（{len(a_rows)} 个，按样本量排序）")
    rep.append("")
    rep.append("| GSE | 样本 | 瘤种% | 日期 | 标题 | 排除特征 |")
    rep.append("|---|---|---|---|---|---|")
    for r in a_rows:
        rep.append(f"| {r['GSE']} | {r['n_samples']} | {r['tumor_pct']:.0f}% | {r['date']} | {r['title'][:50]} | {r['exclude_features'][:30]} |")
    rep.append("")
    rep.append("## A? 疑似（需人工确认）")
    rep.append("")
    for r in [x for x in rows if x["category"].startswith("A?")]:
        rep.append(f"- {r['GSE']} (n={r['n_samples']}, 瘤种{r['tumor_pct']:.0f}%): {r['title'][:55]} — {r['reason']}")
    rep.append("")
    rep.append("## 说明")
    rep.append("- 分类基于样本标题 + gdstype 技术类型，比摘要关键词可靠")
    rep.append("- 生存数据需到各 GSE 页面人工确认")
    rep.append("- 单查: python3 geo_miner.py --gse GSE107559")
    rep.append(f"\n生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    f_rep = os.path.join(OUT_DIR, f"report_{ts}.md")
    with open(f_rep, "w", encoding="utf-8") as f:
        f.write("\n".join(rep))

    # 控制台
    print(f"\n分类统计: {dict(cc)}")
    print(f"\n{'='*100}\nA 类精选（高置信，按样本量排序）:")
    print(f"{'GSE':<12}{'样本':<7}{'瘤种%':<7}{'日期':<12} 标题")
    print("-" * 100)
    for r in a_rows:
        print(f"{r['GSE']:<12}{r['n_samples']:<7}{r['tumor_pct']:.0f}%{'':<3}{r['date']:<12} {r['title'][:55]}")
    print(f"\n✅ 输出:")
    print(f"  {f_all}")
    print(f"  {f_a}")
    print(f"  {f_rep}")


if __name__ == "__main__":
    main()
