#!/usr/bin/env python3
"""快速合并 stringtie gene_abund.tab -> 表达矩阵（向量化 + 多进程，替代慢速 iterrows 版）"""
import pandas as pd, glob, os, sys
from multiprocessing import Pool
from collections import defaultdict

INPUTDIR = sys.argv[1]   # 03_Align_Filter（含 homo 子目录）
OUTDIR   = sys.argv[2]   # 输出目录（自动加时间戳）
QC_CSV   = sys.argv[3] if len(sys.argv) > 3 else None

ANNO_DIRS = ["mRNA/genecode/stringtie", "eRNA/EnhancerAtlas/stringtie",
             "eRNA/Ensembl/stringtie", "eRNA/FANTOM5/stringtie",
             "lncRNA/GENCODE/stringtie", "miRNA/miRBase/stringtie",
             "miRNA/MirGeneDB/stringtie", "lncRNA/NONCODE/stringtie"]

def load_passed(csv_path):
    if not csv_path or not os.path.exists(csv_path):
        return None
    df = pd.read_csv(csv_path)
    passed = set(df[df['Passed'] == 'Yes']['Sample_ID'])
    return passed if passed else None

def merge_one(anno):
    files = glob.glob(os.path.join(INPUTDIR, 'homo', 'GSE*', anno, 'SRR*', 'gene_abund.tab'))
    series = {}
    for f in files:
        srr = f.split(os.sep)[-2]
        if passed is not None and srr not in passed:
            continue
        df = pd.read_csv(f, sep='\t', usecols=['Gene ID', 'TPM'], index_col=0)
        if df.index.duplicated().any():
            df = df[~df.index.duplicated(keep='last')]
        series[srr] = df['TPM']
    if not series:
        return anno, None
    mat = pd.DataFrame(series)          # index=Gene ID, columns=SRR（自动对齐）
    mat = mat.fillna(0.0)
    mat = mat.sort_index(axis=0)
    mat = mat.sort_index(axis=1)
    safe = anno.replace('/', '_')
    out = os.path.join(OUTDIR, f"human_{safe}_matrix.csv")
    mat.to_csv(out, float_format='%.6g')
    return anno, mat.shape

if __name__ == '__main__':
    passed = load_passed(QC_CSV)
    print(f"passed samples: {len(passed) if passed else 'ALL'}", flush=True)
    os.makedirs(OUTDIR, exist_ok=True)
    with Pool(8) as pool:
        for anno, shape in pool.imap_unordered(merge_one, ANNO_DIRS):
            print(f"{anno}: {shape}", flush=True)
    print("DONE", flush=True)
