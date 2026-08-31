
import os
import argparse
import csv
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Step 3.3: Filter Alignment Quality")
    parser.add_argument("--inputdir", required=True, help="Directory containing results")
    parser.add_argument("--outputdir", required=True, help="Output directory for reports")
    parser.add_argument("--cutoff", type=float, default=70.0,
                        help="Alignment rate cutoff (default: 70.0%%)")
    return parser.parse_args()

def extract_alignment_rates(inputdir, outputdir, cutoff=70.0):
    print(f"Extracting alignment rates (cutoff={cutoff}%)...")
    results = []
    seen_samples = set()

    def _add_result(sample_id, log_file):
        """解析单个 QC 日志并加入结果（去重）"""
        if sample_id in seen_samples:
            return
        try:
            with open(log_file, "r") as f:
                content = f.read()
                matches = re.findall(r"(\d+\.\d+)%", content)
                if matches:
                    rate = float(matches[-1])
                    passed = "Yes" if rate >= cutoff else "No"
                    results.append({
                        "Sample_ID": sample_id,
                        "Alignment_Rate": rate,
                        "Passed": passed,
                        "Path": log_file
                    })
                    seen_samples.add(sample_id)
        except Exception as e:
            print(f"Failed to read log {log_file}: {e}")

    # 1) 递归查找所有 QC_results.log（hisat2file 目录，若尚未被清理）
    for root, dirs, files in os.walk(inputdir):
        if "QC_results.log" in files:
            log_file = os.path.join(root, "QC_results.log")
            # 路径结构 .../GSExxx/hisat2file/SampleID/QC_results.log
            sample_id = os.path.basename(root)
            _add_result(sample_id, log_file)

    # 2) 递归查找所有 QC_logs/ 目录（chunk 清理时保留的比对率日志，文件名为 SampleID.log）
    #    ★ 修复：QC_logs 可能位于 inputdir/QC_logs（旧结构）或
    #      每个 GSE 子目录下 inputdir/GSExxx/QC_logs（当前结构），需递归扫描
    for root, dirs, files in os.walk(inputdir):
        if os.path.basename(root) == "QC_logs":
            for fname in sorted(files):
                if fname.endswith(".log"):
                    sample_id = os.path.splitext(fname)[0]
                    _add_result(sample_id, os.path.join(root, fname))

    csv_file = os.path.join(outputdir, "alignment_quality.csv")
    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Sample_ID", "Alignment_Rate", "Passed", "Path"])
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Written {len(results)} records to {csv_file} (cutoff={cutoff}%)")
    return csv_file

def main():
    args = parse_args()
    os.makedirs(args.outputdir, exist_ok=True)
    csv_file = extract_alignment_rates(args.inputdir, args.outputdir, cutoff=args.cutoff)
    
    with open(os.path.join(args.outputdir, "Filter_finished.txt"), "w") as f:
        f.write(f"Filtering finished (cutoff={args.cutoff}%).\n")

if __name__ == "__main__":
    main()
