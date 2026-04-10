import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Map CD-HIT FASTA IDs back to source CSV.")
    parser.add_argument("--fasta", default="cd_hit_process.fasta", help="CD-HIT output FASTA file.")
    parser.add_argument(
        "--csv",
        default="merged_final_optimized_replication_annotated.csv",
        help="Source annotated CSV file.",
    )
    parser.add_argument("--output", default="final_data.csv", help="Output filtered CSV file.")
    parser.add_argument("--id-column", default="ori_id", help="ID column in source CSV.")
    return parser.parse_args()

def main():
    args = parse_args()
    fasta_file = args.fasta
    csv_file = args.csv
    output_file = args.output

    target_ids = set()
    print("正在解析 FASTA 文件...")
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                raw_id = line.strip().lstrip('>')
                clean_id = raw_id.split('|')[0]
                target_ids.add(clean_id)

    print(f"✅ 从 FASTA 中提取到 {len(target_ids)} 个有效 ID")
    print(f"   样本 ID: {list(target_ids)[:5]}")

    print("正在读取 CSV 并匹配...")
    df = pd.read_csv(csv_file)

    id_column = args.id_column
    if id_column not in df.columns:
        raise KeyError(f"Column '{id_column}' not found in CSV: {csv_file}")

    matched_df = df[df[id_column].isin(target_ids)]
    matched_df.to_csv(output_file, index=False)
    print(f"🎉 完成！")
    print(f"   原始 CSV 行数：{len(df)}")
    print(f"   匹配到的行数：{len(matched_df)}")
    print(f"   结果已保存至：{output_file}")


if __name__ == "__main__":
    main()