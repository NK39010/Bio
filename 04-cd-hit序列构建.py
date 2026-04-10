import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Build FASTA input for CD-HIT from annotated CSV.")
    parser.add_argument("--input", default="merged_final_optimized_replication_annotated.csv", help="Input CSV.")
    parser.add_argument("--output", default="ori.fasta", help="Output FASTA.")
    parser.add_argument("--min-len", type=int, default=0, help="Minimum sequence length.")
    parser.add_argument("--max-len", type=int, default=20000, help="Maximum sequence length.")
    return parser.parse_args()

def main():
    args = parse_args()

    print("📥 正在读取数据...")
    df = pd.read_csv(args.input)
    df = df.dropna(subset=['OriC sequence'])

    print(f"🔍 正在根据长度 ({args.min_len}-{args.max_len} bp) 进行筛选...")

    df['seq_len'] = df['OriC sequence'].astype(str).str.len()
    df_filtered = df[(df['seq_len'] >= args.min_len) & (df['seq_len'] <= args.max_len)]

    print(f"   原始序列数：{len(df)}")
    print(f"   筛选后序列数：{len(df_filtered)}")

    df_filtered = df_filtered.drop_duplicates(subset=['OriC sequence', 'ori_id'])
    print(f"   去重后序列数：{len(df_filtered)}")

    print(f"💾 正在写入 {args.output} ...")
    with open(args.output, "w") as f:
        for i, row in df_filtered.iterrows():
            ori_id = row['ori_id']
            seq = row['OriC sequence']
            fasta_id = f"{ori_id}|idx={i}"
            f.write(f">{fasta_id}\n{seq}\n")
    print("✅ 任务完成！")


if __name__ == "__main__":
    main()