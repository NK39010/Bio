import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# ====================== 读取你的 CSV ======================
df = pd.read_csv(r"C:\Users\MoSo\Desktop\data\merged_final_optimized.csv", low_memory=False)

# ====================== 只保留目标字段 ======================
keep_cols = [
    "OriC sequence",
    "ori_id",
    "plasmid_id",
    "species",
    "source",
    "pfamid_fast",
    "rep_id",
    "Rep_type_fast",
    "rep_seq",
    "rep_dna_seq",
    "full_replicon_seq",
    "split",          # 只保留，不填充
    "__index_level_0__"
]

df_out = df[keep_cols].copy()

# ====================== 序列标准化（核苷酸小写 / 氨基酸大写） ======================
def clean_seq(seq, is_aa=False):
    if pd.isna(seq):
        return ""
    seq = str(seq).strip()
    if is_aa:
        return seq.upper()
    else:
        return seq.lower()

df_out["OriC sequence"] = df_out["OriC sequence"].apply(lambda x: clean_seq(x, is_aa=False))
df_out["rep_seq"] = df_out["rep_seq"].apply(lambda x: clean_seq(x, is_aa=True))
df_out["rep_dna_seq"] = df_out["rep_dna_seq"].apply(lambda x: clean_seq(x, is_aa=False))
df_out["full_replicon_seq"] = df_out["full_replicon_seq"].apply(lambda x: clean_seq(x, is_aa=False))

# ====================== 不填充 split！原样保留 ======================
# 这里什么都不做

# ====================== 输出为 parquet ======================
table = pa.Table.from_pandas(df_out, preserve_index=False)
pq.write_table(table, r"C:\Users\MoSo\Desktop\data\merged_final_optimized.csv", compression="snappy")

print("✅ 转换完成！已输出：final_data95.parquet")
print(f"📊 共 {len(df_out)} 条数据 | split 列保持原样，未填充")