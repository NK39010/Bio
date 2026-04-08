import pandas as pd

# === 配置参数 ===
INPUT_FILE = "merged_final_optimized.csv"
OUTPUT_FILE = "ori.fasta"
MIN_LEN = 0
MAX_LEN = 20000

# === 1. 读取数据 ===
print("📥 正在读取数据...")
df = pd.read_csv(INPUT_FILE)

# === 2. 基础清洗 ===
# 去掉空序列
df = df.dropna(subset=['OriC sequence'])

# === 3. 计算长度并筛选 (核心步骤) ===
print(f"🔍 正在根据长度 ({MIN_LEN}-{MAX_LEN} bp) 进行筛选...")

# 计算序列长度（转为字符串防止报错）
df['seq_len'] = df['OriC sequence'].astype(str).str.len()

# 筛选：保留长度在 200 到 5000 之间的行
df_filtered = df[(df['seq_len'] >= MIN_LEN) & (df['seq_len'] <= MAX_LEN)]

print(f"   原始序列数: {len(df)}")
print(f"   筛选后序列数: {len(df_filtered)}")

# === 4. 去重 (可选，推荐用于 CD-HIT 前处理) ===
# 注意：这里基于筛选后的数据进行去重
df_filtered = df_filtered.drop_duplicates(subset=['OriC sequence', 'ori_id'])
print(f"   去重后序列数: {len(df_filtered)}")

# === 5. 写入 FASTA ===
print(f"💾 正在写入 {OUTPUT_FILE} ...")
with open(OUTPUT_FILE, "w") as f:
    for i, row in df_filtered.iterrows():
        ori_id = row['ori_id']
        seq = row['OriC sequence']
        
        # ⭐ 关键：保证唯一ID
        fasta_id = f"{ori_id}|idx={i}"
        
        f.write(f">{fasta_id}\n{seq}\n")

print("✅ 任务完成！")