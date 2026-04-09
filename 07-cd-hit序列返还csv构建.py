import pandas as pd

# --- 配置 ---
fasta_file = "cd_hit_process.fasta"      # CD-HIT 输出的 fasta 文件
csv_file = "merged_final_optimized_replication_annotated.csv"  # 你的原始大表格
output_file = "final_data.csv"   # 最终结果

# --- 第一步：读取 FASTA 文件，提取 ID ---
target_ids = set()

print("正在解析 FASTA 文件...")
with open(fasta_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            # 1. 去掉 '>'
            raw_id = line.strip().lstrip('>')

            # 2. 核心清洗逻辑：
            # 你的截图显示 ID 格式为：NZ_OZ249305_19700_20312|idx=330500
            # 我们需要去掉 |idx=... 这一部分，只保留前面的部分进行匹配
            clean_id = raw_id.split('|')[0]

            target_ids.add(clean_id)

print(f"✅ 从 FASTA 中提取到 {len(target_ids)} 个有效 ID")
# 打印前 5 个看看样子对不对
print(f"   样本 ID: {list(target_ids)[:5]}")

# --- 第二步：读取 CSV 并筛选 ---
print("正在读取 CSV 并匹配...")
# 你的截图显示 CSV 第一列没有列名，或者是默认的 0, 1, 2...
# 如果第一行是列名，header=0；如果第一行就是数据，header=None
df = pd.read_csv(csv_file)

# 假设 CSV 的第一列列名就是 'ori_id' (根据你截图的第一行推测)
# 如果报错 KeyError，请尝试改成 df.iloc[:, 0] 来强制选择第一列
id_column = 'ori_id'

# 筛选：如果 CSV 中的 ori_id 在我们提取的 target_ids 集合中，就保留
matched_df = df[df[id_column].isin(target_ids)]

# --- 第三步：保存 ---
matched_df.to_csv(output_file, index=False)
print(f"🎉 完成！")
print(f"   原始 CSV 行数：{len(df)}")
print(f"   匹配到的行数：{len(matched_df)}")
print(f"   结果已保存至：{output_file}")