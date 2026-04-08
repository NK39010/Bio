import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ================= 配置区域 =================
clstr_file = "cd_hit_process95.fasta.clstr"  # 你的 .clstr 文件
output_file = "cdhit_redundancy_fixed95.png"

# 检查文件是否存在
if not os.path.exists(clstr_file):
    print(f"❌ 错误: 找不到文件 {clstr_file}")
    exit()

# ================= 1. 解析 .clstr 文件 =================
print(">>> 正在解析 .clstr 文件...")

cluster_sizes = []

with open(clstr_file, 'r') as f:
    current_size = 0
    for line in f:
        if line.startswith('>Cluster'):
            # 遇到新簇的标题行，保存上一个簇的大小
            if current_size > 0:
                cluster_sizes.append(current_size)
            current_size = 0  # 重置计数器
        else:
            # 序列行，计数 +1
            current_size += 1

    # 别忘了保存最后一个簇
    if current_size > 0:
        cluster_sizes.append(current_size)

# 转换为 Pandas Series 方便处理
sizes_series = pd.Series(cluster_sizes)
total_sequences = sizes_series.sum()
total_clusters = len(sizes_series)

print(f"解析完成: 共 {total_clusters} 个簇, 包含 {total_sequences} 条序列")

# ================= 2. 绘图 =================
sns.set_theme(style="whitegrid")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12

# 创建画布：1行2列
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# --- 图1: 簇大小分布 (柱状图) ---
# 统计每个大小出现的频次 (例如：大小为1的有500个，大小为2的有100个...)
size_counts = sizes_series.value_counts().sort_index()

# 为了图形美观，如果最大值太大，只显示前20个或者截取一部分
# 这里我们直接画，但如果数据很散，建议用 log 轴
sns.barplot(x=size_counts.index, y=size_counts.values, ax=ax1, color='#4c72b0')
ax1.set_title("Distribution of Cluster Sizes", fontsize=14)
ax1.set_xlabel("Sequences per Cluster (Redundancy Level)", fontsize=12)
ax1.set_ylabel("Number of Clusters", fontsize=12)
ax1.set_yscale('log')  # 使用对数坐标，防止大簇把小簇压扁看不见
ax1.set_xlabel("Cluster Size (1 = Unique)")

# --- 图2: 累积冗余度曲线 ---
# 按大小降序排列，看大簇的影响力
sorted_sizes = sizes_series.sort_values(ascending=False)
cumulative_sequences = sorted_sizes.cumsum()
cumulative_pct = (cumulative_sequences / total_sequences) * 100

ax2.plot(cumulative_pct.index, cumulative_pct.values, color='#55a868', linewidth=2)
ax2.fill_between(cumulative_pct.index, cumulative_pct.values, color='#55a868', alpha=0.2)

ax2.set_title("Cumulative Redundancy Impact", fontsize=14)
ax2.set_xlabel("Top N Largest Clusters", fontsize=12)
ax2.set_ylabel("Cumulative Data Coverage (%)", fontsize=12)

# --- 修复点：标记关键点 ---
# 计算前 10% 的簇包含了多少数据
top_10_pct_idx = int(len(cumulative_pct) * 0.1)  # 这里修正了变量名
y_val = cumulative_pct.iloc[top_10_pct_idx]

# 画点
ax2.scatter([top_10_pct_idx], [y_val], color='red', s=50, zorder=5)

# 加注释
ax2.annotate(f'Top 10% Clusters\ncontain {y_val:.1f}% data',
             xy=(top_10_pct_idx, y_val),
             xytext=(top_10_pct_idx + 50, y_val - 20),
             arrowprops=dict(arrowstyle='->', color='gray'),
             fontsize=10,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

plt.tight_layout()
plt.savefig(output_file, dpi=300)
print(f"✅ 图表已保存为: {output_file}")
plt.show()