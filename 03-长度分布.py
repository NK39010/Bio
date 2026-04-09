import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 1. 读取数据
df = pd.read_csv("merged_final_optimized.csv")
df = df.dropna(subset=['OriC sequence'])
df['length'] = df['OriC sequence'].astype(str).str.len()

# 2. 绘图设置
fig, axes = plt.subplots(1, 2, figsize=(20, 7))

# --- 图 1：对数坐标直方图 (推荐) ---
# 使用对数坐标，能看清 100-10000bp 之间的细节
sns.histplot(data=df, x='length', bins=100, log_scale=(True, False), ax=axes[0], color='skyblue', edgecolor='black')
axes[0].set_title('Distribution of OriC Length (Log Scale X-axis)', fontsize=14)
axes[0].set_xlabel('Length (bp) - Log Scale', fontsize=12)
axes[0].set_ylabel('Count', fontsize=12)
axes[0].axvline(df['length'].median(), color='red', linestyle='--', label=f'Median: {df["length"].median()}')
axes[0].legend()

# --- 图 2：截断坐标箱线图 ---
# 强制把 Y 轴限制在 0-5000bp，把那些几百万的异常值切掉，只看主体
sns.boxplot(data=df, y='length', ax=axes[1], color='lightgreen')
axes[1].set_title('Boxplot (Zoomed in: 0-5000bp)', fontsize=14)
axes[1].set_ylabel('Length (bp)', fontsize=12)
axes[1].set_ylim(0, 10000) # ⭐ 关键：限制坐标轴范围，强制放大主体区域

plt.tight_layout()
plt.show()

# 3. 打印详细统计信息，辅助决策
print("📊 详细长度统计：")
print(df['length'].describe())