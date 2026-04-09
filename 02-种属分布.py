import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from collections import Counter

# ================= 配置区域 =================
csv_file = "merged_final_optimized.csv"
output_excel = "csv_species_statistics.xlsx"
output_image = "csv_species_analysis.png"

# 绘图展示的数量限制
TOP_N_GENUS_PLOT = 30
TOP_N_SPECIES_PLOT = 50
# ===========================================

# 设置绘图风格
sns.set(style="whitegrid")
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans'] 
plt.rcParams['axes.unicode_minus'] = False

# ================= 1. 数据读取与去重 =================
print(f"正在读取 CSV 文件：{csv_file} ...")

try:
    # 读取 CSV
    df_raw = pd.read_csv(csv_file)
    print(f"原始数据行数 (包含重复 Ori): {len(df_raw)}")

    # 检查是否有 'ori_id' 列
    if 'ori_id' not in df_raw.columns:
        print("错误：找不到 'ori_id' 列，无法进行去重。请检查输入文件。")
        exit()

    # 按照 'ori_id' 去重，保留第一条出现的记录
    df = df_raw.drop_duplicates(subset=['ori_id'], keep='first')
    
    print(f"去重后数据行数 (唯一 Ori): {len(df)}")
    print(f"去除了 {len(df_raw) - len(df)} 行重复数据。")

    # 检查是否有 'species' 列
    if 'species' not in df.columns:
        print(f"错误：未找到 'species' 列。现有列：{df.columns.tolist()}")
        exit()
        
    # 直接获取物种列表，不做清洗
    raw_species_list = df['species'].dropna().tolist() 

except FileNotFoundError:
    print(f"错误：找不到文件 '{csv_file}'")
    exit()

# ================= 2. 数据统计 (无清洗) =================
print("正在统计物种（保留原始名称）...")

# 直接使用原始列表进行统计
species_counter = Counter(raw_species_list)
total_count = len(raw_species_list) 

# 区分属 (Genus) 和 种 (Species)
# 这里仅为了统计属，临时分割一下字符串，不影响原始物种名的显示
genus_counter = Counter()
for species_name, count in species_counter.items():
    if pd.isna(species_name):
        continue
    # 提取第一个单词作为属名
    genus = str(species_name).split()[0]
    genus_counter[genus] += count

print(f"统计完成。共识别出 {len(genus_counter)} 个属，{len(species_counter)} 个种。")

# ================= 3. 构建 DataFrame =================

# 构建属统计表
df_genus = pd.DataFrame.from_dict(genus_counter, orient='index', columns=['Count'])
df_genus = df_genus.sort_values(by='Count', ascending=False).reset_index()
df_genus.columns = ['Name', 'Count']
df_genus['Percentage'] = (df_genus['Count'] / total_count * 100).round(4)
df_genus['Cumulative_Count'] = df_genus['Count'].cumsum()
df_genus['Cumulative_Percentage'] = (df_genus['Cumulative_Count'] / total_count * 100).round(4)
df_genus['Rank'] = df_genus.index + 1

# 构建种统计表 (使用原始名称)
df_species = pd.DataFrame.from_dict(species_counter, orient='index', columns=['Count'])
df_species = df_species.sort_values(by='Count', ascending=False).reset_index()
df_species.columns = ['Name', 'Count']
df_species['Percentage'] = (df_species['Count'] / total_count * 100).round(4)
df_species['Cumulative_Count'] = df_species['Count'].cumsum()
df_species['Cumulative_Percentage'] = (df_species['Cumulative_Count'] / total_count * 100).round(4)
df_species['Rank'] = df_species.index + 1

# 构建概览
summary_data = {
    'Metric': [
        'Total Unique OriC',
        'Unique Genera', 
        'Unique Species', 
        'Top 1 Genus Name',
        'Top 1 Genus Count', 
        'Top 1 Genus %', 
        'Top 10 Genera %', 
        'Top 1 Species Name',
        'Top 1 Species Count', 
        'Top 1 Species %', 
        'Top 10 Species %',
        'Top 50 Species %'
    ],
    'Value': [
        total_count,
        len(df_genus),
        len(df_species),
        df_genus.iloc[0]['Name'] if len(df_genus) > 0 else 'N/A',
        df_genus.iloc[0]['Count'] if len(df_genus) > 0 else 0,
        f"{df_genus.iloc[0]['Percentage']:.2f}%" if len(df_genus) > 0 else '0%',
        f"{df_genus.head(10)['Percentage'].sum():.2f}%" if len(df_genus) > 0 else '0%',
        df_species.iloc[0]['Name'] if len(df_species) > 0 else 'N/A',
        df_species.iloc[0]['Count'] if len(df_species) > 0 else 0,
        f"{df_species.iloc[0]['Percentage']:.2f}%" if len(df_species) > 0 else '0%',
        f"{df_species.head(10)['Percentage'].sum():.2f}%" if len(df_species) > 0 else '0%',
        f"{df_species.head(50)['Percentage'].sum():.2f}%" if len(df_species) > 0 else '0%'
    ]
}
df_summary = pd.DataFrame(summary_data)

# ================= 4. 导出 Excel =================
print(f"正在生成 Excel 报告：{output_excel} ...")

with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
    df_summary.to_excel(writer, sheet_name='Summary_概览', index=False)
    df_genus.to_excel(writer, sheet_name='Genus_属统计', index=False)
    df_species.to_excel(writer, sheet_name='Species_种统计', index=False)
    df_species.head(50).to_excel(writer, sheet_name='Top50_核心摘要', index=False)

    # 自动调整列宽
    workbook = writer.book
    for sheet_name in writer.sheets:
        worksheet = writer.sheets[sheet_name]
        for column in worksheet.columns:
            max_length = 0
            column_letter = column[0].column_letter
            for cell in column:
                try:
                    if len(str(cell.value)) > max_length:
                        max_length = len(str(cell.value))
                except:
                    pass
            adjusted_width = (max_length + 2)
            adjusted_width = min(adjusted_width, 50)
            worksheet.column_dimensions[column_letter].width = adjusted_width

print(f"Excel 报告已保存！")

# ================= 5. 生成详细图表 =================
print("正在绘制详细图表...")

h_genus = 10
h_species = max(12, TOP_N_SPECIES_PLOT * 0.4) 

fig = plt.figure(figsize=(20, h_genus + h_species + 4))
gs = fig.add_gridspec(2, 2, height_ratios=[h_genus, h_species])

# --- 图 1: 属 - 详细柱状图 ---
ax1 = fig.add_subplot(gs[0, 0])
data_g = df_genus.head(TOP_N_GENUS_PLOT)
sns.barplot(data=data_g, x='Count', y='Name', ax=ax1, palette="viridis")
ax1.set_title(f'Top {TOP_N_GENUS_PLOT} Host Genera (Unique OriC)', fontsize=16, weight='bold')
ax1.set_xlabel('Count of Unique OriC')
ax1.set_ylabel('Genus')
for i, v in enumerate(data_g['Count']):
    ax1.text(v + max(data_g['Count'])*0.01, i, f"{int(v):,}", va='center', fontsize=10)

# --- 图 2: 属 - 环形图 (Top 10) ---
ax2 = fig.add_subplot(gs[0, 1])
top_10_g = df_genus.head(10)
other_g = df_genus.iloc[10:]['Count'].sum()
plot_data_g = pd.concat([top_10_g, pd.DataFrame({'Name': ['Other'], 'Count': [other_g]})])
colors_g = sns.color_palette("Set2", len(plot_data_g))
wedges, texts, autotexts = ax2.pie(plot_data_g['Count'], labels=plot_data_g['Name'], 
                                   autopct='%1.1f%%', startangle=90, colors=colors_g,
                                   textprops={'fontsize': 11})
centre_circle = plt.Circle((0,0),0.70,fc='white')
ax2.add_artist(centre_circle)
ax2.set_title(f'Host Genera Proportion (Top 10 + Other)', fontsize=16, weight='bold')
for autotext in autotexts:
    autotext.set_color('black')
    autotext.set_fontweight('bold')

# --- 图 3: 种 - 超详细横向柱状图 ---
ax3 = fig.add_subplot(gs[1, 0])
data_s = df_species.head(TOP_N_SPECIES_PLOT)
sns.barplot(data=data_s, x='Count', y='Name', ax=ax3, palette="magma")
ax3.set_title(f'Top {TOP_N_SPECIES_PLOT} Host Species (Unique OriC)', fontsize=16, weight='bold')
ax3.set_xlabel('Count of Unique OriC')
ax3.set_ylabel('Species')
for i, v in enumerate(data_s['Count']):
    ax3.text(v + max(data_s['Count'])*0.01, i, f"{int(v):,}", va='center', fontsize=9)

# --- 图 4: 种 - 帕累托曲线 ---
ax4 = fig.add_subplot(gs[1, 1])
data_s_pareto = df_species.head(TOP_N_SPECIES_PLOT).copy()
bars = ax4.bar(data_s_pareto.index, data_s_pareto['Count'], color='#e67e22', alpha=0.6, label='Count')
ax4_twin = ax4.twinx()
line = ax4_twin.plot(data_s_pareto.index, data_s_pareto['Cumulative_Percentage'], color='red', marker='o', linewidth=2, markersize=6, label='Cumulative %')

ax4.set_title(f'Species Distribution & Cumulative % (Pareto)', fontsize=16, weight='bold')
ax4.set_xlabel(f'Rank (Top {TOP_N_SPECIES_PLOT})')
ax4.set_ylabel('Sequence Count')
ax4_twin.set_ylabel('Cumulative Percentage (%)')
ax4.set_xticks(data_s_pareto.index)
ax4.set_xticklabels(data_s_pareto['Name'], rotation=45, ha='right', fontsize=8)
ax4.grid(axis='y', linestyle='--', alpha=0.7)
ax4.legend(loc='upper left')
ax4_twin.legend(loc='upper right')

# 标注 80% 点
for i, pct in enumerate(data_s_pareto['Cumulative_Percentage']):
    if pct >= 80 and i > 0:
        ax4_twin.annotate(f'{pct:.1f}%\n@ Rank {i+1}', xy=(i, pct), xytext=(5, 5), 
                          textcoords='offset points', color='red', fontsize=10, fontweight='bold')
        break

plt.tight_layout()
plt.savefig(output_image, dpi=300, bbox_inches='tight')
print(f"图表已保存：{output_image}")
plt.show()

print("\n=== 任务完成 ===")
print(f"注意：所有统计均基于去重后的唯一 OriC 序列，且保留了原始物种名称。")
print(f"1. Excel 统计表：{os.path.abspath(output_excel)}")
print(f"2. 可视化图表：{os.path.abspath(output_image)}")