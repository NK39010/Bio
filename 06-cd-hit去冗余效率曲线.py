import os
import matplotlib.pyplot as plt
import numpy as np

# --- 1. 配置区域 ---
# 请确保文件名与你截图中的完全一致
files_to_compare = {
    "95% Identity": "cd_hit_process95.fasta.clstr",
    "90% Identity": "cd_hit_process.fasta.clstr",  # 假设这是默认的90%
    "85% Identity": "cd_hit_process85.fasta.clstr"
}

output_image = "redundancy_comparison.png"

# --- 2. 数据解析函数 ---
def parse_clstr_stats(filepath):
    """
    解析 .clstr 文件，返回簇大小的列表
    """
    if not os.path.exists(filepath):
        print(f"⚠️ 警告: 文件未找到 - {filepath}")
        return None

    cluster_sizes = []
    with open(filepath, 'r', encoding='utf-8') as f:
        current_size = 0
        for line in f:
            if line.startswith('>Cluster'):
                if current_size > 0:
                    cluster_sizes.append(current_size)
                current_size = 0
            else:
                current_size += 1
        # 别忘了最后一个簇
        if current_size > 0:
            cluster_sizes.append(current_size)
    
    return cluster_sizes

# --- 3. 主逻辑 ---
plt.figure(figsize=(10, 6))

total_sequences_all = 0 # 用于验证总数是否一致

for label, filename in files_to_compare.items():
    sizes = parse_clstr_stats(filename)
    
    if sizes is None:
        continue
        
    # 排序：从大到小，这是画累积曲线的关键
    sizes.sort(reverse=True)
    
    # 计算总序列数（所有簇的大小之和）
    total_seqs = sum(sizes)
    total_sequences_all = total_seqs
    
    # 计算累积和
    cum_sum = np.cumsum(sizes)
    
    # X轴：簇的数量 (1, 2, 3... N)
    x = range(1, len(sizes) + 1)
    
    # Y轴：累积覆盖百分比 (0% - 100%)
    y = [ (val / total_seqs) * 100 for val in cum_sum ]
    
    # 绘图
    plt.plot(x, y, label=label, linewidth=2)

# --- 4. 图表美化 ---
plt.title("Comparison of Redundancy Removal Efficiency", fontsize=14)
plt.xlabel("Number of Top Clusters (Sorted by Size)", fontsize=12)
plt.ylabel("Cumulative Data Coverage (%)", fontsize=12)
plt.legend(title="CD-HIT Threshold", loc='lower right')
plt.grid(True, linestyle='--', alpha=0.6)

# 添加说明文字
plt.text(0.05, 0.05, f"Total Sequences: {total_sequences_all}", transform=plt.gca().transAxes, fontsize=10, color='gray')

# 保存并显示
plt.tight_layout()
plt.savefig(output_image, dpi=300)
print(f"✅ 图表已生成: {output_image}")
plt.show()