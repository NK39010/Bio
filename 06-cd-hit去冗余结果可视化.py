import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Visualize CD-HIT redundancy from .clstr.")
    parser.add_argument("--clstr", default="cd_hit_process.fasta.clstr", help="Input .clstr file.")
    parser.add_argument("--output", default="cdhit_redundancy_fixed90.png", help="Output image file.")
    return parser.parse_args()

def main():
    args = parse_args()
    clstr_file = args.clstr
    output_file = args.output

    if not os.path.exists(clstr_file):
        print(f"❌ 错误：找不到文件 {clstr_file}")
        return

    print(">>> 正在解析 .clstr 文件...")
    cluster_sizes = []

    with open(clstr_file, 'r') as f:
        current_size = 0
        for line in f:
            if line.startswith('>Cluster'):
                if current_size > 0:
                    cluster_sizes.append(current_size)
                current_size = 0
            else:
                current_size += 1

        if current_size > 0:
            cluster_sizes.append(current_size)

    sizes_series = pd.Series(cluster_sizes)
    total_sequences = sizes_series.sum()
    total_clusters = len(sizes_series)

    print(f"解析完成：共 {total_clusters} 个簇，包含 {total_sequences} 条序列")

    sns.set_theme(style="whitegrid")
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 12

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    size_counts = sizes_series.value_counts().sort_index()
    sns.barplot(x=size_counts.index, y=size_counts.values, ax=ax1, color='#4c72b0')
    ax1.set_title("Distribution of Cluster Sizes", fontsize=14)
    ax1.set_xlabel("Sequences per Cluster (Redundancy Level)", fontsize=12)
    ax1.set_ylabel("Number of Clusters", fontsize=12)
    ax1.set_yscale('log')
    ax1.set_xlabel("Cluster Size (1 = Unique)")

    sorted_sizes = sizes_series.sort_values(ascending=False)
    cumulative_sequences = sorted_sizes.cumsum()
    cumulative_pct = (cumulative_sequences / total_sequences) * 100

    ax2.plot(cumulative_pct.index, cumulative_pct.values, color='#55a868', linewidth=2)
    ax2.fill_between(cumulative_pct.index, cumulative_pct.values, color='#55a868', alpha=0.2)

    ax2.set_title("Cumulative Redundancy Impact", fontsize=14)
    ax2.set_xlabel("Top N Largest Clusters", fontsize=12)
    ax2.set_ylabel("Cumulative Data Coverage (%)", fontsize=12)

    top_10_pct_idx = int(len(cumulative_pct) * 0.1)
    if len(cumulative_pct) == 0:
        print("⚠️ .clstr 解析结果为空。")
        return

    top_10_pct_idx = max(0, min(top_10_pct_idx, len(cumulative_pct) - 1))
    y_val = cumulative_pct.iloc[top_10_pct_idx]
    ax2.scatter([top_10_pct_idx], [y_val], color='red', s=50, zorder=5)
    ax2.annotate(
        f'Top 10% Clusters\ncontain {y_val:.1f}% data',
        xy=(top_10_pct_idx, y_val),
        xytext=(top_10_pct_idx + 50, y_val - 20),
        arrowprops=dict(arrowstyle='->', color='gray'),
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"✅ 图表已保存为：{output_file}")
    plt.show()


if __name__ == "__main__":
    main()