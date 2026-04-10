import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="OriC length distribution visualization.")
    parser.add_argument("--input", default="merged_final_optimized.csv", help="Input CSV file.")
    parser.add_argument("--output-image", default="length_distribution.png", help="Output image path.")
    parser.add_argument("--max-box-y", type=int, default=10000, help="Y max for boxplot zoom.")
    return parser.parse_args()

def main():
    args = parse_args()
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False

    df = pd.read_csv(args.input)
    df = df.dropna(subset=['OriC sequence'])
    df['length'] = df['OriC sequence'].astype(str).str.len()

    fig, axes = plt.subplots(1, 2, figsize=(20, 7))

    sns.histplot(
        data=df,
        x='length',
        bins=100,
        log_scale=(True, False),
        ax=axes[0],
        color='skyblue',
        edgecolor='black',
    )
    axes[0].set_title('Distribution of OriC Length (Log Scale X-axis)', fontsize=14)
    axes[0].set_xlabel('Length (bp) - Log Scale', fontsize=12)
    axes[0].set_ylabel('Count', fontsize=12)
    axes[0].axvline(df['length'].median(), color='red', linestyle='--', label=f'Median: {df["length"].median()}')
    axes[0].legend()

    sns.boxplot(data=df, y='length', ax=axes[1], color='lightgreen')
    axes[1].set_title(f'Boxplot (Zoomed in: 0-{args.max_box_y}bp)', fontsize=14)
    axes[1].set_ylabel('Length (bp)', fontsize=12)
    axes[1].set_ylim(0, args.max_box_y)

    plt.tight_layout()
    plt.savefig(args.output_image, dpi=300, bbox_inches='tight')
    plt.show()

    print(f"图已保存到: {args.output_image}")
    print("📊 详细长度统计：")
    print(df['length'].describe())


if __name__ == "__main__":
    main()