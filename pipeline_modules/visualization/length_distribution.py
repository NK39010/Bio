from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Visualize OriC length distribution.")
    parser.add_argument("--input", required=True, help="Input CSV file path.")
    parser.add_argument("--output-image", required=True, help="Output image path.")
    parser.add_argument("--max-box-y", type=int, default=10000, help="Y max for boxplot zoom.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_path = Path(args.output_image)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams["font.sans-serif"] = ["SimHei", "Arial Unicode MS", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    df = pd.read_csv(args.input, low_memory=False)
    if "OriC sequence" not in df.columns:
        raise KeyError("Input CSV must contain 'OriC sequence'.")

    df = df.dropna(subset=["OriC sequence"]).copy()
    df["length"] = df["OriC sequence"].astype(str).str.len()
    if df.empty:
        raise ValueError("No non-empty 'OriC sequence' rows found for plotting.")

    fig, axes = plt.subplots(1, 2, figsize=(20, 7))
    sns.histplot(
        data=df,
        x="length",
        bins=100,
        log_scale=(True, False),
        ax=axes[0],
        color="skyblue",
        edgecolor="black",
    )
    axes[0].set_title("Distribution of OriC Length (Log Scale X-axis)")
    axes[0].set_xlabel("Length (bp) - Log Scale")
    axes[0].set_ylabel("Count")
    axes[0].axvline(df["length"].median(), color="red", linestyle="--", label=f"Median: {df['length'].median()}")
    axes[0].legend()

    sns.boxplot(data=df, y="length", ax=axes[1], color="lightgreen")
    axes[1].set_title(f"Boxplot (Zoomed: 0-{args.max_box_y}bp)")
    axes[1].set_ylabel("Length (bp)")
    axes[1].set_ylim(0, args.max_box_y)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close("all")

    print(f"Saved length distribution figure: {output_path}")
    print(df["length"].describe())


if __name__ == "__main__":
    main()

