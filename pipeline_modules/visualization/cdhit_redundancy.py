from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Visualize CD-HIT redundancy from .clstr.")
    parser.add_argument("--clstr", required=True, help="Input .clstr file.")
    parser.add_argument("--output", required=True, help="Output image file.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    clstr_path = Path(args.clstr)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not clstr_path.exists():
        raise FileNotFoundError(f"CLSTR file not found: {clstr_path}")

    cluster_sizes: list[int] = []
    current_size = 0
    with clstr_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">Cluster"):
                if current_size > 0:
                    cluster_sizes.append(current_size)
                current_size = 0
            else:
                current_size += 1
    if current_size > 0:
        cluster_sizes.append(current_size)

    if not cluster_sizes:
        raise ValueError("No clusters parsed from the .clstr file.")

    sizes = pd.Series(cluster_sizes)
    total_sequences = int(sizes.sum())
    total_clusters = int(len(sizes))

    sns.set_theme(style="whitegrid")
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 12
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    size_counts = sizes.value_counts().sort_index()
    sns.barplot(x=size_counts.index, y=size_counts.values, ax=ax1, color="#4c72b0")
    ax1.set_title("Distribution of Cluster Sizes")
    ax1.set_xlabel("Cluster Size (1 = Unique)")
    ax1.set_ylabel("Number of Clusters")
    ax1.set_yscale("log")

    sorted_sizes = sizes.sort_values(ascending=False)
    cumulative_sequences = sorted_sizes.cumsum()
    cumulative_pct = (cumulative_sequences / total_sequences) * 100
    ax2.plot(cumulative_pct.index, cumulative_pct.values, color="#55a868", linewidth=2)
    ax2.fill_between(cumulative_pct.index, cumulative_pct.values, color="#55a868", alpha=0.2)
    ax2.set_title("Cumulative Redundancy Impact")
    ax2.set_xlabel("Top N Largest Clusters")
    ax2.set_ylabel("Cumulative Data Coverage (%)")

    marker_idx = max(0, min(int(len(cumulative_pct) * 0.1), len(cumulative_pct) - 1))
    marker_val = cumulative_pct.iloc[marker_idx]
    ax2.scatter([marker_idx], [marker_val], color="red", s=50, zorder=5)
    ax2.annotate(
        f"Top 10% Clusters\ncontain {marker_val:.1f}% data",
        xy=(marker_idx, marker_val),
        xytext=(marker_idx + 10, max(0.0, marker_val - 20)),
        arrowprops=dict(arrowstyle="->", color="gray"),
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close("all")

    print(f"Saved CD-HIT redundancy figure: {out_path}")
    print(f"Clusters: {total_clusters}, sequences: {total_sequences}")


if __name__ == "__main__":
    main()

