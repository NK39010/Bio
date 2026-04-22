from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Species distribution report and visualization.")
    parser.add_argument("--input", required=True, help="Input CSV file.")
    parser.add_argument("--output-excel", required=True, help="Output Excel file.")
    parser.add_argument("--output-image", required=True, help="Output image file.")
    parser.add_argument("--top-genus", type=int, default=30, help="Top N genera for bar plot.")
    parser.add_argument("--top-species", type=int, default=50, help="Top N species for plots.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    excel_path = Path(args.output_excel)
    image_path = Path(args.output_image)
    excel_path.parent.mkdir(parents=True, exist_ok=True)
    image_path.parent.mkdir(parents=True, exist_ok=True)

    plt.rcParams["font.sans-serif"] = ["SimHei", "Arial Unicode MS", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    df = pd.read_csv(args.input, low_memory=False)
    if "species" not in df.columns:
        raise KeyError("Input CSV must contain 'species'.")

    raw_species = df["species"].dropna().astype(str).str.strip()
    raw_species = raw_species[raw_species != ""]
    if raw_species.empty:
        raise ValueError("No non-empty species rows found for plotting.")

    species_counter = Counter(raw_species.tolist())
    total_count = len(raw_species)

    genus_counter = Counter()
    for species_name, count in species_counter.items():
        genus = str(species_name).split()[0]
        genus_counter[genus] += count

    df_genus = pd.DataFrame.from_dict(genus_counter, orient="index", columns=["Count"])
    df_genus = df_genus.sort_values(by="Count", ascending=False).reset_index()
    df_genus.columns = ["Name", "Count"]
    df_genus["Percentage"] = (df_genus["Count"] / total_count * 100).round(4)
    df_genus["Cumulative_Count"] = df_genus["Count"].cumsum()
    df_genus["Cumulative_Percentage"] = (df_genus["Cumulative_Count"] / total_count * 100).round(4)
    df_genus["Rank"] = df_genus.index + 1

    df_species = pd.DataFrame.from_dict(species_counter, orient="index", columns=["Count"])
    df_species = df_species.sort_values(by="Count", ascending=False).reset_index()
    df_species.columns = ["Name", "Count"]
    df_species["Percentage"] = (df_species["Count"] / total_count * 100).round(4)
    df_species["Cumulative_Count"] = df_species["Count"].cumsum()
    df_species["Cumulative_Percentage"] = (df_species["Cumulative_Count"] / total_count * 100).round(4)
    df_species["Rank"] = df_species.index + 1

    summary_data = {
        "Metric": [
            "Total Rows",
            "Unique Genera",
            "Unique Species",
            "Top 1 Genus Name",
            "Top 1 Genus Count",
            "Top 1 Genus %",
            "Top 10 Genera %",
            "Top 1 Species Name",
            "Top 1 Species Count",
            "Top 1 Species %",
            "Top 10 Species %",
            "Top 50 Species %",
        ],
        "Value": [
            total_count,
            len(df_genus),
            len(df_species),
            df_genus.iloc[0]["Name"] if len(df_genus) > 0 else "N/A",
            df_genus.iloc[0]["Count"] if len(df_genus) > 0 else 0,
            f"{df_genus.iloc[0]['Percentage']:.2f}%" if len(df_genus) > 0 else "0%",
            f"{df_genus.head(10)['Percentage'].sum():.2f}%" if len(df_genus) > 0 else "0%",
            df_species.iloc[0]["Name"] if len(df_species) > 0 else "N/A",
            df_species.iloc[0]["Count"] if len(df_species) > 0 else 0,
            f"{df_species.iloc[0]['Percentage']:.2f}%" if len(df_species) > 0 else "0%",
            f"{df_species.head(10)['Percentage'].sum():.2f}%" if len(df_species) > 0 else "0%",
            f"{df_species.head(50)['Percentage'].sum():.2f}%" if len(df_species) > 0 else "0%",
        ],
    }
    df_summary = pd.DataFrame(summary_data)

    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        df_summary.to_excel(writer, sheet_name="Summary", index=False)
        df_genus.to_excel(writer, sheet_name="Genus", index=False)
        df_species.to_excel(writer, sheet_name="Species", index=False)
        df_species.head(50).to_excel(writer, sheet_name="Top50", index=False)

    top_genus_n = args.top_genus
    top_species_n = args.top_species
    fig = plt.figure(figsize=(20, 18))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.2])

    ax1 = fig.add_subplot(gs[0, 0])
    data_g = df_genus.head(top_genus_n)
    ax1.barh(data_g["Name"][::-1], data_g["Count"][::-1], color="#4c72b0")
    ax1.set_title(f"Top {top_genus_n} Genera (Raw Rows)")
    ax1.set_xlabel("Count")
    ax1.set_ylabel("Genus")

    ax2 = fig.add_subplot(gs[0, 1])
    top_10_g = df_genus.head(10)
    other_g = df_genus.iloc[10:]["Count"].sum()
    plot_data_g = pd.concat([top_10_g, pd.DataFrame({"Name": ["Other"], "Count": [other_g]})])
    ax2.pie(plot_data_g["Count"], labels=plot_data_g["Name"], autopct="%1.1f%%", startangle=90)
    ax2.set_title("Host Genera Proportion (Top 10 + Other)")

    ax3 = fig.add_subplot(gs[1, 0])
    data_s = df_species.head(top_species_n)
    ax3.barh(data_s["Name"][::-1], data_s["Count"][::-1], color="#e67e22")
    ax3.set_title(f"Top {top_species_n} Species (Raw Rows)")
    ax3.set_xlabel("Count")
    ax3.set_ylabel("Species")

    ax4 = fig.add_subplot(gs[1, 1])
    data_s_pareto = df_species.head(top_species_n).copy()
    ax4.bar(data_s_pareto.index, data_s_pareto["Count"], color="#e67e22", alpha=0.6)
    ax4_twin = ax4.twinx()
    ax4_twin.plot(data_s_pareto.index, data_s_pareto["Cumulative_Percentage"], color="red", marker="o")
    ax4.set_title("Species Distribution & Cumulative %")
    ax4.set_xlabel(f"Rank (Top {top_species_n})")
    ax4.set_ylabel("Count")
    ax4_twin.set_ylabel("Cumulative Percentage (%)")

    plt.tight_layout()
    plt.savefig(image_path, dpi=300, bbox_inches="tight")
    plt.close("all")

    print(f"Saved species report: {excel_path}")
    print(f"Saved species figure: {image_path}")


if __name__ == "__main__":
    main()

