import argparse
import os
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Species distribution report and visualization.")
    parser.add_argument("--input", default="merged_final_optimized.csv", help="Input CSV file.")
    parser.add_argument("--output-excel", default="csv_species_statistics.xlsx", help="Output Excel file.")
    parser.add_argument("--output-image", default="csv_species_analysis.png", help="Output image file.")
    parser.add_argument("--top-genus", type=int, default=30, help="Top N genera for bar plot.")
    parser.add_argument("--top-species", type=int, default=50, help="Top N species for plots.")
    return parser.parse_args()


def main():
    args = parse_args()

    plt.rcParams["font.sans-serif"] = ["SimHei", "Arial Unicode MS", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    print(f"正在读取 CSV 文件：{args.input} ...")
    try:
        df_raw = pd.read_csv(args.input)
        print(f"原始数据行数: {len(df_raw)}")
        print("不再按 ori_id 去重，直接基于原始行统计物种分布。")

        if "species" not in df_raw.columns:
            print(f"错误：未找到 'species' 列。现有列：{df_raw.columns.tolist()}")
            return

        df = df_raw.copy()
        raw_species_list = df["species"].dropna().tolist()
    except FileNotFoundError:
        print(f"错误：找不到文件 '{args.input}'")
        return

    print("正在统计物种（保留原始名称）...")
    species_counter = Counter(raw_species_list)
    total_count = len(raw_species_list)

    genus_counter = Counter()
    for species_name, count in species_counter.items():
        genus = str(species_name).split()[0]
        genus_counter[genus] += count

    print(f"统计完成。共识别出 {len(genus_counter)} 个属，{len(species_counter)} 个种。")

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

    excel_path = Path(args.output_excel)
    print(f"正在生成 Excel 报告：{args.output_excel} ...")
    try:
        with pd.ExcelWriter(args.output_excel, engine="openpyxl") as writer:
            df_summary.to_excel(writer, sheet_name="Summary", index=False)
            df_genus.to_excel(writer, sheet_name="Genus", index=False)
            df_species.to_excel(writer, sheet_name="Species", index=False)
            df_species.head(50).to_excel(writer, sheet_name="Top50", index=False)

            for sheet_name in writer.sheets:
                worksheet = writer.sheets[sheet_name]
                for column in worksheet.columns:
                    max_length = 0
                    column_letter = column[0].column_letter
                    for cell in column:
                        try:
                            max_length = max(max_length, len(str(cell.value)))
                        except Exception:
                            pass
                    worksheet.column_dimensions[column_letter].width = min(max_length + 2, 50)
        print("Excel 报告已保存！")
    except Exception as exc:
        fallback_prefix = excel_path.with_suffix("")
        df_summary.to_csv(fallback_prefix.as_posix() + "_summary.csv", index=False, encoding="utf-8-sig")
        df_genus.to_csv(fallback_prefix.as_posix() + "_genus.csv", index=False, encoding="utf-8-sig")
        df_species.to_csv(fallback_prefix.as_posix() + "_species.csv", index=False, encoding="utf-8-sig")
        df_species.head(50).to_csv(fallback_prefix.as_posix() + "_top50.csv", index=False, encoding="utf-8-sig")
        print(f"Excel 写入失败，已改为导出 CSV：{exc}")

    print("正在绘制详细图表...")
    top_genus_n = args.top_genus
    top_species_n = args.top_species

    fig = plt.figure(figsize=(20, 18))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.2])

    ax1 = fig.add_subplot(gs[0, 0])
    data_g = df_genus.head(top_genus_n)
    ax1.barh(data_g["Name"][::-1], data_g["Count"][::-1], color=plt.cm.viridis(
        [i / max(len(data_g) - 1, 1) for i in range(len(data_g) - 1, -1, -1)]
    ))
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
    ax3.barh(data_s["Name"][::-1], data_s["Count"][::-1], color=plt.cm.magma(
        [i / max(len(data_s) - 1, 1) for i in range(len(data_s) - 1, -1, -1)]
    ))
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
    ax4.set_xticks(data_s_pareto.index)
    ax4.set_xticklabels(data_s_pareto["Name"], rotation=45, ha="right", fontsize=8)

    plt.tight_layout()
    plt.savefig(args.output_image, dpi=300, bbox_inches="tight")
    print(f"图表已保存：{args.output_image}")
    plt.close("all")

    print("\n=== 任务完成 ===")
    print("注意：所有统计均基于原始行，不再按 ori_id 去重。")
    print(f"1. Excel 统计表：{os.path.abspath(args.output_excel)}")
    print(f"2. 可视化图表：{os.path.abspath(args.output_image)}")


if __name__ == "__main__":
    main()
