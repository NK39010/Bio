from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd


DEFAULT_OUTPUT_DIR = "id_relationship_report"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize and visualize plasmid_id, ori_id, and rep_id relationships."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input CSV or Parquet file.",
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Directory to store all relationship outputs.",
    )
    parser.add_argument("--plasmid-col", default="plasmid_id", help="Plasmid ID column name.")
    parser.add_argument("--ori-col", default="ori_id", help="Ori ID column name.")
    parser.add_argument("--rep-col", default="rep_id", help="Rep ID column name.")
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top entities to show in plots and text summaries.",
    )
    return parser.parse_args()


def load_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix == ".csv":
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"Unsupported input format: {path.suffix}. Use .csv or .parquet")


def normalize_ids(df: pd.DataFrame, cols: Iterable[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col not in out.columns:
            raise KeyError(f"Missing required column: {col}. Available columns: {list(out.columns)}")
        out[col] = out[col].fillna("").astype(str).str.strip()
        out = out[out[col] != ""]
    return out


def summarize_groups(df: pd.DataFrame, group_col: str, target_col: str) -> pd.DataFrame:
    grouped = (
        df.groupby(group_col, dropna=False)[target_col]
        .agg(
            n_unique="nunique",
            n_rows="size",
            examples=lambda s: ";".join(pd.Series(s.dropna().astype(str).unique()).head(8)),
        )
        .reset_index()
        .sort_values(["n_unique", "n_rows", group_col], ascending=[False, False, True])
    )
    return grouped


def format_markdown_table(df: pd.DataFrame, max_rows: int = 20) -> str:
    sample = df.head(max_rows).copy()
    if sample.empty:
        return "_No rows found._"

    sample = sample.fillna("")
    headers = list(sample.columns)
    rows = [headers]
    rows.extend(sample.astype(str).values.tolist())

    widths = [max(len(str(row[i])) for row in rows) for i in range(len(headers))]

    def format_row(row: list[str]) -> str:
        cells = [str(cell).ljust(widths[idx]) for idx, cell in enumerate(row)]
        return "| " + " | ".join(cells) + " |"

    lines = [format_row(headers)]
    lines.append("| " + " | ".join("-" * w for w in widths) + " |")
    for row in rows[1:]:
        lines.append(format_row(row))
    return "\n".join(lines)


def build_pair_counts(df: pd.DataFrame, left_col: str, right_col: str) -> pd.DataFrame:
    pair_counts = (
        df.groupby([left_col, right_col], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["count", left_col, right_col], ascending=[False, True, True])
    )
    return pair_counts


def build_triplet_counts(df: pd.DataFrame, plasmid_col: str, ori_col: str, rep_col: str) -> pd.DataFrame:
    triplets = (
        df.groupby([plasmid_col, ori_col, rep_col], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["count", plasmid_col, ori_col, rep_col], ascending=[False, True, True, True])
    )
    return triplets


def plot_top_relationships(summary_df: pd.DataFrame, entity_col: str, value_col: str, title: str, output_path: Path) -> None:
    plot_df = summary_df.head(20).sort_values(value_col, ascending=True)
    if plot_df.empty:
        return

    plt.figure(figsize=(12, max(6, 0.35 * len(plot_df) + 2)))
    plt.barh(plot_df[entity_col], plot_df[value_col], color="#2E86AB")
    plt.title(title)
    plt.xlabel(value_col)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def write_markdown_report(
    output_dir: Path,
    source_path: Path,
    total_rows: int,
    unique_plasmids: int,
    unique_oris: int,
    unique_reps: int,
    multi_plasmids: int,
    multi_oris: int,
    multi_reps: int,
    plasmid_summary: pd.DataFrame,
    ori_summary: pd.DataFrame,
    rep_summary: pd.DataFrame,
    triplets: pd.DataFrame,
    top_n: int,
) -> None:
    lines = [
        "# ID Relationship Report",
        "",
        f"- Source: `{source_path}`",
        f"- Rows: `{total_rows}`",
        f"- Unique plasmids: `{unique_plasmids}`",
        f"- Unique ori IDs: `{unique_oris}`",
        f"- Unique rep IDs: `{unique_reps}`",
        f"- Unique plasmid-ori-rep triplets: `{len(triplets)}`",
        f"- Plasmids with multiple ori IDs: `{multi_plasmids}`",
        f"- Ori IDs with multiple rep IDs: `{multi_oris}`",
        f"- Rep IDs with multiple ori IDs: `{multi_reps}`",
        "",
        "## Most connected plasmids",
        "",
        format_markdown_table(plasmid_summary, max_rows=top_n),
        "",
        "## Most connected ori IDs",
        "",
        format_markdown_table(ori_summary, max_rows=top_n),
        "",
        "## Most connected rep IDs",
        "",
        format_markdown_table(rep_summary, max_rows=top_n),
        "",
        "## Top triplets",
        "",
        format_markdown_table(triplets, max_rows=top_n),
        "",
    ]
    (output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_table(input_path)
    df = normalize_ids(df, [args.plasmid_col, args.ori_col, args.rep_col])

    df = df.copy()
    df["triple_id"] = (
        df[args.plasmid_col].astype(str)
        + " || "
        + df[args.ori_col].astype(str)
        + " || "
        + df[args.rep_col].astype(str)
    )

    triplets = build_triplet_counts(df, args.plasmid_col, args.ori_col, args.rep_col)
    plasmid_ori = build_pair_counts(df, args.plasmid_col, args.ori_col)
    plasmid_rep = build_pair_counts(df, args.plasmid_col, args.rep_col)
    ori_rep = build_pair_counts(df, args.ori_col, args.rep_col)

    plasmid_summary = summarize_groups(df, args.plasmid_col, args.ori_col).rename(
        columns={args.plasmid_col: "plasmid_id", "n_unique": "n_ori", "n_rows": "n_rows", "examples": "ori_examples"}
    )
    ori_summary = summarize_groups(df, args.ori_col, args.rep_col).rename(
        columns={args.ori_col: "ori_id", "n_unique": "n_rep", "n_rows": "n_rows", "examples": "rep_examples"}
    )
    rep_summary = summarize_groups(df, args.rep_col, args.ori_col).rename(
        columns={args.rep_col: "rep_id", "n_unique": "n_ori", "n_rows": "n_rows", "examples": "ori_examples"}
    )

    multi_plasmids = plasmid_summary[plasmid_summary["n_ori"] > 1].copy()
    multi_oris = ori_summary[ori_summary["n_rep"] > 1].copy()
    multi_reps = rep_summary[rep_summary["n_ori"] > 1].copy()

    plasmid_summary.to_csv(output_dir / "plasmid_summary.csv", index=False, encoding="utf-8-sig")
    ori_summary.to_csv(output_dir / "ori_summary.csv", index=False, encoding="utf-8-sig")
    rep_summary.to_csv(output_dir / "rep_summary.csv", index=False, encoding="utf-8-sig")
    triplets.to_csv(output_dir / "triplet_counts.csv", index=False, encoding="utf-8-sig")
    plasmid_ori.to_csv(output_dir / "plasmid_ori_counts.csv", index=False, encoding="utf-8-sig")
    plasmid_rep.to_csv(output_dir / "plasmid_rep_counts.csv", index=False, encoding="utf-8-sig")
    ori_rep.to_csv(output_dir / "ori_rep_counts.csv", index=False, encoding="utf-8-sig")
    multi_plasmids.to_csv(output_dir / "plasmids_with_multiple_oris.csv", index=False, encoding="utf-8-sig")
    multi_oris.to_csv(output_dir / "oris_with_multiple_reps.csv", index=False, encoding="utf-8-sig")
    multi_reps.to_csv(output_dir / "reps_with_multiple_oris.csv", index=False, encoding="utf-8-sig")

    plot_top_relationships(
        plasmid_summary,
        entity_col="plasmid_id",
        value_col="n_ori",
        title="Top plasmids by number of unique ori IDs",
        output_path=output_dir / "top_plasmid_ori_counts.png",
    )
    plot_top_relationships(
        ori_summary,
        entity_col="ori_id",
        value_col="n_rep",
        title="Top ori IDs by number of unique rep IDs",
        output_path=output_dir / "top_ori_rep_counts.png",
    )
    plot_top_relationships(
        rep_summary,
        entity_col="rep_id",
        value_col="n_ori",
        title="Top rep IDs by number of unique ori IDs",
        output_path=output_dir / "top_rep_ori_counts.png",
    )

    write_markdown_report(
        output_dir=output_dir,
        source_path=input_path,
        total_rows=len(df),
        unique_plasmids=df[args.plasmid_col].nunique(),
        unique_oris=df[args.ori_col].nunique(),
        unique_reps=df[args.rep_col].nunique(),
        multi_plasmids=len(multi_plasmids),
        multi_oris=len(multi_oris),
        multi_reps=len(multi_reps),
        plasmid_summary=plasmid_summary,
        ori_summary=ori_summary,
        rep_summary=rep_summary,
        triplets=triplets,
        top_n=args.top_n,
    )

    print(f"Report written to: {output_dir.resolve()}")
    print("Key files:")
    print(f"  - {output_dir / 'report.md'}")
    print(f"  - {output_dir / 'plasmid_summary.csv'}")
    print(f"  - {output_dir / 'ori_summary.csv'}")
    print(f"  - {output_dir / 'rep_summary.csv'}")
    print(f"  - {output_dir / 'triplet_counts.csv'}")
    print(f"  - {output_dir / 'plasmids_with_multiple_oris.csv'}")
    print(f"  - {output_dir / 'oris_with_multiple_reps.csv'}")
    print(f"  - {output_dir / 'reps_with_multiple_oris.csv'}")


if __name__ == "__main__":
    main()
