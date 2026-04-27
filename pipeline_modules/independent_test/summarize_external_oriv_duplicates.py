from __future__ import annotations

"""Summarize duplicate oriV sequences across external parquet files."""

import argparse
from pathlib import Path

try:
    import pandas as pd
except ModuleNotFoundError:
    pd = None


ORIV_ALIASES = ["oriv_sequence", "oriV_sequence", "ori_sequence", "OriC sequence", "sequence"]
PLASMID_ALIASES = ["plasmid_id", "plasmid", "accession", "genome_id", "id"]
REP_ALIASES = ["rep_protein", "rep_seq", "rep_sequence", "protein_sequence", "Rep protein"]
HOST_ALIASES = ["host_species", "species", "host", "organism"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Count duplicate oriV sequences in external parquet files.")
    parser.add_argument(
        "--external-parquet",
        nargs="+",
        default=[
            "test-00000-of-00001.parquet",
            "train-00000-of-00001 (1).parquet",
            "validation-00000-of-00001 (1).parquet",
        ],
        help="External parquet files to inspect.",
    )
    parser.add_argument(
        "--output-dir",
        default="output-lab/independent-test/oriv-duplicate-summary",
        help="Directory for duplicate summary CSV files.",
    )
    return parser.parse_args()


def find_col(df: pd.DataFrame, aliases: list[str]) -> str | None:
    lowered = {str(col).strip().lower(): col for col in df.columns}
    for alias in aliases:
        hit = lowered.get(alias.strip().lower())
        if hit is not None:
            return hit
    return None


def col_or_empty(df: pd.DataFrame, aliases: list[str]) -> pd.Series:
    col = find_col(df, aliases)
    if col is None:
        return pd.Series([""] * len(df), index=df.index, dtype="object")
    return df[col]


def main() -> None:
    args = parse_args()
    if pd is None:
        raise RuntimeError("This script requires pandas. Install it with: python -m pip install pandas pyarrow")

    frames = []
    for path_text in args.external_parquet:
        path = Path(path_text)
        df = pd.read_parquet(path)
        oriv_col = find_col(df, ORIV_ALIASES)
        if oriv_col is None:
            raise ValueError(f"No oriV sequence column found in {path}")

        frame = pd.DataFrame(
            {
                "source_file": path.name,
                "source_row": range(len(df)),
                "oriv_sequence": df[oriv_col].fillna("").astype(str).str.strip().str.lower(),
                "plasmid_id": col_or_empty(df, PLASMID_ALIASES).fillna("").astype(str).str.strip(),
                "rep_protein": col_or_empty(df, REP_ALIASES).fillna("").astype(str).str.strip().str.upper(),
                "host_species": col_or_empty(df, HOST_ALIASES).fillna("").astype(str).str.strip(),
            }
        )
        frames.append(frame)

    all_df = pd.concat(frames, ignore_index=True)
    all_df = all_df[all_df["oriv_sequence"] != ""].copy()

    counts = all_df.groupby("oriv_sequence", dropna=False).agg(
        row_count=("oriv_sequence", "size"),
        source_file_count=("source_file", "nunique"),
        plasmid_id_count=("plasmid_id", lambda x: x[x != ""].nunique()),
        rep_protein_count=("rep_protein", lambda x: x[x != ""].nunique()),
        host_species_count=("host_species", lambda x: x[x != ""].nunique()),
        source_files=("source_file", lambda x: ";".join(sorted(set(x)))),
    )
    counts = counts.reset_index().sort_values(["row_count", "source_file_count"], ascending=[False, False])
    duplicated_orivs = counts[counts["row_count"] > 1].copy()
    duplicated_rows = all_df[all_df["oriv_sequence"].isin(set(duplicated_orivs["oriv_sequence"]))].copy()
    duplicated_rows = duplicated_rows.merge(
        duplicated_orivs[["oriv_sequence", "row_count", "source_file_count"]],
        on="oriv_sequence",
        how="left",
    ).sort_values(["row_count", "oriv_sequence", "source_file", "source_row"], ascending=[False, True, True, True])

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_dir / "oriv_counts.csv", index=False, encoding="utf-8-sig")
    duplicated_orivs.to_csv(output_dir / "duplicated_oriv_counts.csv", index=False, encoding="utf-8-sig")
    duplicated_rows.to_csv(output_dir / "duplicated_oriv_rows.csv", index=False, encoding="utf-8-sig")

    if counts.empty:
        top_rows = all_df.iloc[0:0].copy()
        top_count = 0
    else:
        top_oriv = counts.iloc[0]["oriv_sequence"]
        top_count = int(counts.iloc[0]["row_count"])
        top_rows = all_df[all_df["oriv_sequence"] == top_oriv].copy()
        top_rows.to_csv(output_dir / "top_duplicated_oriv_rows.csv", index=False, encoding="utf-8-sig")
        (output_dir / "top_duplicated_oriv_sequence.txt").write_text(str(top_oriv) + "\n", encoding="utf-8")

    total_rows = len(all_df)
    unique_orivs = all_df["oriv_sequence"].nunique()
    duplicated_oriv_count = len(duplicated_orivs)
    rows_with_duplicated_oriv = len(duplicated_rows)
    redundant_rows_if_oriv_dedup = total_rows - unique_orivs

    print(f"total_rows_with_oriv: {total_rows}")
    print(f"unique_oriv_sequences: {unique_orivs}")
    print(f"duplicated_oriv_sequences: {duplicated_oriv_count}")
    print(f"rows_with_duplicated_oriv: {rows_with_duplicated_oriv}")
    print(f"redundant_rows_if_dedup_by_oriv: {redundant_rows_if_oriv_dedup}")
    print(f"max_rows_for_one_oriv: {int(counts['row_count'].max()) if not counts.empty else 0}")
    print(f"top_duplicated_oriv_rows: {top_count}")
    print(f"top_duplicated_oriv_output: {output_dir / 'top_duplicated_oriv_rows.csv'}")
    print(f"output_dir: {output_dir}")


if __name__ == "__main__":
    main()
