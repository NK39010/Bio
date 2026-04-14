from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ID_CANDIDATES = ["plasmid_id", "plsmid_id", "plasmidid", "plasmid"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare plasmid IDs between two parquet or csv files."
    )
    parser.add_argument("--left", required=True, help="First input file (.parquet or .csv).")
    parser.add_argument("--right", required=True, help="Second input file (.parquet or .csv).")
    parser.add_argument(
        "--left-col",
        default=None,
        help="Plasmid ID column in the left file. If omitted, auto-detect.",
    )
    parser.add_argument(
        "--right-col",
        default=None,
        help="Plasmid ID column in the right file. If omitted, auto-detect.",
    )
    parser.add_argument(
        "--output-dir",
        default="relationships/compare_output",
        help="Directory to store comparison results.",
    )
    return parser.parse_args()


def load_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix == ".csv":
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"Unsupported file type: {path.suffix}")


def detect_column(df: pd.DataFrame, preferred: str | None) -> str:
    if preferred:
        if preferred not in df.columns:
            raise KeyError(f"Column '{preferred}' not found. Available: {list(df.columns)}")
        return preferred

    lowered = {col.lower(): col for col in df.columns}
    for candidate in ID_CANDIDATES:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]

    raise KeyError(
        "Could not auto-detect plasmid ID column. "
        f"Tried: {ID_CANDIDATES}. Available columns: {list(df.columns)}"
    )


def normalize_ids(df: pd.DataFrame, col: str) -> pd.Series:
    return df[col].fillna("").astype(str).str.strip()


def main() -> None:
    args = parse_args()
    left_path = Path(args.left)
    right_path = Path(args.right)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    left_df = load_table(left_path)
    right_df = load_table(right_path)

    left_col = detect_column(left_df, args.left_col)
    right_col = detect_column(right_df, args.right_col)

    left_ids = normalize_ids(left_df, left_col)
    right_ids = normalize_ids(right_df, right_col)

    left_unique = set(left_ids[left_ids != ""])
    right_unique = set(right_ids[right_ids != ""])
    common_ids = sorted(left_unique & right_unique)

    left_common_rows = left_ids.isin(common_ids).sum()
    right_common_rows = right_ids.isin(common_ids).sum()

    left_common_df = left_df[left_ids.isin(common_ids)].copy()
    right_common_df = right_df[right_ids.isin(common_ids)].copy()

    left_common_df["__plasmid_id__"] = left_ids[left_ids.isin(common_ids)].values
    right_common_df["__plasmid_id__"] = right_ids[right_ids.isin(common_ids)].values

    left_counts = (
        left_common_df.groupby("__plasmid_id__", dropna=False)
        .size()
        .reset_index(name="left_row_count")
    )
    right_counts = (
        right_common_df.groupby("__plasmid_id__", dropna=False)
        .size()
        .reset_index(name="right_row_count")
    )

    comparison = left_counts.merge(right_counts, on="__plasmid_id__", how="outer").fillna(0)
    comparison["left_row_count"] = comparison["left_row_count"].astype(int)
    comparison["right_row_count"] = comparison["right_row_count"].astype(int)
    comparison["count_diff"] = comparison["left_row_count"] - comparison["right_row_count"]
    comparison = comparison.sort_values(
        ["left_row_count", "right_row_count", "__plasmid_id__"],
        ascending=[False, False, True],
    )

    unique_summary = pd.DataFrame(
        [
            {
                "file": str(left_path),
                "column": left_col,
                "rows": len(left_df),
                "unique_plasmid_ids": len(left_unique),
            },
            {
                "file": str(right_path),
                "column": right_col,
                "rows": len(right_df),
                "unique_plasmid_ids": len(right_unique),
            },
        ]
    )

    summary = pd.DataFrame(
        [
            {"metric": "left_rows", "value": len(left_df)},
            {"metric": "right_rows", "value": len(right_df)},
            {"metric": "left_unique_plasmid_ids", "value": len(left_unique)},
            {"metric": "right_unique_plasmid_ids", "value": len(right_unique)},
            {"metric": "common_unique_plasmid_ids", "value": len(common_ids)},
            {"metric": "left_rows_with_common_ids", "value": left_common_rows},
            {"metric": "right_rows_with_common_ids", "value": right_common_rows},
        ]
    )

    summary.to_csv(output_dir / "summary.csv", index=False, encoding="utf-8-sig")
    unique_summary.to_csv(output_dir / "file_summary.csv", index=False, encoding="utf-8-sig")
    comparison.to_csv(output_dir / "common_plasmid_counts.csv", index=False, encoding="utf-8-sig")

    with open(output_dir / "summary.txt", "w", encoding="utf-8") as f:
        f.write(f"Left file:  {left_path}\n")
        f.write(f"Right file: {right_path}\n")
        f.write(f"Left column:  {left_col}\n")
        f.write(f"Right column: {right_col}\n")
        f.write(f"Left rows: {len(left_df)}\n")
        f.write(f"Right rows: {len(right_df)}\n")
        f.write(f"Left unique plasmid IDs: {len(left_unique)}\n")
        f.write(f"Right unique plasmid IDs: {len(right_unique)}\n")
        f.write(f"Common unique plasmid IDs: {len(common_ids)}\n")
        f.write(f"Left rows with common IDs: {left_common_rows}\n")
        f.write(f"Right rows with common IDs: {right_common_rows}\n")

    print(f"Comparison finished: {output_dir.resolve()}")
    print(f"Common unique plasmid IDs: {len(common_ids)}")
    print(f"Left rows with common IDs: {left_common_rows}")
    print(f"Right rows with common IDs: {right_common_rows}")
    print(f"Outputs:")
    print(f"  - {output_dir / 'summary.txt'}")
    print(f"  - {output_dir / 'summary.csv'}")
    print(f"  - {output_dir / 'file_summary.csv'}")
    print(f"  - {output_dir / 'common_plasmid_counts.csv'}")


if __name__ == "__main__":
    main()
