from __future__ import annotations

"""Build unified dataset outputs (all-data CSV + train/validation parquet)."""

import argparse
import logging
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


logger = logging.getLogger(__name__)


PARQUET_COLS = [
    "rep_protein",
    "oriv_sequence",
    "host_species",
    "database_source",
    "plasmid_id",
    "split",
]


def clean_seq(seq, is_aa: bool = False) -> str:
    if pd.isna(seq):
        return ""
    seq_text = str(seq).strip()
    return seq_text.upper() if is_aa else seq_text.lower()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build dataset parquet files and one all-in-one CSV table."
    )
    parser.add_argument("--input", "--input-csv", dest="input", default="final_data.csv", help="Input CSV path.")
    parser.add_argument(
        "--output-dir",
        default="training_data",
        help="Directory for generated CSV/parquet files.",
    )
    parser.add_argument(
        "--validation-ratio",
        type=float,
        default=0.25,
        help="Fallback validation ratio if split column is absent.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed for fallback split.")
    parser.add_argument(
        "--include-mechanism-column",
        action="store_true",
        help="Include replication mechanism column in output dataset files.",
    )
    parser.add_argument(
        "--mechanism-column",
        default="replication_mechanism_term",
        help="Source mechanism column name used when --include-mechanism-column is set.",
    )
    parser.add_argument("--log-level", default="INFO", help="Logging level: DEBUG, INFO, WARNING, ERROR.")
    parser.add_argument("--log-file", default=None, help="Optional log file path.")
    return parser.parse_args()


def configure_logging(level: str, log_file: str | None = None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_path, encoding="utf-8"))

    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
        force=True,
    )


def normalize_split_series(df: pd.DataFrame, validation_ratio: float, seed: int) -> pd.Series:
    split_source = None
    for col in ("split", "cdhit_split", "split_y", "split_x"):
        if col in df.columns:
            split_source = col
            break

    if split_source is None:
        if not 0 <= validation_ratio < 1:
            raise ValueError("--validation-ratio must be >= 0 and < 1.")
        logger.warning("No split column found, using random fallback split.")
        split = pd.Series(["train"] * len(df), index=df.index, dtype="object")
        val_rows = int(round(len(df) * validation_ratio))
        if validation_ratio > 0 and len(df) > 0:
            val_rows = max(1, val_rows)
        if val_rows > 0:
            val_index = df.sample(n=val_rows, random_state=seed).index
            split.loc[val_index] = "validation"
        return split

    raw = df[split_source].fillna("").astype(str).str.strip().str.lower()
    return raw.map(lambda x: "validation" if x in {"validation", "val", "dev", "eval", "test"} else "train")


def build_standardized_table(
    df: pd.DataFrame,
    split_series: pd.Series,
    include_mechanism_column: bool,
    mechanism_column: str,
) -> pd.DataFrame:
    def col_or_empty(name: str) -> pd.Series:
        if name in df.columns:
            return df[name]
        return pd.Series([""] * len(df), index=df.index, dtype="object")

    out = pd.DataFrame(index=df.index)
    out["rep_protein"] = col_or_empty("rep_seq").apply(lambda x: clean_seq(x, is_aa=True))
    out["oriv_sequence"] = col_or_empty("OriC sequence").apply(lambda x: clean_seq(x, is_aa=False))
    out["host_species"] = col_or_empty("species").fillna("").astype(str).str.strip()
    out["database_source"] = "lab"
    out["plasmid_id"] = col_or_empty("plasmid_id").fillna("").astype(str).str.strip()
    out["split"] = split_series.fillna("train").astype(str).str.strip()
    if include_mechanism_column:
        out["replication_mechanism_term"] = (
            col_or_empty(mechanism_column).fillna("unresolved").astype(str).str.strip().replace("", "unresolved")
        )

    # Enforce the required parquet/csv column order.
    ordered_cols = PARQUET_COLS.copy()
    if include_mechanism_column:
        ordered_cols.append("replication_mechanism_term")
    return out[ordered_cols].copy()


def write_parquet(df: pd.DataFrame, path: Path) -> None:
    table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(table, path, compression="snappy")


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level, args.log_file)

    logger.info("Loading input table from %s", args.input)
    df = pd.read_csv(args.input, low_memory=False)
    logger.info("Loaded %d rows and %d columns", len(df), len(df.columns))

    split_series = normalize_split_series(df, args.validation_ratio, args.seed)
    dataset_df = build_standardized_table(
        df,
        split_series,
        include_mechanism_column=args.include_mechanism_column,
        mechanism_column=args.mechanism_column,
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # CSV dataset: one total table containing all rows.
    all_csv = output_dir / "all_data.csv"
    dataset_df.to_csv(all_csv, index=False, encoding="utf-8-sig")

    train_df = dataset_df[dataset_df["split"] == "train"].copy()
    val_df = dataset_df[dataset_df["split"] == "validation"].copy()

    train_parquet = output_dir / "train.parquet"
    val_parquet = output_dir / "validation.parquet"
    write_parquet(train_df, train_parquet)
    write_parquet(val_df, val_parquet)

    total_out = len(dataset_df)
    train_pct = (len(train_df) / total_out * 100) if total_out else 0
    val_pct = (len(val_df) / total_out * 100) if total_out else 0

    logger.info("all csv: %s", all_csv)
    logger.info("train parquet: %s", train_parquet)
    logger.info("validation parquet: %s", val_parquet)
    logger.info("output columns: %s", ", ".join(PARQUET_COLS))
    logger.info("train rows: %d (%.2f%%)", len(train_df), train_pct)
    logger.info("validation rows: %d (%.2f%%)", len(val_df), val_pct)


if __name__ == "__main__":
    main()
