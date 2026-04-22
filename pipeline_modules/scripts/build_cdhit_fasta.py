from __future__ import annotations

"""Build CD-HIT FASTA input from the annotated CSV table."""

import argparse
import logging
from pathlib import Path

import pandas as pd


logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build FASTA input for CD-HIT from annotated CSV.")
    parser.add_argument("--input", "--input-csv", dest="input", default="merged_final_optimized_replication_annotated.csv", help="Input CSV path.")
    parser.add_argument("--output", "--output-fasta", dest="output", default="ori.fasta", help="Output FASTA path.")
    parser.add_argument("--min-len", type=int, default=0, help="Minimum sequence length.")
    parser.add_argument("--max-len", type=int, default=20000, help="Maximum sequence length.")
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


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level, args.log_file)

    if args.min_len < 0:
        raise ValueError("--min-len must be >= 0.")
    if args.max_len < args.min_len:
        raise ValueError("--max-len must be >= --min-len.")

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Loading annotated CSV from %s", input_path)
    df = pd.read_csv(input_path, low_memory=False)
    logger.info("Loaded %d rows", len(df))

    required_cols = {"OriC sequence", "ori_id"}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise KeyError(f"Input CSV is missing required columns: {', '.join(sorted(missing_cols))}")

    before_dropna = len(df)
    df = df.dropna(subset=["OriC sequence"])
    logger.info("Rows with OriC sequence: %d (dropped %d missing sequences)", len(df), before_dropna - len(df))

    logger.info("Filtering sequence length range: %d-%d bp", args.min_len, args.max_len)
    df["seq_len"] = df["OriC sequence"].astype(str).str.len()
    df_filtered = df[(df["seq_len"] >= args.min_len) & (df["seq_len"] <= args.max_len)].copy()
    logger.info("Rows after length filter: %d", len(df_filtered))

    before_dedup = len(df_filtered)
    df_filtered = df_filtered.drop_duplicates(subset=["OriC sequence", "ori_id"])
    logger.info("Rows after deduplication: %d (removed %d duplicates)", len(df_filtered), before_dedup - len(df_filtered))

    logger.info("Writing FASTA to %s", output_path)
    with output_path.open("w", encoding="utf-8") as f:
        for i, row in df_filtered.iterrows():
            ori_id = str(row["ori_id"]).strip()
            seq = str(row["OriC sequence"]).strip()
            fasta_id = f"{ori_id}|idx={i}"
            f.write(f">{fasta_id}\n{seq}\n")

    logger.info("Done. FASTA records written: %d", len(df_filtered))


if __name__ == "__main__":
    main()
