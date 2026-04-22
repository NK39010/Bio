from __future__ import annotations

"""Restore CD-HIT cluster labels back onto the annotated CSV table."""

import argparse
import logging
from pathlib import Path

import pandas as pd


logger = logging.getLogger(__name__)


def sort_cluster_labels(values: pd.Series) -> pd.Series:
    text = values.fillna("").astype(str)
    numeric = pd.to_numeric(text, errors="coerce")
    return numeric.where(numeric.notna(), float("inf"))


def sort_cluster_frame(df: pd.DataFrame, column: str = "cdhit_cluster") -> pd.DataFrame:
    df = df.copy()
    df["_cluster_sort_key"] = sort_cluster_labels(df[column])
    df = df.sort_values(by=["_cluster_sort_key", column]).drop(columns=["_cluster_sort_key"])
    return df.reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge CD-HIT cluster labels back to the annotated CSV without dropping members."
    )
    parser.add_argument("--fasta", "--input-fasta", dest="fasta", default="cd_hit_process.fasta", help="CD-HIT FASTA path.")
    parser.add_argument(
        "--clstr",
        default=None,
        help="CD-HIT cluster file. Default: same path as FASTA with .clstr appended.",
    )
    parser.add_argument(
        "--csv",
        "--input-csv",
        dest="csv",
        default="merged_final_optimized_replication_annotated.csv",
        help="Source annotated CSV path.",
    )
    parser.add_argument("--output", "--output-csv", dest="output", default="final_data.csv", help="Output CSV path.")
    parser.add_argument("--id-column", default="ori_id", help="ID column in source CSV.")
    parser.add_argument(
        "--dedup-mode",
        choices=["keep-all", "representative-only"],
        default="keep-all",
        help=(
            "How to return CD-HIT results to CSV. "
            "'keep-all' keeps every original row and only adds cluster metadata; "
            "'representative-only' keeps rows whose ori_id is the representative of each cluster."
        ),
    )
    parser.add_argument(
        "--validation-ratio",
        type=float,
        default=0.25,
        help="Validation ratio assigned by cluster, not by row. Default: 0.25 for a 75:25 train/validation split.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed for cluster split.")
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


def parse_cluster_file(clstr_path: Path) -> pd.DataFrame:
    if not clstr_path.exists():
        raise FileNotFoundError(
            f"CD-HIT cluster file not found: {clstr_path}. "
            "Run cd-hit-est first, or pass --clstr with the correct .clstr path."
        )

    records: list[dict[str, object]] = []
    current_cluster: str | None = None
    current_members: list[dict[str, object]] = []

    def flush_cluster() -> None:
        nonlocal current_members, current_cluster
        if current_cluster is None:
            return
        if not current_members:
            return
        cluster_size = len(current_members)
        representative = next(
            (
                str(item["ori_id"])
                for item in current_members
                if bool(item["is_representative"])
            ),
            str(current_members[0]["ori_id"]),
        )
        for item in current_members:
            records.append(
                {
                    "ori_id": item["ori_id"],
                    "cdhit_cluster": current_cluster,
                    "cdhit_cluster_size": cluster_size,
                    "cdhit_cluster_representative_ori_id": representative,
                    "cdhit_is_representative": bool(item["is_representative"]),
                }
            )
        current_members = []

    with open(clstr_path, "r", encoding="utf-8", errors="ignore") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">Cluster "):
                flush_cluster()
                current_cluster = line.split(">Cluster ", 1)[1].strip()
                continue

            if current_cluster is None:
                continue

            if ">" not in line or "..." not in line:
                continue

            header_part = line.split(">", 1)[1].split("...", 1)[0]
            ori_id = header_part.split("|idx=", 1)[0].strip()
            is_representative = line.rstrip().endswith("*")
            current_members.append(
                {
                    "ori_id": ori_id,
                    "is_representative": is_representative,
                }
            )

    flush_cluster()
    return pd.DataFrame(records)


def assign_cluster_split(
    cluster_df: pd.DataFrame,
    validation_ratio: float,
    seed: int,
    weight_col: str = "cdhit_cluster_size",
) -> pd.DataFrame:
    if cluster_df.empty:
        cluster_df = cluster_df.copy()
        cluster_df["split"] = pd.Series(dtype=str)
        return cluster_df

    if weight_col not in cluster_df.columns:
        weight_col = "cdhit_cluster_size"

    cluster_sizes = (
        cluster_df.groupby("cdhit_cluster", as_index=False)[weight_col]
        .sum()
    )
    cluster_sizes = cluster_sizes.copy()
    cluster_sizes["cluster_weight"] = cluster_sizes[weight_col]
    cluster_sizes = cluster_sizes.drop(columns=[weight_col])
    cluster_sizes = sort_cluster_frame(cluster_sizes)

    total_rows = int(cluster_sizes["cluster_weight"].sum())
    if validation_ratio <= 0 or total_rows == 0:
        cluster_df = cluster_df.copy()
        cluster_df["split"] = "train"
        return cluster_df

    shuffled = cluster_sizes.sample(frac=1.0, random_state=seed).reset_index(drop=True)
    target_val_rows = max(1, int(round(total_rows * validation_ratio)))

    val_clusters: set[str] = set()
    collected = 0
    for _, row in shuffled.iterrows():
        cluster_id = str(row["cdhit_cluster"])
        val_clusters.add(cluster_id)
        collected += int(row["cluster_weight"])
        if collected >= target_val_rows:
            break

    cluster_df = cluster_df.copy()
    cluster_df["split"] = cluster_df["cdhit_cluster"].astype(str).map(
        lambda x: "validation" if x in val_clusters else "train"
    )
    return cluster_df


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level, args.log_file)

    if not 0 <= args.validation_ratio < 1:
        raise ValueError("--validation-ratio must be >= 0 and < 1.")

    fasta_path = Path(args.fasta)
    clstr_path = Path(args.clstr) if args.clstr else Path(str(fasta_path) + ".clstr")
    csv_path = Path(args.csv)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Parsing CD-HIT clusters from %s", clstr_path)
    cluster_df = parse_cluster_file(clstr_path)
    if cluster_df.empty:
        raise ValueError(f"No cluster records were parsed from: {clstr_path}")
    logger.info(
        "Parsed %d clustered ori IDs across %d clusters",
        len(cluster_df),
        cluster_df["cdhit_cluster"].nunique(),
    )

    logger.info("Loading annotated CSV from %s", csv_path)
    df = pd.read_csv(csv_path, low_memory=False)
    if args.id_column not in df.columns:
        raise KeyError(f"Column '{args.id_column}' not found in CSV: {csv_path}")
    logger.info("Loaded %d annotated rows", len(df))

    df[args.id_column] = df[args.id_column].fillna("").astype(str).str.strip()
    row_weights = (
        df[df[args.id_column] != ""]
        .groupby(args.id_column, as_index=False)
        .size()
        .rename(columns={"size": "source_row_count"})
        .rename(columns={args.id_column: "_cluster_ori_id"})
    )

    cluster_merge_df = cluster_df.rename(
        columns={"ori_id": "_cluster_ori_id", "split": "_source_cluster_split"}
    )
    cluster_merge_df = cluster_merge_df.merge(
        row_weights,
        how="left",
        left_on="_cluster_ori_id",
        right_on="_cluster_ori_id",
    )
    cluster_merge_df["source_row_count"] = cluster_merge_df["source_row_count"].fillna(1).astype(int)
    logger.info(
        "Assigning cluster-aware split with target train:validation ratio %.0f:%.0f",
        (1 - args.validation_ratio) * 100,
        args.validation_ratio * 100,
    )
    cluster_df = assign_cluster_split(cluster_merge_df, args.validation_ratio, args.seed, weight_col="source_row_count")
    cluster_df["cdhit_split"] = cluster_df["split"]
    cluster_summary = (
        cluster_df.groupby("cdhit_cluster", as_index=False)
        .agg(
            cdhit_cluster_size=("_cluster_ori_id", "size"),
            cdhit_cluster_row_count=("source_row_count", "sum"),
            cdhit_cluster_representative_ori_id=("cdhit_cluster_representative_ori_id", "first"),
            cdhit_is_representative=("cdhit_is_representative", "sum"),
            split=("cdhit_split", "first"),
        )
    )
    cluster_summary = sort_cluster_frame(cluster_summary)

    df = df.merge(cluster_df, how="left", left_on=args.id_column, right_on="_cluster_ori_id")
    df = df.drop(columns=["_cluster_ori_id", "source_row_count"], errors="ignore")

    if "cdhit_split" in df.columns:
        df["split"] = df["cdhit_split"]

    if args.dedup_mode == "representative-only":
        before_filter = len(df)
        df = df[df["cdhit_is_representative"].fillna(False)].copy()
        logger.info(
            "Dedup mode '%s': kept %d rows, removed %d non-representative rows",
            args.dedup_mode,
            len(df),
            before_filter - len(df),
        )
    else:
        logger.info("Dedup mode '%s': keeping all rows", args.dedup_mode)

    missing_cluster = df["cdhit_cluster"].isna().sum()
    if missing_cluster:
        logger.warning("%d rows had no CD-HIT cluster label.", missing_cluster)

    cluster_output_df = cluster_df.rename(columns={"_cluster_ori_id": "ori_id"})

    df.to_csv(output_path, index=False, encoding="utf-8-sig")
    cluster_output_df.to_csv(output_path.with_name("cdhit_cluster_map.csv"), index=False, encoding="utf-8-sig")
    cluster_summary.to_csv(output_path.with_name("cdhit_cluster_summary.csv"), index=False, encoding="utf-8-sig")

    train_rows = int((df["cdhit_split"] == "train").sum()) if "cdhit_split" in df.columns else 0
    val_rows = int((df["cdhit_split"] == "validation").sum()) if "cdhit_split" in df.columns else 0
    total_split_rows = train_rows + val_rows
    train_pct = (train_rows / total_split_rows * 100) if total_split_rows else 0
    val_pct = (val_rows / total_split_rows * 100) if total_split_rows else 0

    logger.info("Done.")
    logger.info("Annotated rows: %d", len(df))
    logger.info("Unique ori IDs with cluster labels: %d", cluster_output_df["ori_id"].nunique())
    logger.info("Clusters: %d", cluster_output_df["cdhit_cluster"].nunique())
    logger.info("Train rows: %d (%.2f%%)", train_rows, train_pct)
    logger.info("Validation rows: %d (%.2f%%)", val_rows, val_pct)
    logger.info("Output file: %s", output_path)
    logger.info("Cluster map: %s", output_path.with_name("cdhit_cluster_map.csv"))
    logger.info("Cluster summary: %s", output_path.with_name("cdhit_cluster_summary.csv"))


if __name__ == "__main__":
    main()
