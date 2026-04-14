from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge CD-HIT cluster labels back to the annotated CSV without dropping members."
    )
    parser.add_argument("--fasta", default="cd_hit_process.fasta", help="CD-HIT output FASTA file.")
    parser.add_argument(
        "--clstr",
        default=None,
        help="CD-HIT cluster file. Default: same path as FASTA with .clstr appended.",
    )
    parser.add_argument(
        "--csv",
        default="merged_final_optimized_replication_annotated.csv",
        help="Source annotated CSV file.",
    )
    parser.add_argument("--output", default="final_data.csv", help="Output CSV file.")
    parser.add_argument("--id-column", default="ori_id", help="ID column in source CSV.")
    parser.add_argument(
        "--validation-ratio",
        type=float,
        default=0.05,
        help="Validation ratio assigned by cluster, not by row.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed for cluster split.")
    return parser.parse_args()


def parse_cluster_file(clstr_path: Path) -> pd.DataFrame:
    if not clstr_path.exists():
        raise FileNotFoundError(f"CD-HIT cluster file not found: {clstr_path}")

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


def assign_cluster_split(cluster_df: pd.DataFrame, validation_ratio: float, seed: int) -> pd.DataFrame:
    if cluster_df.empty:
        cluster_df = cluster_df.copy()
        cluster_df["split"] = pd.Series(dtype=str)
        return cluster_df

    cluster_sizes = (
        cluster_df[["cdhit_cluster", "cdhit_cluster_size"]]
        .drop_duplicates(subset=["cdhit_cluster"])
        .copy()
    )
    cluster_sizes = cluster_sizes.sort_values("cdhit_cluster").reset_index(drop=True)

    total_rows = int(cluster_sizes["cdhit_cluster_size"].sum())
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
        collected += int(row["cdhit_cluster_size"])
        if collected >= target_val_rows:
            break

    cluster_df = cluster_df.copy()
    cluster_df["split"] = cluster_df["cdhit_cluster"].astype(str).map(
        lambda x: "validation" if x in val_clusters else "train"
    )
    return cluster_df


def main() -> None:
    args = parse_args()
    fasta_path = Path(args.fasta)
    clstr_path = Path(args.clstr) if args.clstr else Path(str(fasta_path) + ".clstr")
    csv_path = Path(args.csv)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print("Parsing CD-HIT clusters...")
    cluster_df = parse_cluster_file(clstr_path)
    if cluster_df.empty:
        raise ValueError(f"No cluster records were parsed from: {clstr_path}")

    cluster_df = assign_cluster_split(cluster_df, args.validation_ratio, args.seed)
    cluster_summary = (
        cluster_df.groupby("cdhit_cluster", as_index=False)
        .agg(
            cdhit_cluster_size=("ori_id", "size"),
            cdhit_cluster_representative_ori_id=("cdhit_cluster_representative_ori_id", "first"),
            cdhit_is_representative=("cdhit_is_representative", "sum"),
            split=("split", "first"),
        )
        .sort_values("cdhit_cluster", key=lambda s: s.astype(int))
        .reset_index(drop=True)
    )

    print("Loading annotated CSV...")
    df = pd.read_csv(csv_path, low_memory=False)
    if args.id_column not in df.columns:
        raise KeyError(f"Column '{args.id_column}' not found in CSV: {csv_path}")

    df[args.id_column] = df[args.id_column].fillna("").astype(str).str.strip()

    cluster_merge_df = cluster_df.rename(columns={"ori_id": "_cluster_ori_id"})
    df = df.merge(cluster_merge_df, how="left", left_on=args.id_column, right_on="_cluster_ori_id")
    df = df.drop(columns=["_cluster_ori_id"])

    missing_cluster = df["cdhit_cluster"].isna().sum()
    if missing_cluster:
        print(f"Warning: {missing_cluster} rows had no CD-HIT cluster label.")

    df.to_csv(output_path, index=False, encoding="utf-8-sig")
    cluster_df.to_csv(output_path.with_name("cdhit_cluster_map.csv"), index=False, encoding="utf-8-sig")
    cluster_summary.to_csv(output_path.with_name("cdhit_cluster_summary.csv"), index=False, encoding="utf-8-sig")

    print("Done.")
    print(f"Annotated rows: {len(df)}")
    print(f"Unique ori IDs with cluster labels: {cluster_df['ori_id'].nunique()}")
    print(f"Clusters: {cluster_df['cdhit_cluster'].nunique()}")
    print(f"Validation rows: {(df['split'] == 'validation').sum() if 'split' in df.columns else 0}")
    print(f"Output file: {output_path}")


if __name__ == "__main__":
    main()
