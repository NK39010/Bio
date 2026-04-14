from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


EXPECTED_CDHIT_COLS = [
    "cdhit_cluster",
    "cdhit_cluster_size",
    "cdhit_cluster_representative_ori_id",
    "cdhit_is_representative",
    "cdhit_split",
]


def clean_seq(seq, is_aa=False):
    if pd.isna(seq):
        return ""
    seq = str(seq).strip()
    return seq.upper() if is_aa else seq.lower()


def normalize_mechanism_call(value):
    if pd.isna(value):
        text = ""
    else:
        text = str(value).strip()
    if text == "" or text.lower() == "unknown":
        return "unresolved"
    if text == "RCR_like":
        return "RCR"
    if text == "theta_like":
        return "theta"
    return text


def parse_args():
    parser = argparse.ArgumentParser(description="Build split training CSV and parquet files.")
    parser.add_argument("--input", default="final_data.csv", help="Input CSV file.")
    parser.add_argument(
        "--output-dir",
        default="training_data",
        help="Directory where train/validation CSV and parquet files will be written.",
    )
    return parser.parse_args()


def normalize_split_columns(df):
    df = df.copy()

    if "split" not in df.columns:
        if "cdhit_split" in df.columns:
            df["split"] = df["cdhit_split"]
        elif "split_y" in df.columns:
            df["split"] = df["split_y"]
        elif "split_x" in df.columns:
            df["split"] = df["split_x"]

    if "cdhit_split" not in df.columns and "split" in df.columns:
        df["cdhit_split"] = df["split"]

    return df


def split_dataframe(df):
    df = df.copy()

    split_source = None
    for col in ("cdhit_split", "split"):
        if col in df.columns:
            split_source = col
            break

    if split_source is None:
        df["split"] = "train"
        df["cdhit_split"] = "train"
        return df[df["split"] == "train"], df.iloc[0:0].copy()

    split_values = df[split_source].fillna("").astype(str).str.strip().str.lower()
    val_mask = split_values.isin({"validation", "val", "dev", "eval", "test"})
    train_mask = ~val_mask

    df.loc[train_mask, "split"] = "train"
    df.loc[val_mask, "split"] = "validation"
    df["cdhit_split"] = df["split"]

    train_df = df.loc[train_mask].copy()
    val_df = df.loc[val_mask].copy()
    return train_df, val_df


def finalize_table(df):
    df = df.copy()

    keep_cols = [
        "OriC sequence",
        "ori_id",
        "plasmid_id",
        "species",
        "source",
        "pfamid_fast",
        "rep_id",
        "Rep_type_fast",
        "rep_seq",
        "rep_dna_seq",
        "full_replicon_seq",
        "split",
        "cdhit_split",
        "cdhit_cluster",
        "cdhit_cluster_size",
        "cdhit_cluster_representative_ori_id",
        "cdhit_is_representative",
        "pairing_method",
        "ori_midpoint",
        "rep_midpoint",
        "replicon_length",
        "__index_level_0__",
        "replication_mechanism_call",
    ]
    existing_cols = [col for col in keep_cols if col in df.columns]
    df_out = df[existing_cols].copy()

    if "OriC sequence" in df_out.columns:
        df_out["OriC sequence"] = df_out["OriC sequence"].apply(lambda x: clean_seq(x, is_aa=False))
    if "rep_seq" in df_out.columns:
        df_out["rep_seq"] = df_out["rep_seq"].apply(lambda x: clean_seq(x, is_aa=True))
    if "rep_dna_seq" in df_out.columns:
        df_out["rep_dna_seq"] = df_out["rep_dna_seq"].apply(lambda x: clean_seq(x, is_aa=False))
    if "full_replicon_seq" in df_out.columns:
        df_out["full_replicon_seq"] = df_out["full_replicon_seq"].apply(lambda x: clean_seq(x, is_aa=False))

    if "replication_mechanism_call" in df_out.columns:
        df_out["replication_mechanism_call"] = df_out["replication_mechanism_call"].apply(
            normalize_mechanism_call
        )
        df_out["replication_mechanism_term"] = (
            df_out["replication_mechanism_call"]
            .fillna("unresolved")
            .astype(str)
            .str.strip()
            .replace("", "unresolved")
        )
    else:
        df_out["replication_mechanism_term"] = "unresolved"

    return df_out


def write_outputs(df, output_dir: Path, split_name: str):
    csv_path = output_dir / f"{split_name}.csv"
    parquet_path = output_dir / f"{split_name}.parquet"

    df.to_csv(csv_path, index=False, encoding="utf-8-sig")
    table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(table, parquet_path, compression="snappy")
    return csv_path, parquet_path


def main():
    args = parse_args()
    df = pd.read_csv(args.input, low_memory=False)
    df = normalize_split_columns(df)

    missing_cdhit_cols = [col for col in EXPECTED_CDHIT_COLS if col not in df.columns]
    if missing_cdhit_cols:
        print(
            "⚠️ CD-HIT columns missing from input; split files will not contain complete "
            f"cluster metadata: {', '.join(missing_cdhit_cols)}"
        )

    train_df, val_df = split_dataframe(df)
    train_out = finalize_table(train_df)
    val_out = finalize_table(val_df)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    train_csv, train_parquet = write_outputs(train_out, output_dir, "train")
    val_csv, val_parquet = write_outputs(val_out, output_dir, "validation")

    print(f"train csv: {train_csv}")
    print(f"validation csv: {val_csv}")
    print(f"train parquet: {train_parquet}")
    print(f"validation parquet: {val_parquet}")
    print(f"train rows: {len(train_out)}")
    print(f"validation rows: {len(val_out)}")


if __name__ == "__main__":
    main()
