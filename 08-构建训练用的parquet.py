from __future__ import annotations

import argparse

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


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
    parser = argparse.ArgumentParser(description="Build training parquet file.")
    parser.add_argument("--input", default="final_data.csv", help="Input CSV file.")
    parser.add_argument("--parquet", default="final_data.parquet", help="Output parquet file.")
    return parser.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.input, low_memory=False)

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
        "cdhit_cluster",
        "cdhit_cluster_size",
        "cdhit_cluster_representative_ori_id",
        "cdhit_is_representative",
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
        df_out["replication_mechanism_call"] = (
            df_out["replication_mechanism_call"].apply(normalize_mechanism_call)
        )
        df_out["replication_mechanism_term"] = df_out["replication_mechanism_call"].fillna("unresolved").astype(str).str.strip().replace("", "unresolved")
    else:
        df_out["replication_mechanism_term"] = "unresolved"

    table = pa.Table.from_pandas(df_out, preserve_index=False)
    pq.write_table(table, args.parquet, compression="snappy")

    print(f"✅ parquet: {args.parquet}")
    print(f"📊 共 {len(df_out)} 条数据")


if __name__ == "__main__":
    main()
