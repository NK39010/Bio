from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


AA_SET = set("ACDEFGHIKLMNPQRSTVWY")
NT_SET = set("ATCGN")


def parse_args():
    parser = argparse.ArgumentParser(description="Build semantic token JSON files from final CSV.")
    parser.add_argument("--input", default="final_data.csv", help="Input CSV file.")
    parser.add_argument(
        "--token-base-json",
        default="final_data_tokens_without_mechanism.json",
        help="Output token json without mechanism term.",
    )
    parser.add_argument(
        "--token-with-mechanism-json",
        default="final_data_tokens_with_mechanism.json",
        help="Output token json with mechanism term.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.input, low_memory=False)

    if "ori_id" not in df.columns or "OriC sequence" not in df.columns:
        raise KeyError("CSV must contain 'ori_id' and 'OriC sequence' to generate token files.")

    if "replication_mechanism_term" not in df.columns:
        mechanism_col = "replication_mechanism_call"
        if mechanism_col in df.columns:
            df["replication_mechanism_term"] = (
                "replication_mechanism:"
                + df[mechanism_col].fillna("unknown").astype(str).str.strip().replace("", "unknown")
            )
        else:
            df["replication_mechanism_term"] = "replication_mechanism:unknown"

    species_col = "species"
    if species_col in df.columns:
        df["species_term"] = "species:" + df[species_col].fillna("unknown").astype(str).str.strip().replace("", "unknown")
    else:
        df["species_term"] = "species:unknown"

    aa_col = "rep_seq"
    if aa_col not in df.columns:
        df[aa_col] = ""

    token_df = df[
        ["ori_id", "OriC sequence", "replication_mechanism_term", "species_term", "rep_seq"]
    ].dropna(subset=["ori_id", "OriC sequence"]).copy()

    def nt_tokens(seq: str) -> list[str]:
        seq = str(seq).upper()
        return [f"nt:{ch}" for ch in seq if ch in NT_SET]

    def aa_tokens(seq: str) -> list[str]:
        seq = str(seq).upper()
        return [f"aa:{ch}" for ch in seq if ch in AA_SET]

    token_df["nt_tokens"] = token_df["OriC sequence"].apply(nt_tokens)
    token_df["aa_tokens"] = token_df["rep_seq"].apply(aa_tokens)

    token_df["tokens_base"] = token_df.apply(
        lambda row: [row["species_term"]] + row["aa_tokens"] + row["nt_tokens"],
        axis=1,
    )
    token_df["tokens_with_mechanism"] = token_df.apply(
        lambda row: [row["replication_mechanism_term"]] + row["tokens_base"],
        axis=1,
    )

    token_base_records = token_df[["ori_id", "tokens_base"]].rename(
        columns={"tokens_base": "tokens"}
    ).to_dict(orient="records")
    token_with_mech_records = token_df[
        ["ori_id", "species_term", "replication_mechanism_term", "tokens_with_mechanism"]
    ].rename(columns={"tokens_with_mechanism": "tokens"}).to_dict(orient="records")

    Path(args.token_base_json).write_text(
        json.dumps(token_base_records, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    Path(args.token_with_mechanism_json).write_text(
        json.dumps(token_with_mech_records, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(f"✅ tokens(no mechanism): {args.token_base_json}")
    print(f"✅ tokens(with mechanism): {args.token_with_mechanism_json}")
    print(f"📊 共 {len(token_df)} 条词元数据")


if __name__ == "__main__":
    main()