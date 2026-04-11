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
    parser.add_argument(
        "--tokenizer-json",
        default=None,
        help="Optional tokenizer JSON output path.",
    )
    parser.add_argument(
        "--tokenizer-vocab-txt",
        default=None,
        help="Optional tokenizer vocab text output path.",
    )
    return parser.parse_args()


def normalize_mechanism(value: object) -> str:
    if pd.isna(value):
        text = ""
    else:
        text = str(value).strip()
    if text == "" or text.lower() == "unknown":
        return "unresolved"
    if text.startswith("replication_mechanism:"):
        value = text.split(":", 1)[1]
    else:
        value = text
    if value == "RCR_like":
        value = "RCR"
    elif value == "theta_like":
        value = "theta"
    if value == "" or value.lower() == "unknown":
        value = "unresolved"
    return value


def main():
    args = parse_args()
    df = pd.read_csv(args.input, low_memory=False)

    if "ori_id" not in df.columns or "OriC sequence" not in df.columns:
        raise KeyError("CSV must contain 'ori_id' and 'OriC sequence' to generate token files.")

    if "replication_mechanism_term" in df.columns:
        df["replication_mechanism_term"] = df["replication_mechanism_term"].apply(normalize_mechanism)
    else:
        mechanism_col = "replication_mechanism_call"
        if mechanism_col in df.columns:
            df["replication_mechanism_term"] = df[mechanism_col].apply(normalize_mechanism)
        else:
            df["replication_mechanism_term"] = "unresolved"

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
        return [ch.lower() for ch in seq if ch in NT_SET]

    def aa_tokens(seq: str) -> list[str]:
        seq = str(seq).upper()
        return [ch for ch in seq if ch in AA_SET]

    token_df["nt_tokens"] = token_df["OriC sequence"].apply(nt_tokens)
    token_df["aa_tokens"] = token_df["rep_seq"].apply(aa_tokens)

    token_df["tokens_base"] = token_df.apply(
        lambda row: [row["species_term"]] + row["aa_tokens"] + row["nt_tokens"],
        axis=1,
    )
    token_df["tokens_with_mechanism"] = token_df.apply(
        lambda row: [f"replication_mechanism:{row['replication_mechanism_term']}"] + row["tokens_base"],
        axis=1,
    )

    token_base_records = token_df[["ori_id", "tokens_base"]].rename(
        columns={"tokens_base": "tokens"}
    ).to_dict(orient="records")
    token_with_mech_records = token_df[
        ["ori_id", "species_term", "replication_mechanism_term", "tokens_with_mechanism"]
    ].rename(columns={"tokens_with_mechanism": "tokens"}).to_dict(orient="records")

    def build_tokenizer_json(tokens: list[str]) -> dict[str, object]:
        special_tokens = ["[START]", "[END]", "[PAD]", "[UNK]", "[SEP]"]
        vocab = {token: idx for idx, token in enumerate(special_tokens)}
        for token in tokens:
            if token not in vocab:
                vocab[token] = len(vocab)

        return {
            "version": "1.0",
            "truncation": None,
            "padding": None,
            "added_tokens": [
                {
                    "id": idx,
                    "content": token,
                    "single_word": False,
                    "lstrip": False,
                    "rstrip": False,
                    "normalized": False,
                    "special": True,
                }
                for idx, token in enumerate(special_tokens)
            ],
            "normalizer": None,
            "pre_tokenizer": {"type": "Whitespace"},
            "post_processor": None,
            "decoder": None,
            "model": {
                "type": "WordLevel",
                "unk_token": "[UNK]",
                "vocab": vocab,
            },
        }
    
    Path(args.token_base_json).write_text(
        json.dumps(token_base_records, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    Path(args.token_with_mechanism_json).write_text(
        json.dumps(token_with_mech_records, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    if args.tokenizer_json or args.tokenizer_vocab_txt:
        tokenizer_tokens = set()
        for records in (token_base_records, token_with_mech_records):
            for record in records:
                tokenizer_tokens.update(record["tokens"])
        tokenizer_tokens = sorted(tokenizer_tokens)

        if args.tokenizer_json:
            tokenizer_data = build_tokenizer_json(tokenizer_tokens)
            Path(args.tokenizer_json).write_text(
                json.dumps(tokenizer_data, ensure_ascii=False, indent=2), encoding="utf-8"
            )
            print(f"✅ tokenizer JSON: {args.tokenizer_json}")

        if args.tokenizer_vocab_txt:
            Path(args.tokenizer_vocab_txt).write_text(
                "\n".join(tokenizer_tokens), encoding="utf-8"
            )
            print(f"✅ tokenizer vocab: {args.tokenizer_vocab_txt}")

    print(f"✅ tokens(no mechanism): {args.token_base_json}")
    print(f"✅ tokens(with mechanism): {args.token_with_mechanism_json}")
    print(f"📊 共 {len(token_df)} 条词元数据")


if __name__ == "__main__":
    main()