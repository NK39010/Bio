from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


SPECIAL_TOKENS = ["[START]", "[END]", "[PAD]", "[UNK]", "[SEP]"]
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
NUCLEOTIDES = list("atcg")
MECHANISM_TOKENS = [
    "RCR",
    "theta",
    "unresolved",
]
DEFAULT_MODEL_MAX_LENGTH = 10**30


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build tokenizer files from the final dataset.")
    parser.add_argument(
        "--input",
        default="final_data.csv",
        help="Input CSV or Parquet file.",
    )
    parser.add_argument(
        "--tokenizer-json",
        default="tokenizer/tokenizer.json",
        help="Output tokenizer JSON path.",
    )
    parser.add_argument(
        "--species-min-count",
        type=int,
        default=6,
        help="Keep only species names appearing more than 5 times by default.",
    )
    parser.add_argument("--tokenizer-config-json", default=None, help="Optional tokenizer_config.json path.")
    parser.add_argument("--special-tokens-map-json", default=None, help="Optional special_tokens_map.json path.")
    return parser.parse_args()


def load_table(path: str) -> pd.DataFrame:
    file_path = Path(path)
    if file_path.suffix.lower() == ".parquet":
        return pd.read_parquet(file_path)
    return pd.read_csv(file_path, low_memory=False)


def normalize_species_name(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    if not text:
        return ""
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^0-9A-Za-z_.:-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text


def make_token_entry(token: str, token_id: int, special: bool) -> dict[str, object]:
    return {
        "id": token_id,
        "content": token,
        "single_word": False,
        "lstrip": False,
        "rstrip": False,
        "normalized": False if special else True,
        "special": special,
    }


def make_token_meta(token: str, special: bool) -> dict[str, object]:
    return {
        "content": token,
        "single_word": False,
        "lstrip": False,
        "rstrip": False,
        "normalized": False if special else True,
        "special": special,
    }


def make_special_map_entry(token: str) -> dict[str, object]:
    return {
        "content": token,
        "single_word": False,
        "lstrip": False,
        "rstrip": False,
        "normalized": False,
    }


def build_tokenizer_json(base_vocab: list[str], added_tokens: list[str]) -> dict[str, object]:
    vocab = {token: idx for idx, token in enumerate(base_vocab)}
    next_id = len(vocab)
    tokenizer_added_tokens = [make_token_entry(token, idx, True) for idx, token in enumerate(SPECIAL_TOKENS)]
    tokenizer_added_tokens.extend(
        make_token_entry(token, next_id + offset, False) for offset, token in enumerate(added_tokens)
    )

    return {
        "version": "1.0",
        "truncation": None,
        "padding": None,
        "added_tokens": tokenizer_added_tokens,
        "normalizer": None,
        "pre_tokenizer": {"type": "Whitespace"},
        "post_processor": None,
        "decoder": None,
        "model": {
            "type": "BPE",
            "dropout": None,
            "unk_token": "[UNK]",
            "continuing_subword_prefix": None,
            "end_of_word_suffix": None,
            "fuse_unk": False,
            "byte_fallback": False,
            "ignore_merges": False,
            "vocab": vocab,
            "merges": [],
        },
    }


def build_tokenizer_config_json(added_tokens: list[str]) -> dict[str, object]:
    added_tokens_decoder = {}
    for idx, token in enumerate(SPECIAL_TOKENS):
        added_tokens_decoder[str(idx)] = make_token_meta(token, True)
    base_offset = len(SPECIAL_TOKENS)
    for offset, token in enumerate(added_tokens):
        token_id = base_offset + offset
        added_tokens_decoder[str(token_id)] = make_token_meta(token, False)

    return {
        "added_tokens_decoder": added_tokens_decoder,
        "bos_token": "[START]",
        "clean_up_tokenization_spaces": True,
        "eos_token": "[END]",
        "extra_special_tokens": {},
        "model_max_length": DEFAULT_MODEL_MAX_LENGTH,
        "pad_token": "[PAD]",
        "sep_token": "[SEP]",
        "tokenizer_class": "PreTrainedTokenizerFast",
        "unk_token": "[UNK]",
    }


def build_special_tokens_map_json() -> dict[str, object]:
    return {
        "bos_token": make_special_map_entry("[START]"),
        "eos_token": make_special_map_entry("[END]"),
        "pad_token": make_special_map_entry("[PAD]"),
        "sep_token": make_special_map_entry("[SEP]"),
        "unk_token": make_special_map_entry("[UNK]"),
    }


def main() -> None:
    args = parse_args()
    df = load_table(args.input)

    if "species" not in df.columns:
        raise KeyError("Input file must contain a 'species' column.")

    species_series = df["species"].fillna("").astype(str).str.strip()
    species_series = species_series[species_series != ""]
    species_series = species_series[species_series.str.lower() != "unknown"]
    species_counts = species_series.value_counts()
    species_vocab = [
        normalize_species_name(species)
        for species in species_counts[species_counts >= args.species_min_count].index
        if normalize_species_name(species)
    ]
    species_vocab = sorted(dict.fromkeys(species_vocab))

    base_vocab = SPECIAL_TOKENS + AMINO_ACIDS + NUCLEOTIDES
    added_tokens = MECHANISM_TOKENS + species_vocab

    tokenizer_json = build_tokenizer_json(base_vocab, added_tokens)
    tokenizer_config_json = build_tokenizer_config_json(added_tokens)
    special_tokens_map_json = build_special_tokens_map_json()

    tokenizer_json_path = Path(args.tokenizer_json)
    tokenizer_json_path.parent.mkdir(parents=True, exist_ok=True)
    tokenizer_json_path.write_text(
        json.dumps(tokenizer_json, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    tokenizer_config_path = Path(args.tokenizer_config_json) if args.tokenizer_config_json else tokenizer_json_path.parent / "tokenizer_config.json"
    tokenizer_config_path.parent.mkdir(parents=True, exist_ok=True)
    tokenizer_config_path.write_text(
        json.dumps(tokenizer_config_json, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    special_tokens_map_path = (
        Path(args.special_tokens_map_json)
        if args.special_tokens_map_json
        else tokenizer_json_path.parent / "special_tokens_map.json"
    )
    special_tokens_map_path.parent.mkdir(parents=True, exist_ok=True)
    special_tokens_map_path.write_text(
        json.dumps(special_tokens_map_json, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print(f"tokenizer JSON: {tokenizer_json_path}")
    print(f"tokenizer config: {tokenizer_config_path}")
    print(f"special tokens map: {special_tokens_map_path}")
    print(f"special tokens: {len(SPECIAL_TOKENS)}")
    print(f"amino acids: {len(AMINO_ACIDS)}")
    print(f"nucleotides: {len(NUCLEOTIDES)}")
    print(f"mechanism tokens: {len(MECHANISM_TOKENS)}")
    print(f"species tokens (>={args.species_min_count}): {len(species_vocab)}")
    print(f"base vocab: {len(base_vocab)}")
    print(f"added tokens: {len(added_tokens)}")


if __name__ == "__main__":
    main()
