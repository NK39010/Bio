#!/usr/bin/env python
# coding=utf-8
"""
OriGen training script.

This version is refactored for clarity and adds an explicit switch for whether
to include the replication-mechanism token in each training example.

Input format expected by the training parquet:
    - host species column, e.g. `host_species` or `species`
    - Rep protein column, e.g. `rep_protein` or `rep_seq`
    - oriV DNA sequence column, e.g. `oriv_sequence` or `rep_dna_seq`
    - optional replication mechanism column, e.g. `replication_mechanism_term`

Examples:
    python 10-train_origen.py \
      --tokenizer-path tokenizer \
      --dataset-path final_data.parquet \
      --species-col host_species \
      --rep-col rep_protein \
      --seq-col oriv_sequence \
      --include-mechanism-token \
      --mechanism-col replication_mechanism_term

    python 10-train_origen.py \
      --tokenizer-path tokenizer \
      --dataset-path final_data.parquet \
      --species-col host_species \
      --rep-col rep_protein \
      --seq-col oriv_sequence
"""

from __future__ import annotations

import logging
import math
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List

from datasets import DatasetDict, load_dataset
from transformers import (
    AutoModelForCausalLM,
    AutoTokenizer,
    DataCollatorForLanguageModeling,
    HfArgumentParser,
    RoFormerConfig,
    Trainer,
    TrainingArguments,
    set_seed,
)
from transformers.trainer_utils import get_last_checkpoint


logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)


@dataclass
class ModelDataArguments:
    """User-provided data paths and model architecture parameters."""

    tokenizer_path: str = field(
        metadata={"help": "Path to the local tokenizer folder containing tokenizer.json"}
    )
    dataset_path: str = field(
        metadata={"help": "Path to a parquet file, a directory of parquet files, or a glob pattern"}
    )

    species_col: str = field(default="host_species", metadata={"help": "Species column name"})
    rep_col: str = field(default="rep_protein", metadata={"help": "Rep protein column name"})
    seq_col: str = field(default="oriv_sequence", metadata={"help": "oriV DNA sequence column name"})

    include_mechanism_token: bool = field(
        default=False,
        metadata={"help": "Append the replication mechanism token to each example"},
    )
    mechanism_col: str = field(
        default="replication_mechanism_term",
        metadata={"help": "Replication-mechanism column name"},
    )

    model_size: int = field(default=768, metadata={"help": "Hidden size"})
    num_hidden_layers: int = field(default=12, metadata={"help": "Number of transformer layers"})
    num_attention_heads: int = field(default=12, metadata={"help": "Number of attention heads"})
    max_seq_len: int = field(default=1500, metadata={"help": "Maximum sequence length"})

    do_lower_dna: bool = field(default=True, metadata={"help": "Lower-case DNA sequence before tokenization"})
    validation_split_ratio: float = field(
        default=0.05,
        metadata={"help": "Validation split ratio used when no validation set is provided"},
    )


class CustomTrainingArguments(TrainingArguments):
    """Standard Hugging Face training arguments."""


def normalize_text(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() == "nan":
        return ""
    return text


def format_example(example: Dict[str, Any], args: ModelDataArguments) -> str:
    species = normalize_text(example.get(args.species_col, ""))
    rep = normalize_text(example.get(args.rep_col, ""))
    seq = normalize_text(example.get(args.seq_col, ""))

    if not seq:
        return ""

    if args.do_lower_dna:
        seq = seq.lower()

    parts: List[str] = []
    if species:
        parts.append(species)
    if rep:
        parts.append(rep)

    if args.include_mechanism_token:
        mech = normalize_text(example.get(args.mechanism_col, ""))
        if not mech:
            mech = "unresolved"
        parts.append(mech)

    parts.append(seq)
    return " ".join(parts)


def ensure_tokenizer(tokenizer_path: str):
    tokenizer = AutoTokenizer.from_pretrained(
        tokenizer_path,
        trust_remote_code=True,
        padding_side="right",
    )

    if tokenizer.pad_token is None:
        if tokenizer.eos_token is not None:
            tokenizer.pad_token = tokenizer.eos_token
        else:
            tokenizer.add_special_tokens({"pad_token": "[PAD]"})

    if tokenizer.bos_token is None:
        if tokenizer.eos_token is not None:
            tokenizer.bos_token = tokenizer.eos_token
        else:
            tokenizer.add_special_tokens({"bos_token": "[START]"})

    return tokenizer


def resolve_parquet_files(dataset_path: str) -> Dict[str, List[str]]:
    path = Path(dataset_path)

    if path.is_file():
        return {"train": [str(path)]}

    if path.is_dir():
        parquet_files = sorted(str(p) for p in path.glob("*.parquet"))
        if not parquet_files:
            raise FileNotFoundError(f"No .parquet files found in directory: {dataset_path}")

        train_files = [p for p in parquet_files if Path(p).name.startswith(("train", "training"))]
        val_files = [p for p in parquet_files if Path(p).name.startswith(("validation", "val", "test"))]

        if not train_files and not val_files:
            train_files = parquet_files
            val_files = []
        elif not train_files:
            train_files = [p for p in parquet_files if p not in val_files] or parquet_files

        data_files: Dict[str, List[str]] = {"train": train_files}
        if val_files:
            data_files["validation"] = val_files
        return data_files

    if any(ch in dataset_path for ch in ["*", "?", "["]):
        return {"train": [dataset_path]}

    raise FileNotFoundError(f"Dataset path does not exist: {dataset_path}")


def load_raw_datasets(dataset_path: str, seed: int, validation_split_ratio: float) -> DatasetDict:
    data_files = resolve_parquet_files(dataset_path)
    raw = load_dataset("parquet", data_files=data_files)

    for split_column in ("split", "cdhit_split"):
        if "train" in raw and split_column in raw["train"].column_names:
            split_ds = raw["train"]
            train_ds = split_ds.filter(
                lambda example: normalize_text(example.get(split_column, "")).lower() in {"train", "training"}
            )
            val_ds = split_ds.filter(
                lambda example: normalize_text(example.get(split_column, "")).lower()
                in {"validation", "val", "dev", "eval", "test"}
            )
            if len(train_ds) > 0 and len(val_ds) > 0:
                logger.info(
                    "Using existing split column '%s' from dataset: train=%d, validation=%d",
                    split_column,
                    len(train_ds),
                    len(val_ds),
                )
                return DatasetDict({"train": train_ds, "validation": val_ds})

    if "validation" not in raw:
        logger.info(
            "No validation split found, splitting %.2f%% from train.",
            validation_split_ratio * 100,
        )
        split = raw["train"].train_test_split(test_size=validation_split_ratio, seed=seed)
        raw = DatasetDict({"train": split["train"], "validation": split["test"]})

    return raw


def add_text_column(raw_datasets: DatasetDict, args: ModelDataArguments) -> DatasetDict:
    def build_text_batch(examples: Dict[str, List[Any]]) -> Dict[str, List[str]]:
        texts: List[str] = []
        for idx in range(len(next(iter(examples.values())))):
            example = {key: values[idx] for key, values in examples.items()}
            texts.append(format_example(example, args))
        return {"text": texts}

    return raw_datasets.map(build_text_batch, batched=True, desc="Building training text")


def tokenize_datasets(
    raw_datasets: DatasetDict,
    tokenizer,
    model_args: ModelDataArguments,
) -> DatasetDict:
    if "text" not in raw_datasets["train"].column_names:
        raise ValueError("Text column missing before tokenization.")

    non_empty = raw_datasets.filter(lambda example: bool(normalize_text(example["text"])), desc="Filtering empty text")

    def tokenize_batch(examples: Dict[str, List[str]]) -> Dict[str, Any]:
        return tokenizer(
            examples["text"],
            truncation=True,
            max_length=model_args.max_seq_len,
            padding=False,
            add_special_tokens=True,
        )

    tokenized = non_empty.map(
        tokenize_batch,
        batched=True,
        remove_columns=non_empty["train"].column_names,
        desc="Tokenizing",
    )
    return tokenized


def build_model(tokenizer, model_args: ModelDataArguments):
    config = RoFormerConfig(
        vocab_size=len(tokenizer),
        hidden_size=model_args.model_size,
        num_hidden_layers=model_args.num_hidden_layers,
        num_attention_heads=model_args.num_attention_heads,
        intermediate_size=model_args.model_size * 4,
        hidden_act="gelu",
        attention_probs_dropout_prob=0.1,
        hidden_dropout_prob=0.1,
        max_position_embeddings=model_args.max_seq_len + 512,
        type_vocab_size=0,
        is_decoder=True,
        pad_token_id=tokenizer.pad_token_id,
        eos_token_id=tokenizer.eos_token_id,
        bos_token_id=tokenizer.bos_token_id,
        rotary_value=True,
    )

    model = AutoModelForCausalLM.from_config(config)
    if hasattr(model, "resize_token_embeddings"):
        model.resize_token_embeddings(len(tokenizer))
    return model


def main() -> None:
    parser = HfArgumentParser((ModelDataArguments, CustomTrainingArguments))
    if len(sys.argv) == 2 and sys.argv[1].endswith(".json"):
        model_args, training_args = parser.parse_json_file(json_file=os.path.abspath(sys.argv[1]))
    else:
        model_args, training_args = parser.parse_args_into_dataclasses()

    if not os.path.exists(model_args.tokenizer_path):
        raise FileNotFoundError(f"Tokenizer path does not exist: {model_args.tokenizer_path}")

    set_seed(training_args.seed)

    logger.info("Loading tokenizer from %s", model_args.tokenizer_path)
    tokenizer = ensure_tokenizer(model_args.tokenizer_path)
    logger.info("Tokenizer vocab size: %d", len(tokenizer))

    logger.info("Loading dataset from %s", model_args.dataset_path)
    raw_datasets = load_raw_datasets(
        dataset_path=model_args.dataset_path,
        seed=training_args.seed,
        validation_split_ratio=model_args.validation_split_ratio,
    )

    needed_columns = [model_args.species_col, model_args.rep_col, model_args.seq_col]
    if model_args.include_mechanism_token:
        needed_columns.append(model_args.mechanism_col)

    missing = [col for col in needed_columns if col not in raw_datasets["train"].column_names]
    if missing:
        raise KeyError(
            "Missing required columns in dataset: "
            + ", ".join(missing)
            + f". Available columns: {raw_datasets['train'].column_names}"
        )

    logger.info(
        "Dataset sizes - train: %d, validation: %d",
        len(raw_datasets["train"]),
        len(raw_datasets["validation"]),
    )

    with_text = add_text_column(raw_datasets, model_args)
    tokenized_datasets = tokenize_datasets(with_text, tokenizer, model_args)

    logger.info("Initializing RoFormer causal LM")
    model = build_model(tokenizer, model_args)
    param_count = sum(p.numel() for p in model.parameters()) / 1_000_000
    logger.info("Model parameters: %.2fM", param_count)

    data_collator = DataCollatorForLanguageModeling(
        tokenizer=tokenizer,
        mlm=False,
        pad_to_multiple_of=8 if training_args.fp16 or training_args.bf16 else None,
    )

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=tokenized_datasets["train"],
        eval_dataset=tokenized_datasets["validation"],
        tokenizer=tokenizer,
        data_collator=data_collator,
    )

    if training_args.do_train:
        checkpoint = None
        if os.path.isdir(training_args.output_dir):
            checkpoint = get_last_checkpoint(training_args.output_dir)
            if checkpoint is None and os.listdir(training_args.output_dir):
                raise ValueError(
                    f"Output directory ({training_args.output_dir}) already exists and is not empty. "
                    "Use --overwrite_output_dir to continue."
                )
            if checkpoint is not None:
                logger.info("Resuming from checkpoint: %s", checkpoint)

        logger.info("*** Training ***")
        train_result = trainer.train(resume_from_checkpoint=checkpoint)
        trainer.save_model(training_args.output_dir)
        tokenizer.save_pretrained(training_args.output_dir)
        trainer.save_state()
        trainer.log_metrics("train", train_result.metrics)
        trainer.save_metrics("train", train_result.metrics)
        logger.info("Training finished.")

    if training_args.do_eval:
        logger.info("*** Evaluation ***")
        metrics = trainer.evaluate()
        try:
            metrics["perplexity"] = math.exp(metrics["eval_loss"])
        except OverflowError:
            metrics["perplexity"] = float("inf")
        trainer.log_metrics("eval", metrics)
        trainer.save_metrics("eval", metrics)
        logger.info("Evaluation finished. Perplexity: %.2f", metrics["perplexity"])


if __name__ == "__main__":
    main()
