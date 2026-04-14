#!/usr/bin/env python
# coding=utf-8
"""
OriGen training script.

This version is refactored for clarity and adds an explicit switch for whether
to include the replication-mechanism token in each training example.

Input format expected by the training parquet:
    - host species column, e.g. `species` or `host_species`
    - Rep protein column, e.g. `rep_seq` or `rep_protein`
    - oriV DNA sequence column, e.g. `OriC sequence`, `oriv_sequence`, or `rep_dna_seq`
    - optional replication mechanism column, e.g. `replication_mechanism_term`

Examples:
    python 10-train_origen.py `
      --tokenizer-path tokenizer `
      --dataset-path training_data `
      --species-col species `
      --rep-col rep_seq `
      --seq-col "OriC sequence" `
      --include-mechanism-token `
      --mechanism-col replication_mechanism_term

    python 10-train_origen.py `
      --tokenizer-path tokenizer `
      --dataset-path training_data `
      --species-col species `
      --rep-col rep_seq `
      --seq-col "rep_dna_seq"
"""

from __future__ import annotations

import logging
import math
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import torch
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

    species_col: str = field(default="species", metadata={"help": "Species column name"})
    rep_col: str = field(default="rep_seq", metadata={"help": "Rep protein column name"})
    seq_col: str = field(default="OriC sequence", metadata={"help": "oriV DNA sequence column name"})

    include_mechanism_token: bool = field(
        default=False,
        metadata={"help": "Append the replication mechanism token to each example"},
    )
    mechanism_col: str = field(
        default="replication_mechanism_term",
        metadata={"help": "Replication-mechanism column name"},
    )

    model_size: int = field(default=512, metadata={"help": "Hidden size"})
    num_hidden_layers: int = field(default=8, metadata={"help": "Number of transformer layers"})
    num_attention_heads: int = field(default=8, metadata={"help": "Number of attention heads"})
    max_seq_len: int = field(default=1024, metadata={"help": "Maximum sequence length"})

    do_lower_dna: bool = field(default=True, metadata={"help": "Lower-case DNA sequence before tokenization"})
    validation_split_ratio: float = field(
        default=0.05,
        metadata={"help": "Validation split ratio used when no validation set is provided"},
    )
    use_windows: bool = field(
        default=False,
        metadata={"help": "Use sliding windows for sequences longer than max_seq_len"},
    )


@dataclass
class CustomTrainingArguments(TrainingArguments):
    """Standard Hugging Face training arguments."""


def normalize_text(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() == "nan":
        return ""
    return text


def normalize_bracketed_token(value: object) -> str:
    text = normalize_text(value)
    if not text:
        return ""
    text = text.replace("[", "").replace("]", "")
    text = text.replace("_", " ")
    text = re.sub(r"\s+", " ", text).strip()
    if not text:
        return ""
    return f"[{text}]"


def normalize_species_token(value: object) -> str:
    return normalize_bracketed_token(value)


def resolve_column_name(column_names: List[str], preferred: str, aliases: List[str]) -> str:
    if preferred in column_names:
        return preferred
    for alias in aliases:
        if alias in column_names:
            return alias
    return preferred


def apply_lightweight_training_defaults(training_args: TrainingArguments) -> None:
    """Use a smaller first-run profile unless the user explicitly overrode it."""

    argv = set(sys.argv)

    def not_overridden(*names: str) -> bool:
        return all(name not in argv for name in names)

    if not_overridden("--per_device_train_batch_size", "--per-device-train-batch-size"):
        training_args.per_device_train_batch_size = 2
    if not_overridden("--per_device_eval_batch_size", "--per-device-eval-batch-size"):
        training_args.per_device_eval_batch_size = 2
    if not_overridden("--num_train_epochs", "--num-train-epochs"):
        training_args.num_train_epochs = 1.0
    if not_overridden("--gradient_accumulation_steps", "--gradient-accumulation-steps"):
        training_args.gradient_accumulation_steps = 1


def apply_gpu_defaults(model_args: ModelDataArguments, training_args: TrainingArguments) -> None:
    """Tune the run for a small CUDA GPU if the user did not override key knobs."""

    if not torch.cuda.is_available():
        logger.info("CUDA is not available; keeping CPU-safe defaults.")
        return

    props = torch.cuda.get_device_properties(0)
    vram_gb = props.total_memory / (1024**3)
    logger.info("CUDA device detected: %s (%.2f GiB)", props.name, vram_gb)

    if vram_gb < 6:
        model_size = 384
        hidden_layers = 6
        attention_heads = 6
        max_seq_len = 768
        train_batch_size = 1
        eval_batch_size = 1
        grad_accum = 16
        epochs = 3.0
    elif vram_gb < 10:
        model_size = 512
        hidden_layers = 8
        attention_heads = 8
        max_seq_len = 1024
        train_batch_size = 1
        eval_batch_size = 1
        grad_accum = 8
        epochs = 3.0
    elif vram_gb < 18:
        model_size = 768
        hidden_layers = 12
        attention_heads = 12
        max_seq_len = 1536
        train_batch_size = 2
        eval_batch_size = 2
        grad_accum = 4
        epochs = 4.0
    else:
        model_size = 1024
        hidden_layers = 16
        attention_heads = 16
        max_seq_len = 2048
        train_batch_size = 4
        eval_batch_size = 4
        grad_accum = 2
        epochs = 5.0

    argv = set(sys.argv)

    def not_overridden(*names: str) -> bool:
        return all(name not in argv for name in names)

    if not_overridden("--model_size", "--model-size"):
        model_args.model_size = model_size
    if not_overridden("--num_hidden_layers", "--num-hidden-layers"):
        model_args.num_hidden_layers = hidden_layers
    if not_overridden("--num_attention_heads", "--num-attention-heads"):
        model_args.num_attention_heads = attention_heads
    if not_overridden("--max_seq_len", "--max-seq-len"):
        model_args.max_seq_len = max_seq_len

    if not_overridden("--per_device_train_batch_size", "--per-device-train-batch-size"):
        training_args.per_device_train_batch_size = train_batch_size
    if not_overridden("--per_device_eval_batch_size", "--per-device-eval-batch-size"):
        training_args.per_device_eval_batch_size = eval_batch_size
    if not_overridden("--gradient_accumulation_steps", "--gradient-accumulation-steps"):
        training_args.gradient_accumulation_steps = grad_accum
    if not_overridden("--num_train_epochs", "--num-train-epochs", "--max_steps", "--max-steps"):
        training_args.num_train_epochs = epochs
    if not_overridden("--fp16", "--no_fp16", "--bf16", "--no_bf16"):
        training_args.fp16 = True
        training_args.bf16 = False
    if not_overridden("--gradient_checkpointing", "--no_gradient_checkpointing"):
        training_args.gradient_checkpointing = True
    if not_overridden("--dataloader_num_workers", "--dataloader-num-workers"):
        training_args.dataloader_num_workers = 0
    if not_overridden("--save_strategy", "--save-strategy"):
        training_args.save_strategy = "epoch"
    if not_overridden("--eval_strategy", "--evaluation_strategy", "--eval-strategy", "--evaluation-strategy"):
        training_args.eval_strategy = "epoch"
    if not_overridden("--logging_steps", "--logging-steps"):
        training_args.logging_steps = 50
    if not_overridden("--save_total_limit", "--save-total-limit"):
        training_args.save_total_limit = 2
    if not_overridden("--report_to", "--report-to"):
        training_args.report_to = "none"


def enable_train_and_eval_by_default(training_args: TrainingArguments) -> None:
    """Turn on training/eval unless the user explicitly opted out."""

    argv = set(sys.argv)

    if "--no_do_train" not in argv and "--do_train" not in argv:
        training_args.do_train = True
    if "--no_do_eval" not in argv and "--do_eval" not in argv:
        training_args.do_eval = True


def format_example(example: Dict[str, Any], args: ModelDataArguments) -> str:
    species = normalize_species_token(example.get(args.species_col, ""))
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
        mech = normalize_bracketed_token(example.get(args.mechanism_col, ""))
        if not mech:
            mech = "[unresolved]"
        parts.append(mech)

    parts.append(seq)
    return " ".join(parts)


def build_text_sequence(example: Dict[str, Any], args: ModelDataArguments) -> str:
    return format_example(example, args)


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
            texts.append(build_text_sequence(example, args))
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

    bos_id = tokenizer.bos_token_id
    eos_id = tokenizer.eos_token_id
    if bos_id is None or eos_id is None:
        raise ValueError("Tokenizer must define bos_token and eos_token.")

    content_max_len = max(1, model_args.max_seq_len - 2)
    overlap = content_max_len // 2 if model_args.use_windows else 0
    step = max(1, content_max_len - overlap)

    def tokenize_batch(examples: Dict[str, List[str]]) -> Dict[str, Any]:
        input_ids: List[List[int]] = []
        attention_mask: List[List[int]] = []

        for text in examples["text"]:
            if not text:
                continue

            encoding = tokenizer(
                text,
                truncation=False,
                padding=False,
                add_special_tokens=False,
            )
            ids = encoding["input_ids"]
            if not ids:
                continue

            if model_args.use_windows and len(ids) > content_max_len:
                start = 0
                while start < len(ids):
                    chunk = ids[start : start + content_max_len]
                    if not chunk:
                        break
                    window_ids = [bos_id] + chunk + [eos_id]
                    input_ids.append(window_ids)
                    attention_mask.append([1] * len(window_ids))
                    if start + content_max_len >= len(ids):
                        break
                    start += step
            else:
                chunk = ids[:content_max_len]
                window_ids = [bos_id] + chunk + [eos_id]
                input_ids.append(window_ids)
                attention_mask.append([1] * len(window_ids))

        return {
            "input_ids": input_ids,
            "attention_mask": attention_mask,
            "labels": [list(ids) for ids in input_ids],
        }

    tokenized = non_empty.map(
        tokenize_batch,
        batched=True,
        remove_columns=non_empty["train"].column_names,
        desc="Tokenizing",
    )
    return tokenized


def build_model(tokenizer, model_args: ModelDataArguments, training_args: TrainingArguments):
    config = RoFormerConfig()
    config.vocab_size = len(tokenizer)
    config.hidden_size = model_args.model_size
    config.num_hidden_layers = model_args.num_hidden_layers
    config.num_attention_heads = model_args.num_attention_heads
    config.intermediate_size = model_args.model_size * 4
    config.hidden_act = "gelu"
    config.attention_probs_dropout_prob = 0.1
    config.hidden_dropout_prob = 0.1
    config.max_position_embeddings = model_args.max_seq_len + 512
    config.type_vocab_size = 2
    config.is_decoder = True
    config.pad_token_id = tokenizer.pad_token_id
    config.eos_token_id = tokenizer.eos_token_id
    config.bos_token_id = tokenizer.bos_token_id
    config.rotary_value = True
    config.use_cache = not training_args.gradient_checkpointing

    model = AutoModelForCausalLM.from_config(config)
    if hasattr(model, "resize_token_embeddings"):
        model.resize_token_embeddings(len(tokenizer))
    if training_args.gradient_checkpointing and hasattr(model, "gradient_checkpointing_enable"):
        model.gradient_checkpointing_enable()
    return model


def main() -> None:
    parser = HfArgumentParser((ModelDataArguments, CustomTrainingArguments))
    if len(sys.argv) == 2 and sys.argv[1].endswith(".json"):
        model_args, training_args = parser.parse_json_file(json_file=os.path.abspath(sys.argv[1]))
    else:
        model_args, training_args = parser.parse_args_into_dataclasses()

    if not os.path.exists(model_args.tokenizer_path):
        raise FileNotFoundError(f"Tokenizer path does not exist: {model_args.tokenizer_path}")

    apply_gpu_defaults(model_args, training_args)
    enable_train_and_eval_by_default(training_args)
    set_seed(training_args.seed)
    if not torch.cuda.is_available():
        apply_lightweight_training_defaults(training_args)

    if torch.cuda.is_available():
        torch.backends.cuda.matmul.allow_tf32 = True

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

    column_names = raw_datasets["train"].column_names
    model_args.species_col = resolve_column_name(column_names, model_args.species_col, ["host_species"])
    model_args.rep_col = resolve_column_name(column_names, model_args.rep_col, ["rep_protein"])
    model_args.seq_col = resolve_column_name(
        column_names,
        model_args.seq_col,
        ["oriv_sequence", "rep_dna_seq"],
    )
    if model_args.include_mechanism_token:
        model_args.mechanism_col = resolve_column_name(
            column_names,
            model_args.mechanism_col,
            ["replication_mechanism_call"],
        )

    logger.info(
        "Resolved columns - species: %s, rep: %s, seq: %s%s",
        model_args.species_col,
        model_args.rep_col,
        model_args.seq_col,
        f", mechanism: {model_args.mechanism_col}" if model_args.include_mechanism_token else "",
    )

    with_text = add_text_column(raw_datasets, model_args)
    tokenized_datasets = tokenize_datasets(with_text, tokenizer, model_args)

    logger.info("Initializing RoFormer causal LM")
    model = build_model(tokenizer, model_args, training_args)
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
        processing_class=tokenizer,
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
