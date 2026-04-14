# Bio Pipeline (00-09)

This repository builds a full OriC/Rep analysis and training-data pipeline from raw CSV inputs.

## Pipeline Overview

Recommended end-to-end order:

1. `00-合并.py`
   - Reads the raw directory recursively
   - Outputs:
     - `data/RIPs.csv`
     - `data/selected_ori_regions.csv`

2. `01-OriV&RIFs-mix.py`
   - Ori-centric pairing:
     - each OriC keeps the nearest Rep(s) on the same plasmid
     - OriCs without any matching Rep are kept as no-rep rows
   - Output:
     - `merged_final_optimized.csv`

3. `02-种属分布.py` and `03-长度分布.py` before CD-HIT
   - Outputs:
     - `species_before.xlsx`
     - `species_before.png`
     - `length_before.png`
   - `02` now uses raw rows directly and no longer deduplicates by `ori_id`

4. `02-annotate_replication_mechanism.py`
   - Adds replication-mechanism annotations
   - Output:
     - `merged_final_optimized_replication_annotated.csv`

5. `04-cd-hit序列构建.py`
   - Builds the CD-HIT FASTA input
   - Output:
     - `ori.fasta`

6. `05` via `run_pipeline.py`
   - Runs `cd-hit-est`
   - Outputs:
     - `cd_hit_process.fasta`
     - `cd_hit_process.fasta.clstr`

7. `06-cd-hit去冗余结果可视化.py`
   - Output:
     - `cdhit_redundancy.png`

8. `07-cd-hit序列返还csv构建.py`
   - Maps CD-HIT clusters back to the annotated CSV
   - Default validation split is cluster-aware and close to `9:1`
   - Output:
     - `final_data.csv`

9. `02-种属分布.py` after CD-HIT
   - Outputs:
     - `species_after.xlsx`
     - `species_after.png`

10. `03-长度分布.py` after CD-HIT
    - Output:
      - `length_after.png`

11. `08-构建训练用的parquet.py`
    - Splits the final dataset into train / validation files
    - Outputs:
      - `training_data/train.csv`
      - `training_data/validation.csv`
      - `training_data/train.parquet`
      - `training_data/validation.parquet`

12. `09-构建词元.py`
    - Builds tokenizer files from the training CSV
    - Default input:
      - `training_data/train.csv`
    - Outputs:
      - `tokenizer/tokenizer.json`
      - `tokenizer/tokenizer_config.json`
      - `tokenizer/special_tokens_map.json`

## One-Command Run

Windows / WSL example:

```bash
python run_pipeline.py --root-dir "D:\dataprocess\OriVresult" --cdhit-mode wsl
```

If you want to rerun only the later steps with an existing `final_data.csv`:

```bash
python run_pipeline.py --root-dir "D:\dataprocess\OriVresult" --cdhit-mode wsl --skip-00 --skip-01 --skip-02 --skip-03 --skip-annotate --skip-04 --skip-05 --skip-06 --skip-07
```

## Key Arguments

`run_pipeline.py` common arguments:

- `--root-dir`: root folder for step 00 input CSVs
- `--cdhit-mode {mac|wsl}`: how to run `cd-hit-est`
- `--cdhit-c`: CD-HIT identity threshold, default `0.9`
- `--cdhit-n`: word size, default `5`
- `--cdhit-t`: threads, default `8`
- `--final-csv`: merged post-CD-HIT CSV, default `final_data.csv`
- `--training-dir`: split train/validation output directory, default `training_data`
- `--token-json`: tokenizer JSON path, default `tokenizer/tokenizer.json`

Skipping steps:

- `--skip-00` to `--skip-09`

## Training Outputs

`08-构建训练用的parquet.py` now writes both CSV and Parquet files for each split.

Current split is approximately `9:1`:

- train: about 90%
- validation: about 10%

The split is cluster-aware to reduce leakage across train and validation.

## Tokenizer

`09-构建词元.py` supports CSV or Parquet input, but the pipeline now points it to `training_data/train.csv`.

It keeps:

- 20 amino-acid tokens
- 4 nucleotide tokens
- 3 replication-mechanism tokens
- species names that appear more than 5 times

## Training Defaults

`10-train_origen.py` now defaults to the columns produced by `training_data`:

- `species_col=species`
- `rep_col=rep_seq`
- `seq_col="OriC sequence"`

If you want to train on the Rep DNA sequence instead, pass:

- `--seq_col "rep_dna_seq"`

The replication mechanism is optional:

- add it with `--include_mechanism_token --mechanism_col replication_mechanism_term`
- omit those flags to train on the three base fields only

On CUDA, the script auto-tunes by VRAM tier unless you override it:

- under 6 GB: `model_size=384`, `max_seq_len=768`
- 6 to 10 GB: `model_size=512`, `max_seq_len=1024`
- 10 to 18 GB: `model_size=768`, `max_seq_len=1536`
- 18 GB and up: `model_size=1024`, `max_seq_len=2048`

It also scales batch size, gradient accumulation, and epochs upward with more VRAM.
`fp16=True` and `gradient_checkpointing=True` stay on by default unless you change them.

For backward compatibility, it will also fall back to:

- `host_species` for species
- `rep_protein` for Rep
- `oriv_sequence` or `rep_dna_seq` for the DNA sequence

## Dependencies

Suggested environment:

- Python 3.10+
- `pandas`
- `numpy`
- `biopython`
- `matplotlib`
- `seaborn`
- `openpyxl`
- `pyarrow`
- `pyhmmer`
- `tqdm`

System dependency:

- `cd-hit-est`

## Notes

- `01` is Ori-centric by design now.
- `02` reports species distribution from raw rows, not deduplicated OriC records.
- `07` keeps a validation split close to `9:1` while preserving cluster integrity.
- `08` is responsible for splitting the final CSV into train / validation CSV and Parquet files.
