# Bio Pipeline

This repository builds OriV/Rep training data from plasmid sequences. The upstream starting point is the output produced by OriV-Finder after querying plasmid FASTA files: `RIPs.csv` and `selected_ori_regions.csv`. The downstream pipeline then pairs OriV regions with Rep proteins, annotates Rep protein replication mechanisms, performs CD-HIT clustering, exports train/validation parquet datasets, builds tokenizer files, and provides BLASTn-based independent-test construction tools.

## Requirements

Python dependencies used by the pipeline include:

```powershell
python -m pip install pandas numpy biopython pyarrow pyhmmer tokenizers
```

External command-line tools:

- `cd-hit-est` for OriV clustering.
- `blastn` and `makeblastdb` for independent-test de-duplication.
- WSL is supported for running `cd-hit-est` on Windows.

## Quick Start

The full workflow starts from plasmid nucleotide FASTA files, runs OriV-Finder, and then runs this repository's data-building pipeline.

### Step 0: Predict OriV/Rep with OriV-Finder

Put plasmid FASTA files in `OriV-Finder/data/input`, or point the batch runner to an external directory such as `D:\dataprocess\PLSDB`.

Run OriV-Finder on a PLSDB FASTA file:

```powershell
cd C:\Users\MoSo\Desktop\Bio\OriV-Finder

powershell -ExecutionPolicy Bypass -File .\run_batch_orivfinder.ps1 ^
  -InputDir "D:\dataprocess\PLSDB" ^
  -Patterns "sequences.fasta"
```

For each FASTA input, OriV-Finder writes files like:

```text
OriV-Finder/data/output/<input-name>/
  RIPs.csv
  selected_ori_regions.csv
  All_IGSs.csv
  <input-name>.gbff
```

These `RIPs.csv` and `selected_ori_regions.csv` files are the initial files consumed by the main pipeline.

### Step 1: Build Model Data

Run the complete pipeline on an OriV-Finder output folder:

```powershell
cd C:\Users\MoSo\Desktop\Bio

python run_pipeline.py --root-dir "C:\Users\MoSo\Desktop\Bio\OriV-Finder\data\output" --cdhit-mode wsl
```

Write outputs to another directory:

```powershell
python run_pipeline.py ^
  --root-dir "C:\Users\MoSo\Desktop\Bio\OriV-Finder\data\output" ^
  --output-dir output-lab ^
  --cdhit-mode wsl
```

Include replication mechanism labels in the generated dataset and tokenizer:

```powershell
python run_pipeline.py ^
  --root-dir "C:\Users\MoSo\Desktop\Bio\OriV-Finder\data\output" ^
  --output-dir output-lab ^
  --cdhit-mode wsl ^
  --include-mechanism-column ^
  --mechanism-column replication_mechanism_term ^
  --include-mechanism-tokens
```

Skip CD-HIT if cluster output already exists:

```powershell
python run_pipeline.py ^
  --root-dir "C:\Users\MoSo\Desktop\Bio\OriV-Finder\data\output" ^
  --output-dir output-lab ^
  --cdhit-mode wsl ^
  --skip-cdhit
```

## Output Layout

For `--output-dir output-lab`, the main output folders are:

```text
OriV-Finder/
  data/output/         upstream OriV-Finder outputs: RIPs.csv, selected_ori_regions.csv

output-lab/
  data/                 intermediate CSV, FASTA, CD-HIT files
  train-data/           all_data.csv, train.parquet, validation.parquet
  tokenizer/            tokenizer.json, tokenizer_config.json, special_tokens_map.json
  logs/                 step logs
  visualizations/       species/length/CD-HIT figures
  independent-test/     independent test set built against your dataset
```

## Main Runner

`run_pipeline.py` is the end-to-end entry point. It calls the scripts in `pipeline_modules/scripts` in order.

Main options:

- `--root-dir`: raw input root folder containing RIP and selected_ori_regions CSV files.
- `--output-dir`: output root, default `output`.
- `--cdhit-mode native|wsl|mac|linux`: how to run `cd-hit-est`.
- `--cdhit-c`: CD-HIT identity threshold, default `0.9`.
- `--cdhit-dedup-mode keep-all|representative-only`: whether to keep all rows after clustering or only cluster representatives.
- `--validation-ratio`: validation ratio assigned by CD-HIT cluster, default `0.25`.
- `--include-mechanism-column`: include `replication_mechanism_term` in dataset files.
- `--include-mechanism-tokens` / `--no-include-mechanism-tokens`: control tokenizer mechanism tokens.
- `--skip-data-processing`, `--skip-cdhit`, `--skip-dataset-build`, `--skip-viz`: skip selected stages.

Examples:

```powershell
python run_pipeline.py --root-dir "D:\dataprocess\OriVresult" --cdhit-mode wsl
```

```powershell
python run_pipeline.py ^
  --root-dir "D:\dataprocess\OriVresult" ^
  --output-dir output-lab ^
  --cdhit-mode wsl ^
  --cdhit-c 0.9 ^
  --cdhit-n 5 ^
  --cdhit-t 8 ^
  --cdhit-dedup-mode keep-all
```

```powershell
python run_pipeline.py ^
  --root-dir "D:\dataprocess\OriVresult" ^
  --output-dir output-lab ^
  --skip-data-processing ^
  --skip-cdhit
```

## Pipeline Scripts

### `merge_inputs.py`

Merges OriV-Finder result CSV files under a root directory into two normalized files. This is the first script in the downstream data-building pipeline after OriV-Finder has finished predicting OriV/Rep features from plasmid FASTA files.

- RIP protein table: `RIPs.csv`
- Ori region table: `selected_ori_regions.csv`

Default run:

```powershell
python pipeline_modules\scripts\merge_inputs.py --root-dir "D:\dataprocess\OriVresult"
```

Custom output:

```powershell
python pipeline_modules\scripts\merge_inputs.py ^
  --root-dir "D:\dataprocess\OriVresult" ^
  --output-rip output-lab\data\RIPs.csv ^
  --output-ori output-lab\data\selected_ori_regions.csv
```

Append instead of replacing previous outputs:

```powershell
python pipeline_modules\scripts\merge_inputs.py ^
  --root-dir "D:\dataprocess\OriVresult" ^
  --append
```

### `pair_oriv_rep.py`

Pairs each OriV region with nearby Rep proteins and enriches host taxonomy. It consumes the two files created by `merge_inputs.py`.

Default run:

```powershell
python pipeline_modules\scripts\pair_oriv_rep.py
```

Custom input/output:

```powershell
python pipeline_modules\scripts\pair_oriv_rep.py ^
  --rip-csv output-lab\data\RIPs.csv ^
  --ori-csv output-lab\data\selected_ori_regions.csv ^
  --output output-lab\data\merged_final_optimized.csv ^
  --cache-file taxonomy_cache.pkl
```

Use a different Entrez email:

```powershell
python pipeline_modules\scripts\pair_oriv_rep.py ^
  --rip-csv output-lab\data\RIPs.csv ^
  --ori-csv output-lab\data\selected_ori_regions.csv ^
  --entrez-email your_email@example.com
```

### `annotate_replication_mechanism.py`

Annotates Rep proteins with `rep_plasmid_core_db.hmm` and maps Pfam hits to replication mechanism labels such as `RCR`, `theta`, or `unresolved`.

Default-style run:

```powershell
python pipeline_modules\scripts\annotate_replication_mechanism.py ^
  --input output-lab\data\merged_final_optimized.csv ^
  --output output-lab\data\merged_final_optimized_replication_annotated.csv
```

Specify HMM/rule resources and CPU count:

```powershell
python pipeline_modules\scripts\annotate_replication_mechanism.py ^
  --input output-lab\data\merged_final_optimized.csv ^
  --output output-lab\data\merged_final_optimized_replication_annotated.csv ^
  --hmm-db pipeline_modules\resources\replication_annotation\rep_plasmid_core_db.hmm ^
  --rules pipeline_modules\resources\replication_annotation\pfam_replication_mechanism_table.tsv ^
  --cpus 8
```

Use stricter HMM hit filtering:

```powershell
python pipeline_modules\scripts\annotate_replication_mechanism.py ^
  --input output-lab\data\merged_final_optimized.csv ^
  --output output-lab\data\merged_final_optimized_replication_annotated.csv ^
  --evalue-threshold 1e-5
```

### `build_cdhit_fasta.py`

Builds the OriV FASTA file used by `cd-hit-est`. It reads the annotated CSV and writes FASTA records using `ori_id`.

Basic run:

```powershell
python pipeline_modules\scripts\build_cdhit_fasta.py ^
  --input output-lab\data\merged_final_optimized_replication_annotated.csv ^
  --output output-lab\data\ori.fasta
```

Filter by OriV length:

```powershell
python pipeline_modules\scripts\build_cdhit_fasta.py ^
  --input output-lab\data\merged_final_optimized_replication_annotated.csv ^
  --output output-lab\data\ori.fasta ^
  --min-len 100 ^
  --max-len 20000
```

### Running CD-HIT

`run_pipeline.py` can run CD-HIT for you. To run it manually:

```powershell
wsl cd-hit-est -i /mnt/c/Users/MoSo/Desktop/Bio/output-lab/data/ori.fasta -o /mnt/c/Users/MoSo/Desktop/Bio/output-lab/data/cd_hit_process.fasta -c 0.9 -n 5 -M 0 -T 8 -d 0
```

On Linux/macOS or a shell where `cd-hit-est` is directly available:

```bash
cd-hit-est -i output-lab/data/ori.fasta -o output-lab/data/cd_hit_process.fasta -c 0.9 -n 5 -M 0 -T 8 -d 0
```

### `restore_cdhit_to_csv.py`

Reads CD-HIT `.clstr` output and merges cluster labels back into the annotated CSV. It assigns train/validation split by cluster, so similar OriV sequences do not cross split boundaries.

Keep all original rows and add cluster metadata:

```powershell
python pipeline_modules\scripts\restore_cdhit_to_csv.py ^
  --fasta output-lab\data\cd_hit_process.fasta ^
  --csv output-lab\data\merged_final_optimized_replication_annotated.csv ^
  --output output-lab\data\final_data.csv ^
  --dedup-mode keep-all
```

Keep only CD-HIT representatives:

```powershell
python pipeline_modules\scripts\restore_cdhit_to_csv.py ^
  --fasta output-lab\data\cd_hit_process.fasta ^
  --csv output-lab\data\merged_final_optimized_replication_annotated.csv ^
  --output output-lab\data\final_data.csv ^
  --dedup-mode representative-only
```

Change validation ratio:

```powershell
python pipeline_modules\scripts\restore_cdhit_to_csv.py ^
  --fasta output-lab\data\cd_hit_process.fasta ^
  --csv output-lab\data\merged_final_optimized_replication_annotated.csv ^
  --output output-lab\data\final_data.csv ^
  --validation-ratio 0.1 ^
  --seed 42
```

### `build_training_dataset.py`

Converts `final_data.csv` into the model dataset format:

- `all_data.csv`
- `train.parquet`
- `validation.parquet`

Basic run:

```powershell
python pipeline_modules\scripts\build_training_dataset.py ^
  --input output-lab\data\final_data.csv ^
  --output-dir output-lab\train-data
```

Include replication mechanism column:

```powershell
python pipeline_modules\scripts\build_training_dataset.py ^
  --input output-lab\data\final_data.csv ^
  --output-dir output-lab\train-data ^
  --include-mechanism-column ^
  --mechanism-column replication_mechanism_term
```

Fallback random split if no split column exists:

```powershell
python pipeline_modules\scripts\build_training_dataset.py ^
  --input output-lab\data\final_data.csv ^
  --output-dir output-lab\train-data ^
  --validation-ratio 0.25 ^
  --seed 42
```

### `build_tokenizer.py`

Builds tokenizer artifacts from `all_data.csv` or parquet data. It includes amino acid tokens, nucleotide tokens, selected species tokens, and optional mechanism tokens.

Basic run:

```powershell
python pipeline_modules\scripts\build_tokenizer.py ^
  --input output-lab\train-data\all_data.csv ^
  --tokenizer-json output-lab\tokenizer\tokenizer.json
```

Include mechanism tokens:

```powershell
python pipeline_modules\scripts\build_tokenizer.py ^
  --input output-lab\train-data\all_data.csv ^
  --tokenizer-json output-lab\tokenizer\tokenizer.json ^
  --include-mechanism-tokens
```

Write all tokenizer files explicitly:

```powershell
python pipeline_modules\scripts\build_tokenizer.py ^
  --input output-lab\train-data\all_data.csv ^
  --tokenizer-json output-lab\tokenizer\tokenizer.json ^
  --tokenizer-config-json output-lab\tokenizer\tokenizer_config.json ^
  --special-tokens-map-json output-lab\tokenizer\special_tokens_map.json ^
  --species-min-count 6
```

## Standalone Modules

`pipeline_modules/independent_test` is a standalone module at the same level as `pipeline_modules/visualization`. It is not called by `run_pipeline.py`; run these tools manually when constructing or auditing independent test sets.

### `build_independent_test_blastn.py`

Builds an independent test set from external parquet files by filtering against your dataset. This is useful when your dataset is the training reference and the external files are candidate test data.

Default reference:

```text
output-lab/train-data/all_data.csv
```

Default external candidates:

```text
test-00000-of-00001.parquet
train-00000-of-00001 (1).parquet
validation-00000-of-00001 (1).parquet
```

Run with defaults:

```powershell
python pipeline_modules\independent_test\build_independent_test_blastn.py
```

Keep BLAST intermediate files:

```powershell
python pipeline_modules\independent_test\build_independent_test_blastn.py --keep-blast-workdir
```

Change exact filtering and internal duplicate handling:

```powershell
python pipeline_modules\independent_test\build_independent_test_blastn.py ^
  --internal-dedup-key row ^
  --reference-exact-key oriv-plasmid
```

Use stricter BLASTn thresholds:

```powershell
python pipeline_modules\independent_test\build_independent_test_blastn.py ^
  --pident-threshold 97 ^
  --qcovs-threshold 97
```

Outputs:

```text
output-lab/independent-test/independent_test.csv
output-lab/independent-test/independent_test.parquet
output-lab/independent-test/dedup_report.txt
output-lab/independent-test/rejected_exact.csv
output-lab/independent-test/rejected_blastn.csv
output-lab/independent-test/blast_hits.tsv
```

### `build_external_model_independent_test.py`

Builds an independent test set for a model trained on the external `train + validation` parquet files. The training reference is the external train/validation data; the candidate independent-test pool is your dataset plus the external test parquet.

Default training reference:

```text
train-00000-of-00001 (1).parquet
validation-00000-of-00001 (1).parquet
```

Default candidates:

```text
output-lab/train-data/all_data.csv
test-00000-of-00001.parquet
```

Run with defaults:

```powershell
python pipeline_modules\independent_test\build_external_model_independent_test.py
```

Custom output directory:

```powershell
python pipeline_modules\independent_test\build_external_model_independent_test.py ^
  --output-dir output-lab\external-model-independent-test
```

Use only your dataset as candidate:

```powershell
python pipeline_modules\independent_test\build_external_model_independent_test.py ^
  --candidate-csv output-lab\train-data\all_data.csv ^
  --candidate-parquet
```

Outputs:

```text
output-lab/external-model-independent-test/independent_test.csv
output-lab/external-model-independent-test/independent_test.parquet
output-lab/external-model-independent-test/dedup_report.txt
```

### `summarize_external_oriv_duplicates.py`

Counts duplicated OriV sequences across external parquet files and exports duplicate details. It also writes the most repeated OriV and all corresponding rows.

Run with defaults:

```powershell
python pipeline_modules\independent_test\summarize_external_oriv_duplicates.py
```

Custom files:

```powershell
python pipeline_modules\independent_test\summarize_external_oriv_duplicates.py ^
  --external-parquet test-00000-of-00001.parquet "train-00000-of-00001 (1).parquet" "validation-00000-of-00001 (1).parquet" ^
  --output-dir output-lab\independent-test\oriv-duplicate-summary
```

Outputs:

```text
oriv_counts.csv
duplicated_oriv_counts.csv
duplicated_oriv_rows.csv
top_duplicated_oriv_rows.csv
top_duplicated_oriv_sequence.txt
```

## De-Duplication Modes

`--internal-dedup-key` controls duplicate handling inside the candidate pool before comparing to the reference:

- `none`: do not remove candidate duplicates.
- `row`: remove only identical standardized rows. This is the default and usually safest.
- `oriv`: remove rows with the same exact `oriv_sequence`; strict and may collapse biologically distinct records.
- `rep-oriv`: remove duplicate `rep_protein + oriv_sequence` pairs.
- `plasmid-oriv`: remove duplicate `plasmid_id + oriv_sequence` pairs.

`--reference-exact-key` controls exact filtering against the reference dataset before BLASTn:

- `oriv`: remove candidates whose `oriv_sequence` exactly appears in the reference.
- `oriv-plasmid`: remove exact `oriv_sequence` matches and same `plasmid_id` matches. This is the default.
- `strict`: also remove exact `rep_protein` matches and `rep_protein + oriv_sequence` matches.

BLASTn filtering then removes candidates whose OriV is highly similar to the reference:

```text
pident > 95
qcovs > 95
```

These defaults follow the OriGen-style independent-test filtering idea: exact duplicates are removed first, then highly similar OriV sequences are removed by BLASTn.

## Resource Files

Replication mechanism annotation resources are stored in:

```text
pipeline_modules/resources/replication_annotation/rep_plasmid_core_db.hmm
pipeline_modules/resources/replication_annotation/pfam_replication_mechanism_table.tsv
pipeline_modules/resources/replication_annotation/pfam_replication_mechanism_notes.md
```

## Useful Checks

Show help for the main runner:

```powershell
python run_pipeline.py --help
```

Show help for an individual script:

```powershell
python pipeline_modules\independent_test\build_independent_test_blastn.py --help
```

Check that BLAST is available:

```powershell
where.exe blastn
where.exe makeblastdb
```

Check that CD-HIT is available in WSL:

```powershell
wsl which cd-hit-est
```


