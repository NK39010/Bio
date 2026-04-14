# Train Workspace

This folder keeps training-related files together.

## Layout

- `train/train_origen.py`: training entrypoint
- `train/outputs/`: model checkpoints and trainer outputs

## Example

```bash
python train/train_origen.py \
  --tokenizer_path tokenizer \
  --dataset_path final_data.parquet \
  --output_dir train/outputs/origen_run1 \
  --species_col host_species \
  --rep_col rep_protein \
  --seq_col oriv_sequence \
  --include_mechanism_token \
  --mechanism_col replication_mechanism_term
```

If you do not want the replication-mechanism token:

```bash
python train/train_origen.py \
  --tokenizer_path tokenizer \
  --dataset_path final_data.parquet \
  --output_dir train/outputs/origen_run2 \
  --species_col host_species \
  --rep_col rep_protein \
  --seq_col oriv_sequence
```
