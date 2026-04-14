# Train Workspace

This folder keeps training-related files together.

## Layout

- `train/train_origen.py`: training entrypoint
- `train/outputs/`: model checkpoints and trainer outputs

## Example

```bash
python train/train_origen.py \
  --tokenizer_path tokenizer \
  --dataset_path training_data \
  --output_dir train/outputs/origen_run1 \
  --species_col species \
  --rep_col rep_seq \
  --seq_col "OriC sequence" \
  --include_mechanism_token \
  --mechanism_col replication_mechanism_term
```

If you do not want the replication-mechanism token:

```bash
python train/train_origen.py \
  --tokenizer_path tokenizer \
  --dataset_path training_data \
  --output_dir train/outputs/origen_run2 \
  --species_col species \
  --rep_col rep_seq \
  --seq_col "rep_dna_seq"
```
