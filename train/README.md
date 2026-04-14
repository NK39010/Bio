# Train Workspace

This folder keeps training-related files together.

## Layout

- `train/train_origen.py`: training entrypoint
- `train/outputs/`: model checkpoints and trainer outputs

## Example

```powershell
python train/train_origen.py `
  --tokenizer_path tokenizer `
  --dataset_path training_data `
  --output_dir train/outputs/origen_run1 `
  --species_col species `
  --rep_col rep_seq `
  --seq_col "OriC sequence" `
  --include_mechanism_token `
  --mechanism_col replication_mechanism_term
```

On a CUDA GPU, the script auto-tunes by VRAM tier unless you override it:

- under 6 GB: `model_size=384`, `max_seq_len=768`
- 6 to 10 GB: `model_size=512`, `max_seq_len=1024`
- 10 to 18 GB: `model_size=768`, `max_seq_len=1536`
- 18 GB and up: `model_size=1024`, `max_seq_len=2048`

Batch size, gradient accumulation, and epoch count also scale upward with more VRAM.
`fp16=True` and `gradient_checkpointing=True` stay on by default unless you change them.

If CUDA is not available, it falls back to a CPU-safe lightweight profile.

Training and evaluation are enabled by default unless you explicitly opt out.

If you do not want the replication-mechanism token:

```powershell
python train/train_origen.py `
  --tokenizer_path tokenizer `
  --dataset_path training_data `
  --output_dir train/outputs/origen_run2 `
  --species_col species `
  --rep_col rep_seq `
  --seq_col "rep_dna_seq"
```
