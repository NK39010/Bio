# Compare Plasmid IDs

This script compares `plasmid_id` values between two files.

It supports:

- `.parquet`
- `.csv`

It also auto-detects common column name variants:

- `plasmid_id`
- `plsmid_id`
- `plasmid`

Example:

```bash
python relationships/compare_plasmid_ids.py \
  --left "D:/dataprocess/数据示例/train-00000-of-00001.parquet" \
  --right "C:/Users/MoSo/Desktop/Bio/final_data.parquet" \
  --output-dir relationships/compare_output
```
