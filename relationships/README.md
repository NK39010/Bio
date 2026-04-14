# ID Relationship Explorer

This folder stores scripts and outputs for exploring the relationship between:

- `plasmid_id`
- `ori_id`
- `rep_id`

## What it generates

- summary tables for each ID
- pair-count tables
- triplet-count table
- PNG bar charts for the most connected IDs
- a Markdown report for quick inspection

## Example

```bash
python relationships/show_id_relations.py --input final_data.parquet --output-dir relationships/output
```

CSV is also supported:

```bash
python relationships/show_id_relations.py --input final_data.csv --output-dir relationships/output_csv
```
