# Bio Pipeline

## 默认行为

- `CD-HIT 返还` 默认 `--cdhit-dedup-mode keep-all`（只聚类标签，不按代表序列去重）。
- `数据集` 默认不加入复制机制列（不传 `--include-mechanism-column`）。
- `tokenizer` 默认不加入复制机制 token（默认等价于 `--no-include-mechanism-tokens`）。

只有显式传参时才会加入复制机制相关内容。

## 关键参数

- `--cdhit-dedup-mode keep-all|representative-only`
- `--include-mechanism-column`
- `--mechanism-column replication_mechanism_term`
- `--include-mechanism-tokens`（默认不加）

## 调用示例

### Windows + WSL（默认：不加复制机制）

```bash
python run_pipeline.py --root-dir "D:\dataprocess\OriVresult" --cdhit-mode wsl
```

### 显式启用复制机制列和机制 token

```bash
python run_pipeline.py ^
  --root-dir "D:\dataprocess\OriVresult" ^
  --cdhit-mode wsl ^
  --include-mechanism-column ^
  --mechanism-column replication_mechanism_term ^
  --include-mechanism-tokens
```

### 不跑 CD-HIT，且不加复制机制

```bash
python run_pipeline.py ^
  --root-dir "D:\dataprocess\OriVresult" ^
  --cdhit-mode wsl ^
  --skip-cdhit ^
  --cdhit-dedup-mode keep-all
```

## 数据库资源位置

复制机制注释数据库位于：

- `pipeline_modules/resources/replication_annotation/rep_plasmid_core_db.hmm`
- `pipeline_modules/resources/replication_annotation/pfam_replication_mechanism_table.tsv`
