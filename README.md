# Bio Pipeline (00-09)

这个仓库用于构建一条从原始 CSV 到去冗余、注释、可视化、训练数据导出的全流程管线。

## 流程总览

主流程按以下步骤执行（可跳步）：

1. `00-合并.py`  
   从初始数据目录递归读取 CSV，合并为：
   - `data/RIPs.csv`
   - `data/selected_ori_regions.csv`

2. `01-OriV&RIFs-mix.py`  
   合并 Ori 与 RIP，产出：
   - `merged_final_optimized.csv`

3. `02-种属分布.py`（CD-HIT 前）  
   产出种属分布统计：
   - `species_before.xlsx`
   - `species_before.png`

4. `03-长度分布.py`（CD-HIT 前）  
   产出长度分布图：
   - `length_before.png`

5. `02-annotate_replication_mechanism.py`  
   注释复制方式，产出：
   - `merged_final_optimized_replication_annotated.csv`

6. `04-cd-hit序列构建.py`  
   生成 CD-HIT 输入 FASTA：
   - `ori.fasta`

7. `05`（通过 `run_pipeline.py` 内置命令执行）  
   执行 `cd-hit-est` 去冗余，产出：
   - `cd_hit_process.fasta`
   - `cd_hit_process.fasta.clstr`

8. `06-cd-hit去冗余结果可视化.py`  
   产出去冗余可视化：
   - `cdhit_redundancy.png`

9. `07-cd-hit序列返还csv构建.py`  
   将去冗余 FASTA 映射回注释表，产出：
   - `final_data.csv`

10. `02-种属分布.py`（CD-HIT 后）  
    - `species_after.xlsx`
    - `species_after.png`

11. `03-长度分布.py`（CD-HIT 后）  
    - `length_after.png`

12. `08-构建训练用的parquet.py`  
    构建训练用 parquet（包含复制方式词条列 `replication_mechanism_term`）：
    - `final_data.parquet`

13. `09-构建词元.py`  
    生成 HF 风格 tokenizer 文件到 `tokenizer/` 目录：
    - `tokenizer/tokenizer.json`
    - `tokenizer/tokenizer_config.json`
    - `tokenizer/special_tokens_map.json`
    - 20 种氨基酸
    - 4 种碱基
    - 3 种复制方式
    - 出现次数大于 5 次的物种名

---

## 一键运行（推荐）

### mac 直接运行 CD-HIT

```bash
python3 run_pipeline.py --root-dir "/你的初始数据文件夹" --cdhit-mode mac
```

### 通过 WSL 运行 CD-HIT

```bash
python3 run_pipeline.py --root-dir "/你的初始数据文件夹" --cdhit-mode wsl
```
python run_pipeline.py --root-dir "C:/Users/MoSo/Desktop/ori" --cdhit-mode wsl
---

## 关键参数

`run_pipeline.py` 常用参数：

- `--root-dir`：原始数据目录（必填）
- `--cdhit-mode {mac|wsl}`：CD-HIT 执行方式
- `--cdhit-c`：相似性阈值（默认 `0.9`）
- `--cdhit-n`：word size（默认 `5`）
- `--cdhit-t`：线程数（默认 `8`）
- `--final-csv`：去冗余返还后的 CSV（默认 `final_data.csv`）
- `--final-parquet`：08 的 parquet 输出（默认 `final_data.parquet`）
- `--token-json`：09 输出的 `tokenizer.json` 路径（默认 `tokenizer/tokenizer.json`，其余两个文件会自动同目录生成）

跳步参数：

- `--skip-00` ~ `--skip-09`：跳过对应步骤

示例（只重跑 08 和 09）：

```bash
python3 run_pipeline.py \
  --root-dir "/你的初始数据文件夹" \
  --skip-00 --skip-01 --skip-02 --skip-03 --skip-annotate --skip-04 --skip-05 --skip-06 --skip-07
```

---

## 08 / 09 职责说明

### `08-构建训练用的parquet.py`

- 输入：`final_data.csv`（或你指定的 CSV）
- 输出：`final_data.parquet`
- 会新增列：
  - `replication_mechanism_term`
- 该列来源：
  - 优先读取 `replication_mechanism_call`（来自 `02-annotate` 结果）
  - 生成格式：`RCR` / `theta` / `unresolved`

### `09-构建词元.py`

- 输入：`final_data.csv` 或 `final_data.parquet`
- 输出：
  - `tokenizer/tokenizer.json`
  - `tokenizer/tokenizer_config.json`
  - `tokenizer/special_tokens_map.json`
- 只保留 20 种氨基酸、4 种碱基、3 种复制方式，以及出现次数大于 5 次的物种名

---

## 单步运行示例

只跑 parquet：

```bash
python3 08-构建训练用的parquet.py --input final_data.csv --parquet final_data.parquet
```

只跑词元：

```bash
python3 09-构建词元.py --input final_data.parquet --tokenizer-json tokenizer/tokenizer.json
```

---

## 依赖建议

建议 Python 3.10+，并安装（按脚本实际用到的包）：

- `pandas`
- `numpy`
- `biopython`
- `matplotlib`
- `seaborn`
- `openpyxl`
- `pyarrow`
- `pyhmmer`
- `tqdm`

以及系统命令：

- `cd-hit-est`（`--cdhit-mode mac` 时本机可执行）
- 或 WSL 中可执行 `cd-hit-est`（`--cdhit-mode wsl`）
