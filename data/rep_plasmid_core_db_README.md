# rep_plasmid_core_db

这是从原始 `dataprocess/rep_mini_db.hmm` 重新筛过后的“质粒复制方式核心 Pfam 库”。

## 生成结果

- 保留模型列表：`dataprocess/Pfam/rep_plasmid_core_models.list`
- 删除模型列表：`dataprocess/Pfam/rep_plasmid_removed_models.list`
- 新核心数据库：`dataprocess/Pfam/rep_plasmid_core_db.hmm`
- 判定映射表：`dataprocess/Pfam/pfam_replication_mechanism_table.tsv`

## 筛选原则

- 保留：可直接支持 `RCR`、`theta` 或 `theta/RCR-like` 判定的质粒 Rep 核心家族
- 删除：染色体复制蛋白、ParB 分配系统、古菌/真核 ORC-MCM 相关、噬菌体/病毒 replicase、明显非细菌质粒家族

## 当前保留数量

- 原始模型数：131
- 保留模型数：23
- 删除模型数：108

## 当前保留的核心模型

`Rep_1`, `RepL`, `Rol_Rep_N`, `RepD-like_N`, `Rep_trans`, `RepB-RCR_reg`, `Replitron_HUH`, `Replitron_C`, `Rep3_N`, `Rep3_C`, `RepA_N`, `Bac_RepA_C`, `WHD_RepA_N`, `RepA_C`, `Rep_OBD`, `RepB_C`, `RepC`, `TrfA`, `Rep_assoc_RepA4`, `L_lactis_RepB_C`, `RepB_primase`, `RepB_primase_C`, `RepA1_leader`

## 说明

- 这个库更适合“复制方式判定”，不再追求把所有 replication-associated 模型都收进来
- 它仍然会漏掉 `ColE1-like` 这类以 RNA 为主、缺少典型 Rep 蛋白的复制子
- 如果后面你想做“更宽松的候选发现库”，可以再基于 `conditional` 家族额外做一个扩展版
