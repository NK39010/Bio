# Pfam 家族到质粒复制方式的判定规则

表格文件：`dataprocess/Pfam/pfam_replication_mechanism_table.tsv`

## 1. 推荐判定顺序

1. 先看 `use_for_call = yes` 且 `confidence = high`
2. 若同一蛋白或同一质粒上出现 `Rep3_N + Rep3_C`、`RepA_N + Bac_RepA_C`、`Rep_OBD + RepB_C`、`Replitron_HUH + Replitron_C`、`RepD-like_N + Rep_trans` 这类成对结构域，优先采用组合结论
3. 若只有 `support_marker`，输出 `RCR_like`、`theta_like` 或 `uncertain`，不要硬判
4. 若只命中 `not_mechanism` 或 `use_for_call = no` 的家族，输出 `unresolved`

## 2. 可直接用于机制判定的经验规则

- 命中 `Rep_1`、`RepL`、`Rol_Rep_N`、`RepD-like_N`、`Replitron_HUH` 中任一个高可信家族：
  判为 `RCR`

- 命中 `Rep3_N/Rep3_C`、`RepA_N/Bac_RepA_C/WHD_RepA_N/RepA_C`、`Rep_OBD/RepB_C`、`RepC`、`TrfA`：
  判为 `theta`

- 只命中 `RepB-RCR_reg`：
  输出 `RCR_support`，最好再找对应 initiator 主结构域

- 只命中 `RepA1_leader`：
  输出 `theta_support`，说明可能是 antisense RNA 调控的 RepA 系统

- 只命中 `Rep_N`、`PaRep2a`、`PaRep2b`、`RepB_primase`：
  输出 `replication-associated, mechanism unresolved`

## 3. 不要直接拿来判复制方式的家族

以下几类建议全部排除出“复制方式判定”步骤：

- 染色体复制/解旋相关：
  `DnaA*`、`DnaB*`、`DnaG*`、`Helicase*`

- 古菌/真核复制相关：
  `MCM*`、`ORC*`、`Cdc6*`、`Rep-A_N`、`Rep_fac-*`

- 分配或稳定维持相关：
  `ParB*`、`RepB`

- 明显非细菌质粒相关：
  `Phage_*`、`Viral_Rep`、`Replicase`、`CNV-Replicase_N`、`Rep_1B`、`Rep_4`

## 4. 一个更稳妥的最终输出字段

建议你的脚本最后不要只输出二分类，而是输出下面四类之一：

- `RCR`
- `theta`
- `RCR_like/theta_like`
- `unresolved`

## 5. 这张表的边界

- 这是“Pfam 结构域 -> 复制方式”的经验映射，不等同于实验验证
- `ColE1-like` 这类 RNA 主导复制子常常没有典型 Rep 蛋白，所以单靠这个表会漏检
- 若一条质粒同时有两套 replicon 信号，要考虑多复制子质粒，而不是强制只保留一种机制
