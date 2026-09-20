# 2026-09-20 coherence 列名规范化 + 覆盖层补齐 + 3a 序号

## 背景

`coherence_metrics` 表的列名是一路演进来的，积了两类问题：

1. **说错了的**：`centroid_x`/`centroid_y` 装的是 `geometric_median()`（Weiszfeld
   几何中位数），而 centroid 是均值——名字指错了估计量；`term_size` /
   `n_group_genes` 相邻两列都是整数、都带 size 语义，差一层 DR 过滤却看不出来。
2. **说不清的**：`median_distance` 等五列的 "distance" 没说清是**两两距离**；
   `z_score`/`p_value` 是裸名，暗示表里只有一个检验，但表里明明有第二个且用了显式命名——
   一张表两套后缀约定；`covered_genes` 的 "covered" 读起来像数据覆盖度，实际是"被打分成员"。

另有一处**不是命名问题**：`abundance_cv` / `conservation_cv` 的底层 feature 列由
`first_present_column()` 从候选里挑第一个命中的，表里不记录用了哪个，同一个列名在不同
features 表下代表不同的量。

## 决策

### 1. 漏斗补全（`n_measured_members`）

原表只有两层：`n_group_genes`（注释全集）和 `term_size`（DR<threshold 成员）。
中间少了一层"有拟合值的成员"，于是 `n_group_genes - term_size` 这个 gap 把两个不同机制
混成一个数：

实测（go_macrocomplex，全基因组 4513 个有拟合值 → 1559 个 DR < -0.3）：

| 组 | annotated | measured | scored | 覆盖率 | depleted |
|---|---|---|---|---|---|
| preribosome | 76 | 53 | 44 | 70% | 83% |
| cytosolic large ribosomal subunit | 85 | 46 | 29 | 54% | 63% |
| MCM core complex | 3 | 3 | 3 | 100% | 100% |

preribosome 的 gap（32）里 23 个是覆盖率、9 个是表型。**这直接关系到 D2 归因**：
`attribution.py` 的 docstring 自己承认 "technical issues (sparse insertion coverage,
curve-fit quality) are NOT auto-labelled"，而 `measured/annotated` 正是那个缺失的信号——
有了它才能把"覆盖不足导致的假不相干"和"真的异质性"分开。

实现上一次读取 `fitting_results.tsv` 服务两层：全集 = `n_measured_members` 层，
再用 `DR < threshold` 派生 permutation 零分布的点云。

### 2. 命名规则

两条规则覆盖整张表：

- **观察值列名 = `ZSCORE_METHODS` 的方法 key**（`median_pairwise_distance` 等），
  检验列 = `{method}_z` / `{method}_p`。于是整张表的几何部分可从方法列表推导，不再手写。
- **`n_*` 前缀专用于漏斗**（annotated ⊇ measured ⊇ scored），三者必须单调。

### 3. 映射表

| 旧 | 新 | 理由 |
|---|---|---|
| `centroid_x` / `centroid_y` | `geom_median_DR` / `geom_median_DL` | 是几何中位数不是质心；顺带标明轴 |
| `term_size` | `n_scored_members` | 它是 permutation 的 draw size，不是"组的大小" |
| `n_group_genes` | `n_annotated_members` | 与上者形成漏斗对照（长表契约同步改） |
| — | `n_measured_members` | **新增**，漏斗中间层 |
| `covered_genes` | `scored_member_names` | "covered" 误导为覆盖率 |
| `median_distance` | `median_pairwise_distance` | 明确是 pdist 归约 |
| `mean_distance` | `mean_pairwise_distance` | 同上 |
| `std_distance` | `std_pairwise_distance` | 同上 |
| `min_distance` | `min_pairwise_distance` | 同上 |
| `max_distance` | `max_pairwise_distance` | 同上（几何上 = 组直径） |
| `mpd` | **删除** | 与 `median_distance` 完全同值，同一张表两个名 |
| `z_score` | `median_pairwise_distance_z` | 统一 `{method}_z` 后缀 |
| `p_value` | `median_pairwise_distance_p` | 统一 `{method}_p` 后缀 |
| `mean_pairwise_distance_zscore` | `mean_pairwise_distance_z` | 同上 |
| `mean_pairwise_distance_p_value` | `mean_pairwise_distance_p` | 同上 |
| `p_fdr` | `q_value` | BH 校正后是 q 值不是 p 值；只有一个 FDR 列，不需要方法前缀 |
| `shared_fraction` | `frac_shared_members` | 说清是谁的比例 |
| — | `abundance_cv_feature` | **新增**，记录 CV 实际取自哪个 feature 列 |
| — | `conservation_cv_feature` | **新增**，同上 |

**不改的**（名字已准确）：`source` / `group_id` / `group_name` / `n_permutations` /
`abundance_cv` / `conservation_cv` / `paralog_fraction` / 去冗余的
`dag_depth` / `redundancy_cluster` / `cluster_size` / `is_representative` / `representative_*`。

### 4. `3a` 序号

`coherence.smk` → `3a_coherence.smk`，前缀同步带进 `results/3a_coherence/` 与
`logs/3a_coherence/`（`<chapter><section>` 约定，见 Snakefile 的 Includes 块）。
第 3 章 = "group-level analysis in fitness space"。clustering 及其下游仍未加前缀。

**前缀只到 rule 文件和 results/logs 目录**：`workflow/scripts/coherence/` 与
`workflow/src/coherence/` 保持无前缀（与 `pcr_qc`/`coverage`/`verification` 等脚本目录一致）。

## 迁移注意

**脚本不在 Snakemake 的 `input:` 里**（规则用 `shell:` 调脚本），所以改脚本 Snakemake
不会重跑。旧 schema 的 parquet 会被静默复用，下游读新列名直接 KeyError。
本次已删除 `results/3a_coherence/HD_DIT_HAP/` 下的全部派生产物，需重跑：

```bash
mamba activate snakemake
snakemake --use-conda --cores 8 results/3a_coherence/HD_DIT_HAP/coherence_terms_representatives.tsv
```

验证侧 notebook（`d1` / `d2` / `d2b` / `a2` / `x_closed_loop`）的列名已同步：
`{method}_zscore` → `{method}_z`、`covered_genes` → `scored_member_names`、
`centroid_*` → `geom_median_*`。这些 notebook 的产物同样是旧 schema，需重跑。

## 顺带修掉的隐患

`3a_coherence.smk` 里 `compute_coherence` 的 fallback 默认值仍是符号翻转前的
`_COH_CFG.get("dr_threshold", 0.3)`（config 与脚本默认都是 `-0.3`）。若该 key 有一天
从 config 缺失，DAG 会静默地筛选出**相反**的基因集（DR > +0.3 而非 DR < -0.3）。
已改为 `-0.3`。

## 验证

- `pytest tests/ -k coherence`
- 端到端：`prepare_annotation.py` → `compute_coherence.py`（真实 HD_DIT_HAP 数据）
  产出 173 组，与改名前一致；漏斗三层单调（173/173 行无违反）；
  preribosome 76/53/44、cytosolic LSU 85/46/29，与手工核算逐值相符。
- `--features` 路径下 `abundance_cv_feature` 正确记录
  `copies_per_cell_EMM_Proliferating_Cell`（第一个命中的候选）。
