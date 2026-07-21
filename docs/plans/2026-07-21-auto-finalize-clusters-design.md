# 自动 finalize 聚类阶段设计

**日期：** 2026-07-21
**范围：** 中等 — 给核心链的手工步骤 `notebooks/clustering/finalize_gene_clusters.ipynb` 增加一条纯 Snakemake、无需人工判断的全自动替代路径，使 `final_clusters.tsv` 不再只能由笔记本产出。手工笔记本原样保留为可选微调入口。同时统一最终簇列名契约为 `cluster`。

## 1. 背景与目标

原设计（`2026-07-15-DIT-HAP-analysis-design.md` §5）刻意把 64→9 的簇合并留作人工笔记本：它靠肉眼在 (DR, DL) 特征空间上把 64 个候选簇手挑映射到 9 个最终簇，没有任何可复现的阈值规则。这带来两个问题：

- 核心链因此不是单一 DAG——每次要产出 `final_clusters.tsv` 都得手动开笔记本。
- 笔记本里的两个 merge 字典（`reformat_cluster`、`reorder_reformat_cluster`）按当前 kmeans 的具体候选整数 id 硬编码，换 dataset 或换算法即失效。

**目标：** 提供一条确定性、无人工判断、纯 Snakemake 的路径，从 `prepare_clustering_data` 的中间产物直接产出符合下游契约的 `final_clusters.tsv`。这不是"复现"人工手挑映射（不可能），而是另起一条确定性路线（直接聚类到 k=9）。手工路径继续保留，仅在需要人工微调时使用。

## 2. 架构与 DAG

复用现有 `prepare_clustering_data` spine（它产出的 `annotated_data.pkl` + `scaled_data.pkl` 与 64-候选那套共享），但绕开 64-候选的四方法 fan-out。自动路径只跑 kmeans 到 k=9，再确定性编号。

```
现状（手工）：
  prepare_clustering_data → cluster_one_method×4(k=64) → select → candidate_clusters.tsv
                                                              ↓ [人工笔记本]
                                                    resources/curated/final_clusters.tsv → enrichment/ml

新增（自动，默认）：
  prepare_clustering_data → auto_finalize_clusters(kmeans k=9 + DR 编号)
                                    → results/clustering/final/{dataset}/final_clusters.tsv → enrichment/ml
```

- 64-候选 pipeline + 笔记本原样保留，一行不删——只是不再是唯一入口。
- 自动路径只吃 `prepare_clustering_data` 的两个 pkl，所以自动模式下不需要跑 4×64 的 fan-out，反而更轻。
- 下游 `enrichment.smk` / `ml.smk` 里硬编码的 `final_clusters="resources/curated/final_clusters.tsv"` 改成 input 函数，按 `config.clustering.finalize_mode` 选择路径，支持按 dataset 覆盖。

## 3. 自动 finalize 逻辑

新共享函数 `auto_finalize()` 放进 `workflow/src/clustering/candidates.py`（与现有确定性逻辑同处，便于测试）；新脚本 `workflow/scripts/clustering/auto_finalize_clusters.py` 作为薄 CLI 包装。

`auto_finalize(annotated, scaled, n_clusters=9, random_state, wt_cluster=9)`：

1. `labels = cluster_one_method("kmeans", scaled, n_clusters, random_state)` — 复用现有函数，0-based 原始标签，与 scaled 行对齐。
2. **确定性编号**（替代人工肉眼 merge）：
   - 按原始标签算每簇 `annotated["DR"]` 均值。
   - 均值最低的簇判定为 WT，赋 `wt_cluster` id（9）。
   - 其余 8 簇按 DR 均值升序编号 1..8（1 = 次低 DR）。
   - tiebreak：DR 均值相同时用 (mean DR, mean DL) 二级排序 + 原始标签号做最终 tiebreak，保证跨运行完全确定。
3. 编号写进最终 `cluster` 列（见 §4 契约）。

**WT 语义依据：** 低 DR = viable = WT 样（如 gls2 DR=0.059），高 DR = essential（如 taf11 DR=1.06）。最低 DR = WT = 9 与现有 `wt_cluster=9`、"< wt_cluster 做 nonWT 比较"约定完全兼容，不用改 config 的 wt 值。

## 4. 列契约统一

最终契约：`final_clusters.tsv`（auto 或 manual）最终簇列一律叫 `cluster`。

| | `cluster`（最终簇 1..9） | `raw_cluster`（合并前原始标签） |
|---|---|---|
| **auto** | kmeans k=9 → 按 DR 均值确定性编号 | 无（没有合并步骤） |
| **manual** | 笔记本 merge+reorder 后的 9 簇 | k=64 候选原始标签 |

两种产出都是完整注释表（index 名 `Systematic ID`，`index=True` 写出）+ `A, DR, DL` + 最终 `cluster` 列。schema 逐列一致，下游据此读取。

## 5. config 新增

`config/analysis.yaml` 的 `clustering:` 下：

```yaml
finalize_mode: auto            # auto | manual
finalize_mode_overrides: {}    # 如 {HD_DIT_HAP: manual}
final_n_clusters: 9
# wt_cluster 已有（enrichment/ml 下），auto 编号沿用它 = 9
```

## 6. 改动清单

**新增：**
- `workflow/src/clustering/candidates.py`：`auto_finalize()` 函数。
- `workflow/scripts/clustering/auto_finalize_clusters.py`：CLI 脚本。
- `workflow/rules/clustering.smk`：`auto_finalize_clusters` 规则（输出 `results/clustering/final/{dataset}/final_clusters.tsv`）。
- `tests/test_clustering.py`（或新测试）：`test_auto_finalize`——编号确定性、最低 DR=9、tiebreak。

**下游读取改动（`revised_cluster` → `cluster`，最终契约重命名）：**
- `workflow/src/enrichment/cluster_enrichment.py:71`：`CLUSTER_COLUMN = "revised_cluster"` → `"cluster"`。
- `workflow/src/ml/data.py:57`：`revised_cluster` → `cluster`（仍映射成 `DIT_HAP_cluster`）。
- `workflow/scripts/ml/prepare_features_targets.py:189`：同上。
- 三处及相关脚本 docstring 里 `revised_cluster` 字样一并更新。

**选择机制改动：**
- `workflow/rules/enrichment.smk` / `workflow/rules/ml.smk`：把硬编码的 `final_clusters=` 改成按 `finalize_mode`（+ overrides）返回路径的 input 函数。

**manual 笔记本改动：**
- `notebooks/clustering/finalize_gene_clusters.ipynb`：读 `candidate_clusters.tsv`（其列名为 `cluster`=64 原始）→ 先 `rename(cluster → raw_cluster)` → merge/reorder 结果写进新 `cluster` 列。产出为完整注释表 + `cluster` + `raw_cluster`。

**测试改动（造数据列名同步）：**
- `tests/test_cluster_enrichment.py`、`tests/test_prepare_features_targets.py`、`tests/test_train_automl.py`：造数据的 `revised_cluster` → `cluster`。

## 7. 验证

- 新增单测覆盖 `auto_finalize` 的编号确定性（同一 seed 两次结果一致）、最低 DR 簇被赋 9、tiebreak 稳定。
- 全量 `pytest` 通过（含改名后的三个下游测试）。
- 端到端 dry-run（`snakemake -n`）确认 auto 模式下 enrichment/ml 的 DAG 指向 `results/clustering/final/...`、manual 模式指向 `resources/curated/...`。
- 条件允许时对 `default_dataset` 实跑 auto 路径，确认产出 `final_clusters.tsv` 的 9 个簇、`cluster` 列 1..9、WT=9。
