# 聚类 finalize 多变体（variant）系统设计

**日期：** 2026-07-21
**范围：** 中等偏大 — 把今天早些时候合并进 `main` 的二元 `finalize_mode: auto|manual` 开关升级为**命名变体系统**：给定多种 finalize 策略，每种各自产出一份 `final_clusters.tsv`，enrichment 对所有变体平行跑（便于对比），ml/thesis 图表用其中被选中的那一份（`selected_variant`）。
**取代：** `docs/plans/2026-07-21-auto-finalize-clusters-design.md`（该文档描述的 auto/manual 二元开关是本设计的一个特例，见 §7 迁移）。

## 1. 背景与动机

今天早些时候（合并提交 `fc0f18f`）给核心链加了一条全自动 finalize 路径（kmeans 直接聚到 k=9 + 按 DR 编号），用 `config.clustering.finalize_mode: auto|manual` 二选一切换，并把最终簇列名统一成 `cluster`（手工版额外保留 `raw_cluster`）。

用户的真实需求比"二选一"更宽：希望**同时**跑多种 finalize 策略，每种都产出结果并平行跑下游富集，再根据指标或人工判断挑一个用于最终。需要覆盖的具体情况：

1. k-means 直接设定 cluster=9
2. k-means 设定 cluster=64 然后**自动**合并到 9
3. k-means 设定 cluster=64 然后**手工**合并到 9
4. 换其他方法（HC / GMM）重复 1–3
5. 手工在二维空间打格子分成 9 份
6. 以后可能新增的方法

二元 `finalize_mode` 装不下这个矩阵（方法 × 策略），所以升级为可扩展的命名变体系统。列契约（`cluster` / `raw_cluster`）保持不变。

## 2. 核心概念：variant

一个 **variant** = 一种具体的 finalize 策略实例，有唯一名字（如 `kmeans_direct9`）。每个 variant 声明一个 `type` 和该 type 需要的参数。四种 type 覆盖上面 6 种情况：

| type | 覆盖情况 | 机制 | 用到 64-候选 fan-out？ | 输出含 `raw_cluster`？ |
|---|---|---|---|---|
| `direct` | 1, 4a | `cluster_one_method(method, scaled, k=9)` → 按 DR 编号 | 否 | 否 |
| `auto_merge` | 2, 4b | 对 64 候选簇的**簇心**在 scaled (DR,DL) 空间做 ward 层级合并到 9 组 → 按 DR 编号 | 是（复用已有 `{method}_labels.pkl`） | 是 |
| `manual_merge` | 3, 4c | 人工在笔记本里把 64 原始 id 映射到合并组；最终编号由共享函数自动算，**人不再手工指定最终 id** | 是 | 是 |
| `grid` | 5 | `dr_cuts` / `dl_cuts` 轴切分组合成网格，落到原始 (DR,DL) 上；格子数必须 == `final_n_clusters` → 按 DR 编号 | 否 | 否 |

**统一编号规则（关键）：** 每种 type 最后都调用同一个共享函数 `renumber_by_dr(annotated, raw_labels, wt_cluster)`——按每组 `annotated["DR"]` 均值升序编号，最低 DR 组赋 `wt_cluster`（=9），其余按 DR 升序填 1..8，tiebreak 用 (mean DR, mean DL, 组原始 id)。**没有任何 type 手工指定最终 id**，包括 grid 和 manual_merge。这样格子/合并只负责"分组"，编号完全一致，跨变体可比。

## 3. 各 type 的确定性逻辑

放进 `workflow/src/clustering/candidates.py`（与现有确定性逻辑同处，便于单测）。

### 3.0 共享编号（从今天的 `auto_finalize` 抽出）
```
renumber_by_dr(annotated, raw_labels: pd.Series, n_clusters, wt_cluster) -> pd.Series(final cluster, index=raw_labels.index)
```
逐组算 mean DR / mean DL，稳定排序，最低 DR→wt_cluster，其余升序。组数 != n_clusters 时抛 ValueError（防止静默错编号）。

### 3.1 `direct`（取代今天的 `auto_finalize`，泛化出 method 参数）
```
finalize_direct(annotated, scaled, method, n_clusters=9, random_state, wt_cluster) -> annotated + `cluster`
```
= 今天 `auto_finalize` 的行为，只是把硬编码的 `BEST_METHOD` 换成显式 `method` 入参（default 仍 kmeans，保证 `kmeans_direct9` 与今天 `finalize_mode: auto` 逐字节一致）。

### 3.2 `auto_merge`
```
finalize_auto_merge(annotated, scaled, raw_labels_64: pd.Series, n_clusters=9, wt_cluster) -> annotated + `cluster` + `raw_cluster`
```
1. 输入是某方法在 k=64 时的原始标签（复用 `cluster_one_method` 已产出的 `_work/{method}_labels.pkl`，**不新增 fan-out**）。
2. 算每个 64-簇在 scaled (DR,DL) 空间的簇心（各簇均值），得 64×2 矩阵。
3. 对这 64 个簇心做 `AgglomerativeClustering(n_clusters=9, linkage="ward")`，得到 64→9 的合并组标签（簇心不按簇大小加权——最简单可辩护的默认，记为假设）。
4. 把合并组标签映射回每个基因（gene→raw64→mergegroup），交给 `renumber_by_dr` 编号。
5. 输出保留 `raw_cluster` = 原始 64 标签。

### 3.3 `grid`
```
finalize_grid(annotated, scaled, dr_cuts, dl_cuts, n_clusters=9, wt_cluster) -> annotated + `cluster`
```
1. 用 `dr_cuts` 把 DR 轴切成 len(dr_cuts)+1 段，`dl_cuts` 同理切 DL 轴；组合成 (len(dr_cuts)+1)×(len(dl_cuts)+1) 个矩形格子。
2. 校验格子数 == `final_n_clusters`（否则抛 ValueError，提示调 cuts）。
3. 每个 scaled 基因按 (DR,DL) 落格（`np.digitize`），得格子 id → 交给 `renumber_by_dr` 编号。
4. 落格用**原始 scaled 值**（即已 cap/除过的 DR/DL，与聚类同一空间），cuts 语义即 scaled 空间阈值——在设计文档和 config 注释里写清。

### 3.4 `manual_merge`（笔记本，人工判断）
不写 Snakemake 规则（人工步骤）。笔记本泛化：
- 参数 `METHOD`（选哪个方法的 64 候选）+ `VARIANT`（决定输出路径）。
- 直接读 `_work/{method}_labels.pkl` 拿 64 原始标签 + `annotated_data.pkl` 拿注释表（不再依赖 candidate_clusters.tsv 的 `cluster` 列，那列只在 kmeans 下才等于 64 候选）。
- 人工只维护**一个** dict：`raw64 -> merge_group`（合并组任意整数标签）。
- 最终编号调 `renumber_by_dr`（去掉今天那第二个手工 reorder dict——人不再手工排号）。
- 写到 §4 的 per-variant curated 路径，输出含 `cluster` + `raw_cluster`。

## 4. 列契约（不变）

`final_clusters.tsv` 最终簇列一律 `cluster`（1..9，WT=9）；有合并步骤的（auto_merge / manual_merge）额外保留 `raw_cluster`（合并前原始标签）。完整注释表 + `A, DR, DL` + `cluster`(+`raw_cluster`)，index 名 `Systematic ID`，`index=True` 写出。下游据此逐列一致读取。

## 5. 路径约定

- 可构建变体（`direct` / `auto_merge` / `grid`）：`results/clustering/{dataset}/{variant}/final_clusters.tsv`（Snakemake 产出，可删可重跑）。
- curated 变体（`manual_merge`）：`resources/curated/final_clusters/{dataset}/{variant}.tsv`（人工维护，版本控制，不可重跑）。

## 6. config 形态（`config/analysis.yaml` 的 `clustering:` 下）

**破坏性变更：** 删掉今天加的 `finalize_mode` / `finalize_mode_overrides`（今天才合并，无外部消费者）。新增：
```yaml
  final_n_clusters: 9            # 所有变体的最终簇数 k（沿用）
  variants:
    kmeans_direct9:       {type: direct,       method: kmeans}
    kmeans_merge9_auto:   {type: auto_merge,   method: kmeans}
    kmeans_merge9_manual: {type: manual_merge, method: kmeans}
    # 追加方法示例：
    # hc_direct9:  {type: direct,     method: hierarchical_agg}
    # grid9:       {type: grid,       dr_cuts: [0.3, 0.6, 0.9], dl_cuts: [2.0, 5.0]}
  selected_variant: kmeans_direct9        # ml/thesis 用哪一个；默认等价今天的 finalize_mode: auto
  selected_variant_overrides: {}          # 按 dataset 覆盖，如 {HD_DIT_HAP: kmeans_merge9_manual}
```
`wt_cluster` 已有（`enrichment` 下），编号沿用它 = 9。

## 7. 从今天 auto/manual 二元开关的迁移

- `finalize_mode: auto`  ≡ 只配一个 `direct` 变体并设为 selected。
- `finalize_mode: manual` ≡ 只配一个 `manual_merge` 变体并设为 selected。
- `final_clusters_path(dataset)` → `final_clusters_path(dataset, variant)`（多一个 variant 维度）。
- `selected_variant(dataset)` 新增：给 ml.smk 选唯一变体。

## 8. Snakemake 接线

`clustering.smk`：
- 保留 prepare spine + cluster_one_method(×4, k=64) + select（一行不动，`auto_merge`/`manual_merge` 复用其 `{method}_labels.pkl`）。
- 新 helper：`final_clusters_path(dataset, variant)`、`selected_variant(dataset)`、按 type 分组的变体名列表（供 wildcard 约束）。
- 新规则 `finalize_direct` / `finalize_auto_merge` / `finalize_grid`，各带 `variant` 通配符，约束到 config 里对应 type 的变体子集（仿现有 `method` 通配符约束写法）。`manual_merge` 无规则，沿用"缺输入即报错、去跑笔记本"的语义，只是按 variant 分文件。

`enrichment.smk`：`prepare_genesets` 及下游 ontology/finalize 规则新增 `variant` 通配符维度，输出目录多一层 `.../{variant}/...`，`final_clusters` 输入用 `final_clusters_path(wc.dataset, wc.variant)`。这样每个变体各自跑一套富集，可对比。

`ml.smk`：`prepare_ml_data` 用 `final_clusters_path(dataset, selected_variant(dataset))`——单变体，不 fan-out。

## 9. 改动清单

**新增/改 src：** `candidates.py`：抽 `renumber_by_dr`；`auto_finalize`→`finalize_direct`(+method)；加 `finalize_auto_merge`、`finalize_grid`。
**新增/改脚本：** `auto_finalize_clusters.py`→`finalize_direct_clusters.py`；加 `finalize_auto_merge_clusters.py`、`finalize_grid_clusters.py`。
**改规则：** clustering.smk（helper + 3 规则）、enrichment.smk（variant 维度）、ml.smk（selected_variant）。
**改 config：** 见 §6。
**改笔记本：** finalize_gene_clusters.ipynb（§3.4）。
**改测试：** rename/扩 test_clustering.py（renumber_by_dr / finalize_direct / auto_merge / grid + 各 driver）；enrichment/ml 测试若涉及路径维度同步。
**改文档：** 本文档 + CLAUDE.md + README。

## 10. 验证

- 单测：`renumber_by_dr` 编号确定性 + 最低 DR=9 + tiebreak；`finalize_direct` 与今天 `auto_finalize` 行为一致；`finalize_auto_merge` 能从合成数据恢复已知 9 组；`finalize_grid` 落格正确 + 格子数不符时 ValueError；各 driver 写出正确 schema。
- 全量 `pytest` 通过。
- dry-run 确认 enrichment 按 variant fan-out、ml 指向 selected_variant。
- 对 `HD_DIT_HAP` 实跑每种 buildable type 至少一个变体，确认产出 9 簇、`cluster` 1..9、WT=9 最低 DR、`auto_merge`/manual 带 `raw_cluster`。
