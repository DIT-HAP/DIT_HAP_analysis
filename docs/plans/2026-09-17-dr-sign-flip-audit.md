# DR / LFC 符号翻转 —— 下游影响审计

日期：2026-09-17
状态：**审计已完成，改动已应用；按用户要求未跑任何程序测试** —— 用户将逐 rule 审核分析后再做详细测试。
适用仓库：`DIT_HAP_analysis`（消费 `DIT_HAP_snakemake` 的 `release/`）

---

## 1. 上游变更（实测确认）

对比 `projects/HD_DIT_HAP/arc/release/`（2026-07，旧）与 `projects/HD_DIT_HAP/release/`（2026-09-17，新）：

| 量 | 旧 | 新 | 关系 |
|---|---|---|---|
| `DR`（gene_level） | -0.168 … **+2.091**，均值 +0.303 | **-2.091** … +0.168，均值 -0.301 | 取反 |
| `LFC`（各时间点均值） | +0.515 / +1.213 / +1.666 / +2.042 | -0.515 / -1.213 / -1.666 / -2.042 | 取反（精确） |
| `DL` | 0 … 9.511 | 0 … 9.511 | **未翻转** |

- DR 有 87.6%（3952/4513）是逐位精确取反，与取反的相关性 0.9993。剩下 561 个都是接近 0 的值，差异来自拟合本身重跑，不是符号问题。
- DL 有 94.8% 逐位相同、`|Δ|` 中位数为 0；但均值 0.922 → 0.885，**拟合确实重跑过**，所以即便不论符号，所有结果都需要重新生成。

**新约定：DR / LFC 越负，depletion 越严重。DR ≈ 0 表示不耗竭（WT 侧）。**

新 DR 分位数（n=4513）：0% -2.091 / 5% -1.071 / 25% -0.637 / 50% -0.065 / 75% 0.001 / 95% 0.029 / 100% 0.168。

---

## 2. 会静默产出错误科学结论的（最高优先级）

这三处不报错、不警告，只是结果反了。

### 2.1 `comparison` —— 与 gRNA 研究的相关性符号翻转

- **位置**：`workflow/src/comparison/core.py:153,185-193`；曲线在 `prepare_fitness_table.py`
- **机制**：该阶段把 DIT-HAP 的 `DR` 与 curated 的 `260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv` 的 `um` 做相关/回归。那张 gRNA 表是**冻结的 curated 文件，不会随上游重跑**，仍是旧约定（`um` ∈ [-0.167, +1.805]，均值 +0.331）。
- **实测**：配对 4465 个基因，`corr(DIT-HAP DR, gRNA um)` 由 **+0.920 → -0.920**。
- **后果**：这一阶段存在的全部意义就是"两个研究是否一致"，改完会报出强**反**相关，而实际一致度是 r=+0.92。
- **建议**：在读入 gRNA 表时把 `um`（必要时含 `lam` 相关列）取反，使其对齐到新约定；或者在比较前把 DIT-HAP 侧取反。**二选一并写清注释**，不要两边都翻。
- **验证**：重跑后 `fitness_correlation_stats.tsv` 的相关系数应回到 ≈ +0.92。

### 2.2 `clustering` —— WT 簇与最耗竭簇标签对调

- **位置**：`workflow/src/clustering/candidates.py:191-195`（`renumber_by_dr`）
- **现状**：docstring 写"the group with the **lowest** mean DR becomes wt_cluster (WT)"。
- **新数据下**：最低 mean DR = 最严重 depletion，于是 WT 拿到最耗竭那一组，而真正的 WT 组拿到 9。
- **后果**：`final_clusters.tsv` 的 `cluster` 列 1..9 全部错位。下游 `enrichment`（按 cluster 号分基因集）、`ml`（`wt_cluster: 9` 定义 nonWT split）、`comparison`（读 curated cluster 表）全部跟着错，且不会报错。
- **建议**：改为"**最高** mean DR = WT"。这是整套下游标签的地基，**先改这个再改别的**，否则重跑出来的中间产物都要再废一次。
- **验证**：改后跑一次，检查 WT 组的 mean DR 是否最接近 0（最大），以及 cluster 9 的 mean DR 是否最负。

### 2.3 `annotate` —— 工作簿里两列 DR 符号相反

- **位置**：`workflow/scripts/annotate/build_annotation_reference.py:12-13,28-29,285`
- **机制**：同一个 workbook 里并排写"Gene-level HD_DIT_HAP depletion (DR/DL)"和"gRNA-level depletion (DR/DL)"。
- **后果**：两列符号相反，读表的人会把一致解读成矛盾。与 2.1 同源，处理方式应与 2.1 保持一致。
- **验证**：生成 `{dataset}_annotated.xlsx` 后，两列在同一基因上的符号应同向。

---

## 3. 阈值类：新数据下选空或全选

以下均为实测（n=4513）。因为 DR 是精确取反，阈值取反后选中的基因集合**基本等价**（差异只来自拟合重跑的微调），所以映射是机械的。

| # | 位置 | 现在 | 新数据命中 | 建议改为 | 取反后命中 |
|---|---|---|---|---|---|
| 1 | `workflow/src/ml/data.py:46,64` | `DR > 0.3` | **0** | `DR < -0.3` | 1559 |
| 2 | `config/analysis.yaml:86` `ml.dr_filter` | `0.3` | — | `-0.3` | — |
| 3 | `config/analysis.yaml:84` `ml.dr_split_threshold` | `0.35` | `>0.35`→**0**，`<=0.35`→4513 | `-0.35` | `<`→1498，`>=`→3015 |
| 4 | `workflow/scripts/ml/prepare_features_targets.py:248-249` | `DR > / <= 0.35` | 同 3 | `DR < / >= -0.35` | 同 3 |
| 5 | `workflow/src/verification/core.py:146` `WT2nonWT` | `DR > 0.35` | **0** | `DR < -0.35` | 1498 |
| 6 | `verification/core.py:147-148` `scE2E`/`sc2E` | `DR > 0.75` | **0** | `DR < -0.75` | 938 |
| 7 | `verification/core.py:149` `E2V` | `DR < 0.35` | **4513**（全选） | `DR > -0.35` | 3015 |
| 8 | `config/analysis.yaml:112` `coherence.dr_threshold` | `0.3` | — | `-0.3` | — |
| 9 | `workflow/scripts/coherence/compute_coherence.py:168,224,349` | `fitting["DR"] > dr_threshold` | **0**（背景集空） | `DR < dr_threshold` | 1559 |
| 10 | `workflow/src/domain_differences/core.py:77,168` | `DR_THRESHOLD = 0.15`，`DR > 0.15` | **6** | `DR < -0.15` | 1805 |

**注意 7 和 8/9 的两点：**

- #7 `E2V` 的 `"sort": "asc"` 方向要跟着一起定：旧语义是"DR 小 = 轻表型"。取反后要重新确认这一组该按哪个方向排。
- #8/#9 的 `dr_threshold` 不只是筛基因，它还定义 coherence 的"**背景集**"（`background = fitting[DR > thr]`）。改成 `DR < -0.3` 后背景集变成 1559 个深度耗竭基因——**这个语义变了**：原来背景是"有意义的耗竭基因"，现在取反后仍然是"有意义的耗竭基因"，所以是对的。但请人工确认一次背景集大小符合预期（原来约 1559 上下）。

另外几处只是文案，但会误导，一并改：`workflow/scripts/coherence/compute_coherence.py:9,343,345`、`plot_coherence.py:213`、`deduplicate_terms.py:20`、`workflow/scripts/ml/prepare_ml_data.py:9,97`、`workflow/src/domain_differences/core.py:8,34,162,194`、`workflow/src/ml/data.py:7,54,66`。

---

## 4. 数值变换类

### 4.1 `DR_CAP`：从上限变下限

- **位置**：`workflow/src/clustering/candidates.py:66`（`DR_CAP = 1.3`）、`:118`（`x if x < dr_cap else dr_cap`）；`config/analysis.yaml:10` `clustering.dr_cap`；`workflow/scripts/clustering/prepare_clustering_data.py:133`（`--dr-cap` help）
- **现状**：裁剪 DR 的上限（抑制极端值）。
- **新数据下**：DR 最大 0.168，**永不触发**（旧数据有 7 个基因被裁到 1.3）。这个 byte-faithful 的 quirk 静默失效。
- **建议**：改为下限 `x if x > -dr_cap else -dr_cap`，配置写 `-1.3`（或把常量改成负值并改名 `DR_FLOOR`，见 §8 命名建议）。
- **验证**：新数据下 `DR < -1.3` 应命中 **7** 个基因，与旧数据 `DR > 1.3` 的 7 个精确对应。

### 4.2 `grid` 变体的 DR 切点：网格塌成一列

- **位置**：`config/analysis.yaml:31` `grid9: {dr_cuts: [0.4, 0.8], ...}`；消费方 `workflow/src/clustering/candidates.py:298`（`np.digitize(scaled["DR"], sorted(dr_cuts))`）
- **现状**：切点是**缩放空间**的阈值。
- **新数据下实测**：缩放后 DR < 0.4 的占 **100%**，`[0.4,0.8)` 与 `>=0.8` 各 **0** 个 → 整个 3×3 网格塌成一列，所有基因落进同一个 cell。
- **建议**：`dr_cuts: [-0.8, -0.4]`。实测分箱为 **847 / 583 / 3083**，与旧数据按 `[0.4, 0.8]` 的 **3079 / 585 / 849** 精确镜像。
- **验证**：重跑 `grid9` variant，三个 DR bin 的基因数应约为 850 / 585 / 3080。

### 4.3 coherence 的归一化：**更正——不需要改数值**

> **2026-09-17 更正**：初版把这里写成"整个 DR-DL 空间越界"，夸大了。复查后：coherence
> 的**所有**指标（median/mean/max pairwise distance、knn distance、geometric median）
> 都是这个空间里的欧氏距离，而坐标轴取反是**等距变换**——所有点对距离不变，因此
> **z-score 与 p 值完全不变**。而且两处"归一化"其实是恒等变换：
> `compute_coherence.py` 的 `_min_max_normalize(DR, 0.0, 1.0)` = `(DR-0)/(1-0)` = DR，
> `attribute_incoherence.py` 的 `_DR_NORM = 1.0` 是除以 1。`metrics.py` 的
> `normalize_dr_dl` / `DR_NORM_MAX` **无任何调用方**（死代码）。

- **位置**：`workflow/src/coherence/metrics.py:17`、`workflow/scripts/coherence/compute_coherence.py:88`、`workflow/scripts/coherence/attribute_incoherence.py:103`
- **实际影响**：只有**图形朝向**（散点图的 DR 轴翻到负半轴）和文档措辞，指标数值不变。
- **已做**：数值一律不动（改数值反而会改变 DR/DL 两轴的相对权重，凭空改变结果）；只更新了文档，说明新 DR 范围与"距离类指标对取反映射不变"。
- **连带**：`workflow/src/coherence/attribution.py` 的示例 DR 值与 "masked (lower DR)" 措辞按新方向改写；GMM 的 core/minor 判定基于离散度，方向无关，逻辑未动。

---

### 4.4 coverage 的 DR 直方图 bin 范围：**初版漏掉的硬伤**

> **2026-09-17 补充**：初版说 coverage 整体不受影响，只对了一半。判"被覆盖"确实与符号
> 无关，但 coverage 还要**画 DR 直方图**，而 bin 边界是写死的。

- **位置**：`workflow/src/coverage/core.py:85`（`DR_BINS`）、`:397` 与 `:427-437`（详细基因表按 DR 排序）
- **现状**：`DR_BINS = np.arange(-0.2, 1.5, 0.05)`，覆盖区间 [-0.2, 1.45]。
- **新数据下实测**：只有 **61.9%**（2795/4513）的基因落在 bin 范围内，**1718 个（38.1%）被静默丢弃**——而且丢的正是最耗竭的那一批，也就是信号本身。旧数据同口径是 99.9%。
- **已做**：`DR_BINS = np.arange(-1.45, 0.25, 0.05)`，即旧边界的**精确取反**，覆盖 99.8%，仍是截断极端耗竭尾（旧的在顶部截 0.1%，新的在底部截约 1%）。详细基因表的 DR 排序由降序改为升序（最耗竭在前）。`DL_BINS` 未动。
- **验证**：重跑后 DR 直方图应有数据铺满横轴，而不是挤在右端。

---



## 5. 出图 / 公式类

- **位置**：`workflow/src/plotting/gene_level.py`
- **这是活代码**：`verification/core.py`、`plot_variant_clusters.py`、`plot_all_variant_clusters.py`、`plot_group_scatter.py` 都在用。
- 涉及点：
  - `:167` `sigmoid_gompertz(A, DR, DL)` —— DR 是**公式参数**不是阈值
  - `:175` `xend = max(DL + A / DR, 1)` —— `A/DR` 随符号变号，曲线终点会跳
  - `:206,237,309,321,344,372` —— DR 作为 x 轴 feature（`visualize_cluster_on_feature_space`、`plot_cluster_on_axis`、`plot_given_genes_on_feature_space`）
- **建议**：已改为在 `sigmoid_gompertz` 与切线段里用 `abs(DR)`——DR 在那里是"最大耗竭速率"，是量值不是带方向的阈值。这样旧数据下 DR ≥ 0 的基因行为逐位不变，同时也修掉了旧约定下 DR 略为负的近 WT 基因（原本会画出倒置曲线）。
- **注意**：这是全清单里唯一一处**不是机械取反**的改动，**请务必出图人眼核对**曲线形状。

---

## 6. 不受影响（列出来防止误改）

| 位置 | 为什么不受影响 |
|---|---|
| `workflow/src/coverage/core.py:289` 判"被覆盖" | 用的是 `DR.notna()`，与符号无关。**但 coverage 并非整体不受影响——见 §4.4。** |
| 所有 `DL` 相关计算 | DL 未翻转（但拟合重跑过，数值有微调）。 |
| `workflow/scripts/annotate/build_annotation_reference.py` 的 `gRNA_DR` | 来自 curated gRNA 表，是**另一个数据源**，不能跟着翻。见 §2.3。 |
| `workflow/src/utr/core.py:14-15,306` 的 `um_ratio` | 是 `insertion DR / gene DR` 的**比值**，分子分母同时翻号，比值不变。 |
| `workflow/src/noncoding_rna/core.py:203,211` 的 DR 降序排序 | 只是显示顺序，不是阈值。可选调整，不属错误。 |
| `data_config.lfc()` / `LFCs` | 声明了 `LFC.tsv` 但全仓库**零调用**，LFC 的翻转暂时咬不到任何代码。 |
| `enrichment` / `comparison` / `utr` / `noncoding_rna` 里的 `um`/`lam` → `DR`/`DL` 改名逻辑 | 纯改名，与符号无关。 |

---

## 7. 顺带发现（与本次无关，但该处理）

- **`workflow/scripts/ml/prepare_features_targets.py` 已无任何规则引用**（`ml.smk` 调的是 `prepare_ml_data.py`），但里面仍有 `DR_gt_p35` / `DR_le_p35` 的切分逻辑，且 `tests/test_prepare_features_targets.py` 还在跑。要么随本次一起改，要么直接删。
- **`tests/test_clustering.py:49`** 断言 `scaled["DR"] == [0.5, 1.3, DR_CAP, 1.29]`，直接钉死了 cap 的旧方向，改 §4.1 时必须同步改。
- **`tests/test_train_automl.py:82`**、**`tests/test_verification.py:137-162`** 都用旧方向的 DR 造数据，会随 §3 一起失败——这是好事，它们就是回归网。

---

## 8. 建议的执行顺序

按"改动影响面 × 是否阻塞下游"排：

1. **clustering（§2.2 + §4.1 + §4.2）** —— cluster 标签是所有下游的输入，先改它避免下游返工。
2. **coherence（§3 #8/#9 + §4.3）** —— 阈值和归一化必须一起改，否则归一化空间仍然是错的。
3. **ml（§3 #1-#4）** —— 阈值改完建模集才能非空。
4. **verification（§3 #5-#7）** —— 阈值 + 排序方向。
5. **domain_differences（§3 #10）**。
6. **comparison + annotate（§2.1 + §2.3）** —— 需要先定"翻哪一边"，是唯一一个需要**产品决策**而不是机械改写的。
7. **plotting/gene_level.py（§5）** —— 单独一轮，出图对比后再动。

### 命名建议

`DR_CAP` / `dr_cap` / `--dr-cap` 在翻转后名不副实（它变成了下限）。建议同时改名为 `DR_CLAMP` / `dr_clamp`，避免下一个人再踩一次。

---

## 9. 每处的统一验证方法

1. `snakemake -n <target>` 确认 DAG 能解析（不报 missing input）。
2. 跑完后检查三个关键数字：
   - `final_clusters.tsv` 里 cluster 9（WT）的 mean DR 应为**最大**（最接近 0）
   - `results/ml/models/.../modeling_data.parquet` 行数应与旧版同量级（约 1559），不是 0
   - `results/comparison/.../fitness_correlation_stats.tsv` 的相关系数应为 **≈ +0.92**，不是 -0.92
3. `pytest` —— `test_clustering.py` / `test_train_automl.py` / `test_verification.py` 是最直接的回归网。

---

## 10. 已应用的改动清单（2026-09-17）

未跑任何程序测试；只做了 Python 语法解析与 YAML 解析校验（均通过）。

### 10.1 clustering（枢纽，标签地基）

| 文件 | 改动 |
|---|---|
| `workflow/src/clustering/candidates.py` | `renumber_by_dr` 排序由降序改升序 → **最高 mean DR = WT**；`DR_CAP = 1.3` 改名 `DR_CLAMP = -1.3` 并由上限改下限（`x if x > clamp else clamp`）；相关 docstring/注释 |
| `config/analysis.yaml` | `clustering.dr_cap: 1.3` → `dr_clamp: -1.3`；`grid9.dr_cuts: [0.4, 0.8]` → `[-0.8, -0.4]` |
| `workflow/rules/clustering.smk` | 参数名与默认值 `dr_clamp=-1.3`，CLI `--dr-clamp` |
| `workflow/scripts/clustering/prepare_clustering_data.py` | 同上重命名，help 文案 "ceiling" → "floor" |
| `tests/test_clustering.py` | 夹具 DR 取反，断言改下限语义，两个测试改名 |

### 10.2 ml

| 文件 | 改动 |
|---|---|
| `workflow/src/ml/data.py` | `DR_FILTER = -0.3`；`query("DR > …")` → `query("DR < …")`；日志与文档示例 |
| `config/analysis.yaml` | `ml.dr_split_threshold: -0.35`、`ml.dr_filter: -0.3` |
| `workflow/scripts/ml/prepare_ml_data.py` | help 与 docstring |
| `workflow/scripts/ml/prepare_features_targets.py` | `DR_SPLIT_THRESHOLD = -0.35`；split 键 `DR_gt_p35/DR_le_p35` → `DR_lt_p35/DR_ge_p35`，比较符同步翻转。**注意：此脚本仍无任何规则引用** |

### 10.3 verification

| 文件 | 改动 |
|---|---|
| `workflow/src/verification/core.py` | `_CRITICAL_GROUPS` 四个 filter 阈值取反，且 `sort` 方向全部镜像（desc↔asc），保持同一批基因同一顺序 |
| `tests/test_verification.py` | 夹具 DR 取反，filter/sort/断言同步 |

### 10.4 coherence

| 文件 | 改动 |
|---|---|
| `config/analysis.yaml` | `coherence.dr_threshold: -0.3` |
| `workflow/scripts/coherence/compute_coherence.py` | `DR > thr` → `DR < thr`；dataclass 默认、CLI 默认、usage 示例、全部 help 文案 |
| `workflow/scripts/coherence/{plot_coherence,deduplicate_terms,attribute_incoherence}.py` | 文案 |
| `workflow/src/coherence/metrics.py` | **数值未动**，仅重写文档说明"距离类指标对坐标轴取反映射不变" |
| `workflow/src/coherence/attribution.py` | 示例 DR 值与 "masked (lower DR)" 措辞 |

### 10.5 domain_differences

| 文件 | 改动 |
|---|---|
| `workflow/src/domain_differences/core.py` | `DR_THRESHOLD = -0.15`；函数改名 `filter_high_dr_genes` → `filter_depleted_genes`，比较 `>` → `<`；调用参数 `high_dr_genes` → `depleted_genes` |
| `workflow/scripts/domain_differences/compute_domain_stats.py` | **校验守卫由"必须非负"改为"必须为负"（原本会 100% 抛错）**；变量/文案 |
| `tests/test_domain_differences.py` | 阈值断言、夹具、两个测试改名 |

### 10.6 comparison + annotate（翻转对象是 gRNA 侧）

| 文件 | 改动 |
|---|---|
| `workflow/src/comparison/core.py` | 新增 `GRNA_METRIC_SIGN = -1.0`，合并后对 `um_gRNA` 取反 |
| `workflow/src/annotation/core.py` | 新增 `GRNA_DR_SIGN = -1.0`，`build_grna_block` 对 `gRNA_DR` 取反 |
| `tests/test_annotation.py` | `gRNA_DR` 断言 1.048 → -1.048 |

### 10.7 coverage

| 文件 | 改动 |
|---|---|
| `workflow/src/coverage/core.py` | `DR_BINS = arange(-1.45, 0.25, 0.05)`（旧边界精确取反）；详细基因表 DR 排序由降序改升序；docstring |

### 10.8 出图公式

| 文件 | 改动 |
|---|---|
| `workflow/src/plotting/gene_level.py` | `sigmoid_gompertz` 的 `alpha` 改用 `abs(DR)`；切线段 `rate = abs(DR)`，`xend`/`y_slope` 同步。**唯一非机械取反处，请出图核对** |

### 10.9 未改动（刻意）

- `workflow/src/coherence/metrics.py` 的 `DR_NORM_MAX` / `normalize_dr_dl`：无调用方的死代码，改它只会误导
- `workflow/scripts/coherence/*` 的归一化数值：见 §4.3
- `tests/test_verification.py` 的 `_make_gene_results` 夹具仍用正的 DR 值——那两个测试只验合并/计数，与符号无关
- `workflow/src/utr/core.py` 的 `um_ratio`：分子分母同时翻号，比值不变
