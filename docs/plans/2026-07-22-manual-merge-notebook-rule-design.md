# 把 manual_merge finalize 笔记本包成可构建 rule（方案 A）

**日期：** 2026-07-22
**范围：** 中等。把 `notebooks/clustering/finalize_gene_clusters.ipynb`（`manual_merge` 变体）
从"人工步骤、Snakemake 不可构建"改成一条真正的 rule：Snakemake 用**原生 `notebook:` 指令**
（仿 `DIT_HAP_pipeline` 的 `extract_genome_region`）无头执行该笔记本，产出
`results/clustering/{dataset}/{variant}/final_clusters.tsv`。同时让"所有 variant 汇总散点图"
把这个变体也画进去。

> **执行机制定案（2026-07-22，用户指定）：** 用 Snakemake 原生 `notebook:` 指令，**不用
> papermill、不写包装脚本**。参考 `DIT_HAP_pipeline/workflow/rules/preparation.smk` 的
> `extract_genome_region`。笔记本直接读注入的 `snakemake` 对象、自写 `final_clusters.tsv`
> + `metrics.tsv`；无 papermill 时 Snakemake 回退到 `jupyter-nbconvert --execute`，故 rule
> 的 conda env 只需加 `nbconvert` + `ipykernel`（不加 papermill）。下面 §2 保留了早期
> papermill 方案的记录，但**已被本决定取代**。

## 1. 决策与后果（用户已选方案 A）

方案 A = 把笔记本本身作为 rule 执行（不改写成 py 脚本）。这意味着笔记本里那份为
HD_DIT_HAP + kmeans 手调的 `merge_groups` dict 被冻结进 DAG——这是方案 A 固有的取舍，
用户已接受。因为 kmeans 用 `random_state=42` 确定性聚类，重跑可复现同一结果，冻结的只有
`merge_groups` 这一份人工判断（仍留在笔记本里，可人工编辑）。

**核心决策：`manual_merge` 升级为"可构建变体"（builder = 笔记本）。** 一旦它有了 rule，
就应当和其它 buildable 变体走同一套路径/汇总逻辑，否则 `final_clusters_path()` 会对同一个
变体给出两个不同路径（curated vs results），非常脆弱。因此：

- `final_clusters_path(manual_merge)` → `results/clustering/{dataset}/{variant}/final_clusters.tsv`
  （与用户当前笔记本写出路径一致）。
- `buildable_variants()` 纳入 manual_merge（此后 == 全部变体）。
- 汇总散点图 / `compare_variants` / enrichment / ml 全部自动统一识别它。

**唯一行为变化（需知晓）：** 默认 `snakemake`（rule all 含
`all_variants_cluster_scatter.pdf` 和 `variant_metrics_comparison.tsv`）此后会**自动执行
该笔记本**。这正是方案 A 的意图（笔记本进入流水线），但和过去"rule all 绝不自动跑人工步骤"
的旧不变量相反——下面文档同步更新。

> 备选（更保守，不推荐）：保留 `buildable_variants()` 只含 direct/auto_merge/grid，另设
> `all_plottable_variants()=buildable+manual` 只喂两个散点图 rule；enrichment/compare 不变。
> 缺点：同一变体两套路径、两个概念，不一致。默认取上面的统一方案。

## 2. 执行机制：Snakemake 原生 `notebook:` 指令（定案）

rule 用 `notebook: "../../notebooks/clustering/finalize_gene_clusters.ipynb"`。Snakemake 在
笔记本首部注入一个 preamble cell：`os.chdir(workdir)`（= 仓库根）+ 把仓库根加进 `sys.path`
+ 定义 `snakemake` 对象（`.input/.output/.params/.wildcards/...`）。执行器优先 papermill，
缺则回退 `jupyter-nbconvert --execute`；执行后的副本写到 `log: notebook=...`。

- **无包装脚本、无 papermill 依赖。** 笔记本本身即 builder，读 `snakemake` 对象、自写两个
  产物。`metrics.tsv` 由笔记本末尾 cell 调 `score_labels` 写出（与 buildable 脚本同一函数）。
- **双模可用。** 笔记本用 `try: snakemake / except NameError` 判定运行模式：无头态读注入对象，
  交互态回退到内置默认值。人在 Jupyter 里打开仍能正常跑、看 review 图、调 `merge_groups`。

### cwd / import 陷阱（已解决）
笔记本原用 `sys.path.insert(0, Path.cwd().parents[1])` 且 `REPO=Path.cwd().parents[1]`，
依赖 cwd 恰为 `notebooks/clustering/`。原生 `notebook:` 执行时 preamble 已 `chdir` 到仓库根，
故无头态直接 `REPO = Path.cwd()`；交互态才回退 `Path.cwd().parents[1]`。import cell 里用
`try: snakemake` 分支选对 REPO 后再 `sys.path.insert`，两模式都正确。

## 3. 具体改动清单

### 3.1 笔记本 `notebooks/clustering/finalize_gene_clusters.ipynb`
双模重排，保持逻辑与 `merge_groups`（用户当前工作区版本）不变：
- imports cell：`try: snakemake; REPO = Path.cwd()` / `except NameError: REPO = Path.cwd().parents[1]`，
  再 `sys.path.insert(0, REPO)` + import（额外 import `score_labels`）。
- config cell：`try: snakemake` 判定 `_SM`。无头态从 `snakemake.wildcards/params/input/output`
  取 `DATASET/VARIANT/METHOD/N_INTERMEDIATE/WT_CLUSTER/FINAL_N_CLUSTERS/RANDOM_STATE` 与
  `ANNOTATED_PATH/SCALED_PATH/OUTPUT/METRICS_OUTPUT`；交互态回退到内置默认值 + `results/` 路径。
  之后读数据 + `cluster_one_method` 聚类（原逻辑）。
- 写出 cell：写 `final_clusters.tsv`，再算 `score_labels`（`variant_type="manual_merge"`）写
  `metrics.tsv`。
- review 图 / `merge_groups` cell 不变。header markdown 更新为双模说明。
- 不打 `parameters` tag（原生 `notebook:` 不需要；靠 `snakemake` 对象注入）。

### 3.2 `workflow/rules/clustering.smk`
- `final_clusters_path()`：所有 type（含 manual_merge）统一走 `results/.../final_clusters.tsv`
  （抽出 `_FINALIZE_TYPES` 常量）。
- `buildable_variants()`：纳入 manual_merge == 全部变体（更新 docstring）。
- 新 rule `finalize_manual_merge`：`wildcard_constraints: variant=_alt(_variants_of_type("manual_merge"))`；
  input = `_work/{annotated,scaled}.parquet`；output = final_clusters.tsv + metrics.tsv；
  params = method/n_intermediate/final_n_clusters/random_state/wt_cluster；
  `conda: statistics_and_figure_plotting.yml`；`log: notebook=...ipynb`（执行后副本）；
  `notebook: "../../notebooks/clustering/finalize_gene_clusters.ipynb"`。
- 文件头注释更新。`plot_variant_clusters` / `plot_all_variants*` / `compare_variants` 无需改
  （已走 `buildable_variants()` / 全局 variant 约束，自动含 manual）。

### 3.3 env `workflow/envs/statistics_and_figure_plotting.yml`
加 `nbconvert`、`ipykernel`（原生 `notebook:` 无 papermill 时回退 `jupyter-nbconvert --execute`
所需；该 env 已含笔记本用到的 matplotlib/seaborn/sklearn/scipy/pandas/pyarrow）。**不加 papermill。**

### 3.4 配置 & 文档（把"manual_merge 不可构建"这条不变量同步掉）
- `config/analysis.yaml`：改 manual_merge 那行注释。
- `docs/plans/2026-07-21-clustering-finalize-variants-design.md`：§2 表格、§3.4、§5、§8 里
  manual 不可构建的描述改为"由 `finalize_manual_merge` rule 无头执行笔记本构建"。
- `Snakefile`、`workflow/rules/enrichment.smk` 头部关于 manual "missing input by design" 的
  注释更新。
- `CLAUDE.md`（git-ignored）：核心链图 + 目录契约里 manual/curated 描述更新。

### 3.5 测试 `tests/test_clustering.py`
无新增脚本 dataclass（native `notebook:` 无包装脚本），故不加 `*Config` 单测。笔记本全流程
执行不做单测（需 kernel + nbconvert）；靠下面的实跑验证。现有 clustering 单测须保持全绿。

## 4. 验证
1. `pytest`：应保持全绿（无回归）。
2. `snakemake -n results/clustering/HD_DIT_HAP/kmeans_merge9_manual/final_clusters.tsv` 和
   `.../all_variants_cluster_scatter.pdf`：确认 DAG 里出现 `finalize_manual_merge`，且汇总
   散点图输入含 manual。
3. 实跑（`--use-conda` 构建含 nbconvert/ipykernel 的 env，`_work/*.parquet` 已在盘上）：
   `snakemake --use-conda --cores 4 results/clustering/HD_DIT_HAP/kmeans_merge9_manual/final_clusters.tsv`
   → 确认产出 final_clusters.tsv（9 簇、`cluster` 1..9、WT=9 最低 DR、含 `raw_cluster`）+
   metrics.tsv + 执行后笔记本；再跑
   `.../all_variants_cluster_scatter.pdf` → 确认 8 个 variant 面板（含 kmeans_merge9_manual）。
