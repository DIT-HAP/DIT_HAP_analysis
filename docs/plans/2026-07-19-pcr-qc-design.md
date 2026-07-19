# PCR quality control 阶段设计

**日期：** 2026-07-19
**范围：** 中等 — 把 `DIT_HAP_pipeline/workflow/notebooks/thesis_figures.ipynb` 的 "PCR quality control" 部分（建库质控 2×2 四联图）迁入本项目，作为一个新的确定性 Snakemake 阶段 `pcr_qc`，图表统一走本项目 `config/DIT_HAP.mplstyle`。

## 1. 背景与分类修正

原设计文档（`2026-07-15-DIT-HAP-analysis-design.md` §5）把 `thesis_figures` 整体归类为"人工（图表迭代）"的 notebook。但其中 "PCR quality control" 一节（源 notebook cell 4-6）产出一张 2×2 四联图，四个面板全是确定性计算（读固定文件 → 散点 + 线性回归 → 存 PDF），无任何人工判断：

- **(a) PBL vs PBR**：同一建库的两个转座子引物（左/右）读数相关性。
- **(b) 技术重复**：同一物理样本 `LD1328-7_0h_YES` 在两个 upstream project（`LD_DIT_HAP_generationRAW` 与 `Spore2YES6_1328`）里各独立跑一遍，比 `Reads`。
- **(c) 生物学重复**：`LD_DIT_HAP_generationRAW` 里 `LD1328-4` vs `LD1328-8`，比 `Reads`。
- **(d) Spike-in 线性度**：稀释比 vs 读数比的对数线性拟合。

按本项目"确定性计算交给 Snakemake"的原则，(a)(b)(c) 应转为脚本+规则，不留在 notebook。面板 (d) 复用的 `spike_in_results.tsv` 概念上属于独立的 `spikein.smk`（Phase 3+ 尚未实现的 deferred 阶段），本次先用占位表，契约写明将来接口。

## 2. 数据接入 — release 契约之外的例外

其它所有阶段只消费 upstream 的 `release/`（打包后的稳定契约）。PCR QC 例外：它读的是 **pre-release 中间产物** `results/8_merged/{sample}_{timepoint}_{condition}.tsv`（strand-resolved 的 PBL/PBR/Reads 表），`release/` 从不打包这些。

因此需要给数据源注册表和 loader 各开一个受控的例外口子，但仍守住"构造 DIT_HAP_snakemake 路径的地方只有 `data_config.py`"这条原则。

### 2.1 `config/datasets.yaml`

给涉及的两个 dataset 加一个可选字段 `results_dir`（与 `release_dir` 平级，指向 project 的中间 `results/` 目录）。只有需要 PCR QC 的 dataset 才加，语义上明说这是 release 契约之外的用途：

```yaml
datasets:
  LD_DIT_HAP_generationRAW:
    label: "LD, raw generations"
    release_dir: projects/LD_DIT_HAP_generationRAW/release
    results_dir: projects/LD_DIT_HAP_generationRAW/results   # 新增：pre-release 中间产物（仅 PCR QC 用）
    has_time_points: true
    has_imputation: true
  Spore2YES6_1328:
    label: "Spore-to-YES6 1328"
    release_dir: projects/Spore2YES6_1328/release
    results_dir: projects/Spore2YES6_1328/results            # 新增
    has_time_points: true
    has_imputation: false
```

### 2.2 `workflow/src/data_config.py`

新增一个模块常量 + 一个 loader（`8_merged` 这个 upstream 编号目录只在此出现一次，将来上游改名只改一处）：

```python
# pre-release 中间产物子目录（合并 PBL/PBR strand 读数的 merge_strand_insertions 输出）。
# release/ 从不打包这些；只有 PCR QC 这类质控需要，且仅限注册了 results_dir 的 dataset。
MERGED_READS_SUBDIR = "8_merged"

def merged_reads_path(dataset: str, sample: str, timepoint: str, condition: str) -> Path:
    """解析 pre-release 的 merged PBL/PBR 读数表路径（results/8_merged/...）。
    只有在 datasets.yaml 里注册了 results_dir 的 dataset 才能解析；否则 raise KeyError。"""
```

`DatasetConfig` 不改（它描述 release 契约）；`merged_reads_path()` 直接从注册表读 `results_dir`，缺字段就报错，与"这是例外用途"的语义一致。

### 2.3 `config/analysis.yaml`

比哪些样本，是实验参数，进 config（与 clustering/enrichment/ml 一致）：

```yaml
# --- PCR / library-prep quality control (pcr_qc.smk) ---
pcr_qc:
  # 面板 (a): 同一建库 PBL vs PBR 引物读数
  pbl_pbr:
    dataset: LD_DIT_HAP_generationRAW
    sample: LD1328-7
    timepoint: 0h
    condition: YES
  # 面板 (b): 同一样本跨两个 upstream project 的技术重复（比 Reads）
  technical_replicate:
    dataset_1: LD_DIT_HAP_generationRAW
    dataset_2: Spore2YES6_1328
    sample: LD1328-7
    timepoint: 0h
    condition: YES
  # 面板 (c): 同一 project 内两个样本的生物学重复（比 Reads）
  biological_replicate:
    dataset: LD_DIT_HAP_generationRAW
    sample_1: LD1328-4
    sample_2: LD1328-8
    timepoint: 0h
    condition: YES
```

## 3. 代码组织 + 统一绘图风格

### 3.1 `workflow/src/plotting/generic.py`（新建）

设计文档 §7 早已规划的通用图模块（此前未建）。放不含生物学概念的通用图，从 `DIT_HAP_pipeline/workflow/src/plot.py` 移植 `create_scatter_correlation_plot`（log-log 散点 + y=x 对角线 + PCC/R²/slope/RMSE 统计框），函数签名保持不变，供面板 (a)(b)(c) 复用。`donut_chart` 一并移植（未来其它图会用）。

### 3.2 `workflow/scripts/pcr_qc/plot_pcr_qc.py`（新建）

遵循 python-script-conventions（loguru、frozen dataclass config、`@logger.catch(reraise=True)`、类型标注、一行 docstring）：
- 从 CLI 接收三组样本对 + spike-in 表路径 + 输出路径。
- 经 `data_config.merged_reads_path()` 解析 (a)(b)(c) 的输入 `.tsv`。
- 拼 2×2 `subplot_mosaic`，(a)(b)(c) 调 `create_scatter_correlation_plot`，(d) 读 spike-in 表做 `linregress` 线性拟合散点（面板 d 逻辑是本图专属，留在脚本，不进 generic）。
- 存 `results/pcr_qc/PCR_quality_control.pdf`。

### 3.3 统一风格

`from workflow.src.plotting.style import AX_WIDTH, AX_HEIGHT, COLORS` —— 唯一入口应用 `config/DIT_HAP.mplstyle`（Arial + base 字号）。**不**复刻源 notebook cell 0 的 rcParams 覆盖（Times New Roman + 放大字号）。这样 PCR QC 图与本项目其它所有产物（如 ml 的 `prediction_and_residuals.pdf`）风格一致——这正是"本项目图表统一风格"的含义。代价是与旧论文里那张图的字体/字号不同，已确认接受。

## 4. spike-in 占位

把 pipeline 已跑出的 `spike_in_results.tsv`（31 行小文件）复制到 `resources/curated/spike_in_results_PLACEHOLDER.tsv`。脚本 config 与 docstring 的 Input 段同时写明：

```
面板 (d) spike-in 表：
  当前：resources/curated/spike_in_results_PLACEHOLDER.tsv（占位）
  将来：results/spikein/spike_in_results.tsv（由 spikein.smk 产出，Phase 3+）
```

## 5. Snakemake 接线

### 5.1 `workflow/rules/pcr_qc.smk`（新建）

单一 rule `plot_pcr_qc`，**无** dataset wildcard（此 QC 是对几个具名 library 的一次性对比，不是 per-dataset 泛化操作，与 `spikein.smk` 同构）。`input` 用 `data_config.merged_reads_path()` 解析出的具体 `.tsv` + spike-in 占位表；`output` 为 `results/pcr_qc/PCR_quality_control.pdf`；`params` 从 `config["pcr_qc"]` 取样本标识；`conda` 复用 `../envs/statistics_and_figure_plotting.yml`；`log` 到 `logs/pcr_qc/`。

### 5.2 `Snakefile`

- 顶部 `include: "workflow/rules/pcr_qc.smk"`。
- `rule all` 加一条注释掉的 `results/pcr_qc/PCR_quality_control.pdf`（沿用"列出但注释"风格）。

## 6. 测试 `tests/test_pcr_qc.py`（仿 test_clustering.py）

- `merged_reads_path()` 对注册了 `results_dir` 的 dataset 解析正确；对未注册的 dataset raise。
- `create_scatter_correlation_plot`：合成数据返回 Axes；`xscale="log"` 过滤非正值；统计文本框存在。
- 脚本 config dataclass 的 `validate()`：缺输入报错、输出目录被创建。
- 全程用 `tmp_path` 造小 tsv，不依赖真实大文件；`matplotlib.use("Agg")`。

## 7. 验证

- `pytest` 全量回归，确认现有测试仍通过 + 新测试通过。
- `snakemake -n results/pcr_qc/PCR_quality_control.pdf` dry-run 确认 DAG 与 shell 插值正确。
- `snakemake --use-conda results/pcr_qc/PCR_quality_control.pdf` 真跑一次出图验证。

## 8. 不做的事

- 不实现 `spikein.smk`（面板 d 用占位表，Phase 3+ 再补真实来源）。
- 不迁移 `thesis_figures.ipynb` 的其它节（代时计算、insertion distribution、normalization 等）——本次只做 PCR quality control 一节。
- 不改 `DatasetConfig`/`InsertionLevelConfig`/`GeneLevelConfig`（它们描述 release 契约；`results_dir` 走独立的 `merged_reads_path()`）。
- 不复刻源 notebook 的 rcParams 覆盖。
