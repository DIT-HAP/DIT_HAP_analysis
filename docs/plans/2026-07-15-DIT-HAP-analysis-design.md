# DIT-HAP Analysis 项目设计

**日期**: 2026-07-15
**范围**: 新建独立项目 `DIT_HAP_analysis`，重新组织 `DIT_HAP_pipeline/workflow/notebooks/` 中的下游分析（enrichment、clustering、ML、thesis 图表等）。
**状态**: 已确认，待实施

---

## 1. 背景与目标

`DIT_HAP_pipeline/workflow/notebooks/` 下现有 17 个笔记本承担全部下游分析（特征收集、聚类、富集分析、ML、结构域差异、thesis 图表等），存在三个问题：

1. **确定性计算与人工判断混在同一个笔记本里**，例如 `gene_level_clustering.ipynb` 既有可重跑的聚类算法，也有手动合并簇的人工决策；再重跑往往会不小心冲掉人工决策的结果。
2. **笔记本之间的依赖是隐式的**——只能通过读代码里硬编码的路径字符串才能发现"这个笔记本依赖另一个笔记本的产出"，例如 `further_analysis_based_on_enrichment.ipynb` 依赖 `comprehensive_enrichment_analysis.ipynb` 产出后经人工/LLM 整理的 xlsx。
3. **共享代码模块职责不清**——`workflow/src/utils.py` 混装不相关函数，`enrichment_functions.py`、`pombe_feature_functions.py` 各自变成大杂烧模块。

同时，上游流水线已完成 `DIT_HAP_snakemake` 的重构，采用 `projects/{project_name}/release/` 作为下游消费的稳定契约（`insertion_level/`、`gene_level/` 下的固定文件名）；`DIT_HAP_streamlit` 已经用一套「YAML 注册表 + frozen dataclass 路径配置 + 缓存加载器」的模式消费这份契约。新项目从这个模式借鉴思路（配置驱动、类型化路径），但作为独立代码库实现，不与 streamlit 项目共享代码——两者的数据形态、分析方式和运行环境不同，硬耦合会增加双方维护成本。

**目标**：产出一个新的、结构清晰、职责边界明确、易于增量修改的下游分析项目 `DIT_HAP_analysis`，将确定性计算交给 Snakemake 编排，需要人工判断或反复迭代的分析保留为笔记本，两者通过明确的输入输出契约衔接。

---

## 2. 顶层结构

```
DIT_HAP_analysis/                  # 新建独立 git 仓库，与 DIT_HAP、DIT_HAP_snakemake、DIT_HAP_streamlit 平级
├── Snakefile                      # 编排确定性分析阶段
├── config/
│   ├── datasets.yaml              # 数据源注册表 — 指向各 DIT_HAP_snakemake project 的 release/ 绝对路径
│   ├── analysis.yaml              # 本项目自己的分析参数（聚类 k 值范围、富集阈值等）
│   └── DIT_HAP.mplstyle           # 复用绘图样式
├── workflow/
│   ├── rules/                     # .smk 规则模块，按分析阶段分组
│   ├── scripts/                   # 确定性阶段的脚本，遵循 python-script-conventions
│   └── src/                       # 共享库：数据读取、绘图、富集统计、特征工程
├── notebooks/                     # 人工判断/迭代类分析，明确写清读取哪个 results/ 阶段的产出
├── results/{stage}/{dataset}/     # Snakemake 产出，语义化命名而非笔记本编号
├── reports/                       # 图表/HTML 报告
└── resources/                     # 外部参考数据库副本 + 本项目人工整理产出
```

与现状的关键差异：
- `results/` 不再是 `18_gene_level_clustering/`、`19_enrichment_analysis/` 这种编号目录，改成语义化命名（`results/clustering/`、`results/enrichment/`），阶段间的依赖关系从目录名就能看出来。
- 不存在物理的 `data/` 目录；数据源通过配置注册表引用，见第 3 节。
- 笔记本按"是否可自动重跑"一分为二：能进 DAG 的转成脚本，不能的留在 `notebooks/`，但输入路径写死指向 `results/{stage}/`，减少隐式依赖链。

---

## 3. 数据源注册表：`config/datasets.yaml`

**要解决的问题**：现状 17 个笔记本到处硬编码 `Path("../../results/HD_DIT_HAP_generationRAW/17_gene_level_curve_fitting/...")`，上游一旦调整路径规范就要挨个笔记本改。注册表把"数据在哪"和"分析怎么做"彻底分开。

**决策**：不拷贝、不建软链接，配置注册表直接引用 `DIT_HAP_snakemake` 里各 project 的 `release/` 绝对路径。两个仓库需在同一台机器、路径相对固定的前提下工作。

```yaml
# config/datasets.yaml
default_dataset: HD_DIT_HAP_generationRAW

snakemake_repo: /data/c/yangyusheng_optimized/DIT_HAP_snakemake   # 上游仓库根，唯一需要改的路径

datasets:
  HD_DIT_HAP_generationRAW:
    label: "HD, raw generations"
    release_dir: projects/HD_DIT_HAP_generationRAW/release   # 相对 snakemake_repo
    has_time_points: true
  LD_DIT_HAP_generationRAW:
    label: "LD, raw generations"
    release_dir: projects/LD_DIT_HAP_generationRAW/release
    has_time_points: true
  Spikein:
    label: "Spike-in calibration"
    release_dir: projects/Spikein/release
    has_time_points: false
  # ... 其余项目同构（HD_DIT_HAP, HD_DIT_HAP_generationPLUS1, LD_DIT_HAP_generationPLUS1,
  #     HD_diploid, LD_haploid, Spore2YES6_1328）

reference:
  pombase_version: "2025-10-01"
```

对应 `workflow/src/data_config.py`：frozen dataclass 描述每类文件的路径集合（`GeneLevelConfig`、`InsertionLevelConfig` 等），`load_dataset_config(name)` 读 YAML 并拼出绝对路径。**这是本项目读取外部仓库数据的唯一入口**——除它之外任何脚本/笔记本都不应该自己拼 `DIT_HAP_snakemake` 的路径。

`has_time_points` 对应上游 `packaging.smk` 里的同名分支逻辑（有些 project 如 `Spikein` 没有 `time_points`，release 里就没有 gene-level curve fitting 产出）；本项目消费时据此判断某个 dataset 能否进入需要 gene-level 曲线拟合结果的分析阶段。

---

## 4. `resources/` 布局

分两层：外部数据库的本地副本（可重新下载/更新）与本项目人工整理产出（版本控制，不可自动重现）。

**决策**：PomBase、STRING、KEGG、BioGrid 等外部数据库，本项目自己完整保留一份，不依赖 `DIT_HAP_snakemake` 里的版本，避免跟上游数据更新耦合。AlphaFold 结构数据体量大，仍引用仓库外部路径，不拷贝。

```
resources/
├── external/                          # 外部数据库本地副本，脚本可重新下载/更新
│   ├── pombase/{release_version}/     # 基因组、注释、ontology（GO/FYPO/MONDO）、gene metadata
│   ├── stringdb/                      # protein_interactions, protein_information, evolutionary_data
│   ├── kegg/                          # brite_json, brite_table, pathways
│   ├── biogrid/                       # PPI/GI tab3 文件
│   └── ensembl/                       # paralog 数据
├── literature/                        # 文献补充数据表（Marguerat 2012、Harigaya 2016 等）
└── curated/                           # 人工整理/校对产出，本项目独有，版本控制
    ├── essentiality_verification.csv           # 人工校验的必需性汇总
    ├── deletion_library_categories.xlsx        # Hayles 2013 整理表
    ├── enrichment_term_categorization.xlsx     # GO term 富集结果的人工/LLM 分类整理
    ├── final_clusters.tsv                      # 人工确认后的基因聚类结果（枢纽产物，见第 6 节）
    └── non_essential_domain_candidates.xlsx    # 结构域差异候选人工审核表
```

AlphaFold 数据保持引用仓库外部路径（如 `../../../resource/AlphaFold_Dataset/...`），由读取该数据的脚本/笔记本在 config 里声明该路径，不进 `resources/`。

外部数据库下载脚本放 `workflow/scripts/reference/fetch_*.sh`，直接搬用 `DIT_HAP_snakemake` 里已有的 PomBase/STRING/KEGG 下载脚本并改路径，作为 Snakemake rule 的输入准备步骤，可重跑更新版本。

---

## 5. 笔记本分类：确定性脚本 vs 人工笔记本

**编排原则**：Snakemake 编排确定性步骤，需要人工判断或反复迭代调参的步骤留作独立 notebook；两者通过第 6 节的输入输出契约衔接，不强制全部串成一条 Snakemake DAG。

| 现有笔记本 | 归类 | 新项目里的形态 |
|---|---|---|
| `pombe_feature_collection` | 确定性 | 脚本 → `results/features/{pombase_version}/` |
| `comprehensive_enrichment_analysis` | 确定性 | 脚本 → `results/enrichment/raw/{dataset}/` |
| `machine_learning_data_preparation` | 确定性 | 脚本 → `results/ml/features_targets/{dataset}/` |
| `machine_learning_analysis` | 确定性（AutoML 训练本身无需人工介入） | 脚本 → `results/ml/models/{dataset}/` |
| `gene_coverage_analysis` | 确定性 | 脚本 → `results/coverage/{dataset}/` |
| `compare_with_deletion_library` | 确定性 | 脚本 → `results/verification/{dataset}/` |
| `compare_with_other_large_scale_studies` | 确定性 | 脚本 → `results/comparison/{dataset}/` |
| `non_coding_RNA_analysis` | 确定性 | 脚本 → `results/noncoding_rna/{dataset}/` |
| `spike_in` | 确定性 | 脚本 → `results/spikein/` |
| `complex_analysis` | 确定性 | 脚本 → `results/complex/{dataset}/` |
| `gene_level_clustering` | **拆分**：聚类算法与评估指标确定性；"手动合并簇"是人工决策 | 脚本产出候选聚类 → `results/clustering/candidates/{dataset}/`；人工合并留作 notebook，产出 `resources/curated/final_clusters.tsv` |
| `upstream_and_downstream_analysis` | **拆分**：UTR 归属统计确定性；图表反复调整属人工 | 脚本产出统计表 → `results/utr/{dataset}/`；notebook 做可视化精修 |
| `genes_with_domain_differences` | **拆分**：核心统计确定性，依赖人工审核的候选表 | 脚本 → `results/domain_differences/{dataset}/` |
| `further_analysis_based_on_enrichment` | 人工（依赖 LLM/人工分类整理的 xlsx） | notebook，读 `results/enrichment/raw/` + `resources/curated/enrichment_term_categorization.xlsx` |
| `thesis_figures` | 人工（图表迭代） | notebook，读各 `results/{stage}/` |
| `visualize_dit_hap_style` | 非分析，样式演示 | 移入 `notebooks/reference/`，不算分析产物 |
| `gff_processing_and_annotation` | 不在本项目范畴 | 已在 `DIT_HAP_snakemake` 重构计划中转为上游脚本 |

`gene_level_clustering` 拆分后的枢纽产物 `resources/curated/final_clusters.tsv`（现状 `kmeans_cluster_result.tsv`，被 7+ 个下游笔记本消费）被显式放进 `resources/curated/`（版本控制）而非 `results/`（可重跑覆盖）——因为它包含不可自动重现的人工判断，语义上是"已审核的输入"而非"计算的输出"。下游脚本读它时能清楚意识到这是需要人来更新的文件。

---

## 6. 人工笔记本的输入输出契约

`notebooks/` 按分析主题分目录（不用原笔记本名，因为同一主题常会拆出"最终确认"和"深挖探索"两类笔记本）：

```
notebooks/
├── clustering/
│   └── finalize_gene_clusters.ipynb       # 读 results/clustering/candidates/ → 写 resources/curated/final_clusters.tsv
├── enrichment/
│   ├── categorize_enrichment_terms.ipynb  # 读 results/enrichment/raw/ → 写 resources/curated/enrichment_term_categorization.xlsx
│   └── deep_dive_enriched_terms.ipynb     # 读 results/enrichment/raw/ + curated 分类表 → 图表
├── domain_analysis/
│   └── review_domain_differences.ipynb    # 读 results/domain_differences/ → 图表/审核
├── thesis_figures/
│   └── figures.ipynb                      # 读多个 results/{stage}/ → reports/thesis/*.pdf
└── reference/
    └── visualize_dit_hap_style.ipynb      # 样式演示，非分析
```

每个笔记本第一个 markdown cell 固定格式，声明输入输出契约：

```markdown
## Inputs
- results/enrichment/raw/{dataset}/go_enrichment_full_filtered.tsv  (from: workflow/scripts/enrichment/run_go_enrichment.py)
- resources/curated/final_clusters.tsv  (manually curated, see notebooks/clustering/finalize_gene_clusters.ipynb)

## Outputs
- resources/curated/enrichment_term_categorization.xlsx  (manually curated — consumed by deep_dive_enriched_terms.ipynb)
```

第二个 cell 才是实际的 `Config` dataclass，路径统一来自 `workflow/src/data_config.py` 暴露的常量/函数，不再手写字符串路径。这样任何人打开笔记本第一屏就能定位它在数据链条里的位置。

---

## 7. `workflow/src/` 重新分组

现状问题：`utils.py` 混装两个不相关函数；`enrichment_functions.py`、`pombe_feature_functions.py` 各自是大杂烧模块，配置类、计算函数、IO 混在一个文件；`protein_structure_functions.py` 只有一个函数，独立成模块意义不大；绘图函数按"谁调用它"而非"做什么"分散在 `plot.py` 和 `subset_visualization.py` 两处。

按职责重新分组，通用逻辑与领域特定逻辑分离：

```
workflow/src/
├── data_config.py          # 数据源注册表加载（第 3 节）
├── io.py                    # 通用文件读取（read_file）、表格拼接（原 concat_tables.py）
├── gene_ids.py               # 基因 ID/名称/同义词映射（原 utils.py::update_sysIDs）
├── plotting/
│   ├── style.py               # mplstyle 加载、AX_WIDTH/AX_HEIGHT/COLORS 常量（原 plot.py 顶部）
│   ├── generic.py              # 通用图表：散点相关图、donut chart（原 plot.py 剩余部分）
│   └── gene_level.py            # 基因层面专用图：depletion 曲线、feature space 散点、基因内插入可视化
│                                 #（合并原 subset_visualization.py + intragenic_insertion_visualization.py）
├── enrichment/
│   ├── ontology.py             # ontology 数据加载与配置（原 enrichment_functions.py 中 OntologyDataConfig/GeneMetaData 等）
│   └── pipeline.py              # 富集计算与富集图（ontology_enrichment_pipeline、stringdb_enrichment、revigo_analysis、customized_enrichment_plot）
└── features/
    ├── genome.py                # 基因组坐标常量与 DNA 层特征（原 pombe_feature_functions.py 一部分）
    └── protein.py                 # 蛋白特征提取、codon usage、pLDDT 统计（合并原 pombe_feature_functions.py 剩余部分 + protein_structure_functions.py）
```

**原则**：`plotting/generic.py` 里的函数不能出现任何生物学专属概念（如 gene systematic ID），保证独立于本项目复用；`gene_level.py` 反过来完全允许硬编码基因组学假设。`enrichment/` 与 `features/` 各自是独立子包，互不引用，避免隐式耦合（现状 `enrichment_functions.py` 意外依赖 `pombe_feature_functions.py`）。

---

## 8. `Snakefile` 与 `workflow/rules/` 规则设计

**Wildcard 设计**：`dataset` 是核心 wildcard，取值约束为 `config/datasets.yaml` 中注册的 project 名。规则分两类：

- **per-dataset 规则**（clustering、utr、verification、coverage、spikein、domain_differences）：输出路径带 `results/{stage}/{dataset}/...`，因为不同 project（HD/LD/Spikein 等）算出来的东西天然分开。
- **dataset 独立规则**（features 等）：`pombe_feature_collection` 算的是基因的静态生物学特征，跟哪个测序项目无关，只跟 PomBase 版本有关，输出 `results/features/{pombase_version}/...`，不带 `dataset` wildcard。

**人工产物作为"不可构建"的输入**：`resources/curated/final_clusters.tsv` 等人工确认产物，Snakemake 规则直接当普通输入文件引用，不写生成它的 rule。如果文件不存在，DAG 会报"missing input, no rule to produce it"——这是有意为之，提醒需要先跑对应笔记本、人工确认后再继续。

```
workflow/rules/
├── features.smk          # pombe_feature_collection → results/features/{pombase_version}/
├── clustering.smk         # gene_level_clustering(确定性部分) → results/clustering/candidates/{dataset}/
├── enrichment.smk          # comprehensive_enrichment(读 resources/curated/final_clusters.tsv) → results/enrichment/raw/{dataset}/
├── ml.smk                  # data_preparation + analysis → results/ml/{features_targets,models}/{dataset}/
├── coverage.smk
├── verification.smk
├── comparison.smk
├── noncoding_rna.smk
├── spikein.smk
├── complex.smk
├── utr.smk
└── domain_differences.smk
```

`Snakefile` 顶部通过 `yaml.safe_load` 直接读取 `config/datasets.yaml`（不是通过 `--configfile` 传入，因为这是数据注册表而非实验参数），校验 `dataset` wildcard 只能取注册表中的 key。`rule all` 沿用现有仓库风格——列出所有 stage 的最终目标但注释掉，按需取消注释跑指定阶段。

跨阶段依赖示例（体现显式 DAG 而非隐式笔记本链）：

```python
rule comprehensive_enrichment:
    input:
        clusters = "resources/curated/final_clusters.tsv",       # 人工产物，非本项目构建
        gene_ontology = lambda wc: reference_paths().gene_ontology_gaf,
    output:
        "results/enrichment/raw/{dataset}/go_enrichment_full_filtered.tsv"
```

---

## 9. 未决事项 / 后续实施注意点

- 各 `.smk` 规则的脚本转换（笔记本 → `workflow/scripts/{stage}/*.py`）需逐个笔记本核对输入输出，遵循 `python-script-conventions`（loguru、`@dataclass`、`@logger.catch`、类型标注、一行 docstring）。
- `config/analysis.yaml` 的具体字段（聚类 k 值范围、富集 FDR 阈值等）待实施阶段从各笔记本 Config 类中提取汇总。
- 实施建议按阶段拆分为多个 PR/commit：先搭骨架 + 数据源注册表，再逐个迁移确定性阶段为脚本，最后整理人工笔记本的输入输出声明。
