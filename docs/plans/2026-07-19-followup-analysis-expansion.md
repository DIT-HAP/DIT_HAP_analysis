# DIT-HAP 下游分析扩展计划

**日期**: 2026-07-19
**分支/worktree**: `worktree-followup-analysis-expansion`
**范围**: 把「实验验证结果进一步分析」+「更多干实验分析」的关键词笔记，扩充成可执行的分析任务清单。
**状态**: §1–§10 为全景任务清单（7 主题）；**§A-D 为经头脑风暴确认的深化设计（主题 A + D，优先落地）**。

---

## 决策摘要（2026-07-19 头脑风暴）

经与用户逐项确认：
- **优先主题**：A（实验验证深挖）+ D（coherence 归因）。其余主题（B/C/E/F/G）保留在 §3–§8 作为后续 backlog。
- **功能分组源**：GO cellular_component 复合体 + GO biological_process + KEGG pathway（三源共用，A 与 D 共享）。
- **产出形态**：先在 notebook 探索出科学结论/例子，稳定后再**选性**下沉为 scripts + `.smk` 规则。
- **A×D 交汇**：走「纠正成员 essentiality 判定 + 重算 coherence」的闭环。
- **major/minor 文献先验**：手工 curated 表（`resources/curated/complex_subunit_roles.xlsx`）+ 数据驱动 GMM 子聚类对照。

**可行性已验证**（用 `data_analysis` env 实测）：验证集 406 基因中 217 个（53%）是 complex 成员，覆盖 160/477 个 complex；97 个 flip 基因分布在 58 个 complex；flip 最密集的 complex（90S preribosome / SNARE / small-subunit processome，各 3 个）恰好也是 coherence 极值复合体——闭环抓手具体。

---

## 0. 本文定位

这是一份**分析任务规划文档**，不是代码实现方案。它把你笔记里的关键词/半成品任务，逐条扩充成：
「用户原始问题 → 现状（已有什么数据/代码）→ 具体子任务（带数据路径与方法）→ 预期产出 → 优先级/工作量」。

执行时，每个分析主题在 `notebooks/<主题>/` 下新建探索性 notebook（读 `results/{stage}/` 与 `resources/curated/`），确定性、可复跑的部分再下沉为 `workflow/scripts/{stage}/*.py` + `.smk` 规则。这与项目既有的「Snakemake 编排确定性步骤 + notebook 保留人工判断」约定一致（见 `docs/plans/2026-07-15-DIT-HAP-analysis-design.md` §5–6）。

---

## 1. 共享前提：数据位置、术语、环境

**术语**（后文统一用）：
- `DR` = depletion rate（耗竭速率，旧代码里叫 `um`）；`DL` = depletion level（耗竭幅度，旧代码里叫 `lam`）；`A` = 曲线 amplitude。三者是 Gompertz/曲线拟合的 gene-level 参数。
- WT/背景簇 = cluster **9**。
- 上游有一套平行的 gRNA 筛选（`gRNA_*` 列），常用来和 DIT-HAP 相互印证。

**数据来源**（你笔记里给的都是 Mac 本地路径，服务器上真实位置如下）：

| 笔记里的文件 | 服务器真实路径 |
|---|---|
| `organized_verification_summary.xlsx` | `DIT_HAP_pipeline/results/HD_DIT_HAP_generationRAW/20_essentiality_verification/organized_verification_summary.xlsx` |
| `all_coding_genes_with_DIT_HAP_clustering.tsv` | `.../18_gene_level_clustering/all_coding_genes_with_DIT_HAP_clustering.tsv` |
| `go_enrichment_categorized_merged_..._GLM.xlsx` | `.../19_enrichment_analysis/go_enrichment_categorized_merged_gemini_and_claude_and_GLM.xlsx` |
| `complex_coherence_metrics.csv` | `.../25_complex_analysis/complex_coherence_metrics.csv` |

本项目消费上游数据的**唯一正规入口**是 `workflow/src/data_config.py`（读 `config/datasets.yaml` 拼 `release/` 路径）。但上面这些是 `DIT_HAP_pipeline` 旧结果树的产物，不在新 `release/` 契约里——短期直接按上述绝对路径读；中期若某个分析要进 DAG，应改为从 `release/gene_level/gene_level_fitting_statistics.tsv`（列含 `DR`/`DL`/`A`/`R2`）重算，而不是读旧结果。

**gene-level 特征矩阵**（ML/coherence/paralog 都要用）：
`results/features/2026-06-01/pombe_coding_gene_protein_features.tsv`，5126 基因 × 102 列。**注意它没有 GO term ID / complex 成员列**——凡是要按 pathway/complex 分组的分析，都得另外 join PomBase ontology / STRING / KEGG。

**环境**（worktree 内已把 `resources/external`、`resources/literature`、`results` 软链到主树，复用已有下载）：
- 探索用 `data_analysis` conda env（`/data/a/yangyusheng/miniforge3/envs/data_analysis`）：pandas 2.2 / numpy 1.26 / sklearn 1.6 / openpyxl 3.1，能读旧 numpy1 pickle 和 xlsx。
- ML 可解释性任务需要额外装 `mljar-supervised` + `shap`（当前所有 env 都没有）——建议新建一个专用 env（见 §6 任务 ML-0）。
- **Snakefile 有硬编码 `workdir: /data/c/.../DIT_HAP_analysis`**：在 worktree 里跑 `snakemake` 仍会作用于**主树**代码。要验证 worktree 里的改动，直接用 `data_analysis` 的 python 跑 `workflow/scripts/**/*.py`，别经 snakemake。

**Phase 状态**：核心链 clustering→enrichment→ml 已迁移（Phase 2 done）；complex / coherence / noncoding / utr / verification 的脚本化在另一个锁定的 worktree `migrate-remaining-notebooks`（Phase 3）里进行中。**本计划是在这些 stage 之上做新分析，不重复迁移工作**——如果某分析依赖的 stage 还没迁完，就先在 notebook 里直接读旧 pipeline 产物起步。

---

## 2. 主题 A — 实验验证结果的进一步分析

**数据**：`organized_verification_summary.xlsx`（1 sheet，411 行 / 406 唯一基因，31 列）。关键列：
`Systematic ID, Name, FYPOviability, DeletionLibrary_essentiality, Category, DR, DL, R2, RMSE, verification_phenotype, verification_essentiality, comments, total_colonies, filtered_count, filtering_fraction, deletion_with_area_count/fraction, median_area_day3..6, mean_area_day3..6`。
另有已 curated 的 `resources/curated/essentiality_verification.csv`（411 行 × 20 列，含 colony-area 原始测量）。

**关键现状**：
- **没有显式的「正确/失败/差异」列**——outcome 是靠 `verification_essentiality`（最终判定：E=169, V=242）对比 `DeletionLibrary_essentiality` 和 DIT-HAP `Category` 推出来的。
- 与缺失文库的一致性：E∩E=87，V∩V=220，**不一致 = 97**（其中「文库记 V → 我们验证 E」80 个，「文库记 E → 我们验证 V」17 个）。DR 能干净地区分：验证-E 平均 DR≈0.81 vs 验证-V≈0.53。
- `comments`（81 条非空，中文）是目前最接近「数据质量存疑/边界」的信号；cell 19 硬编码了 18 个 `noised_genes`（数据质量排除）。
- **文件里完全没有 pathway/complex 注释**——4 个问题全都需要先补一个 gene→complex/pathway 的 join。

### 任务清单

**VER-1｜给验证结果打上明确的 outcome 标签（前置任务，其他都依赖它）**
- 在 notebook 里派生一列 `outcome ∈ {concordant_E, concordant_V, flip_V→E, flip_E→V, ambiguous}`，规则：`verification_essentiality` vs `DeletionLibrary_essentiality`；`ambiguous` 由 `comments` 关键词 + `noised_genes` + 低 `total_colonies`/`filtering_fraction`/低 `R2` 触发。
- 产出：`results/verification/HD_DIT_HAP/verification_outcome_table.tsv`（在 411 行上加 outcome + 质控 flag）。这是主题 A 所有后续分析的基础表。

**VER-2｜正确验证的基因落在哪些 pathway/complex，是否改变对该 complex 的认识（对应问题①）**
- Join 一个 gene→complex 映射：PomBase GO cellular_component macrocomplex（`resources/external/pombase/.../macromolecular_complex_annotation.tsv`，476 term / 1961 基因）+ 可选 KEGG BRITE / STRING cluster。
- 对每个 complex 统计其成员在验证集里的 outcome 分布；重点看「整个 complex 从文库的 V 被我们改判为 E」或反之的 complex。
- 结合主题 C 的 coherence：如果某 complex 之前被判 incoherent，是不是因为个别成员的文库判定错了，验证纠正后 coherence 会不会变好。
- 产出：`complex × outcome` 汇总表 + 每个「认识被改变」的 complex 的一段文字解读。

**VER-3｜未研究基因里我们有可靠新结果的（对应问题②）**
- 筛选：`gene_name` 仍是 Systematic-ID 形式（如 `SPCC594.04c`）的未命名/未表征基因 ∩ 高 `R2` ∩ `comments` 干净 ∩ outcome 明确。
- 「文库 V → 验证 E」的 80 个是最强的新颖性候选集，优先在里面找未表征基因。
- 产出：候选新基因表，按 `R2`、colony-area 一致性、DIT-HAP/gRNA 双平台一致性排序。

**VER-4｜仍验证失败的：是我们的数据问题还是别的原因（对应问题③）**
- 对 `ambiguous`/矛盾的基因，用 `filtering_fraction`、`deletion_with_area_fraction`、`total_colonies`（低 = 我们数据覆盖不足）、`R2`/`RMSE`（拟合差）逐个归因；对照 18 个 `noised_genes`。
- 生物学解释走 paralog 补偿路线（旧 notebook cell 26 已发现：paralog 高表达代偿 → 单敲除无表型），可复用主题 B 的 paralog 数据。
- 产出：失败基因的归因分类（数据质量 / 生物学冗余 / 条件依赖 / 其他）。

**VER-5｜差异如何解释（对应问题④，97 个不一致基因）**
- 交叉 `FYPOviability`（41 个 condition-dependent，可能解释条件依赖差异）、DR/DL、colony-area day3–6 曲线形状。
- 用 DIT-HAP vs gRNA 双平台 depletion 曲线（旧 notebook cells 34–41）确认我们的判定站得住。
- 产出：97 个差异基因的逐条解释表 + 汇总类别。

**优先级**：VER-1（必做前置）→ VER-2/VER-5（科学价值最高）→ VER-3 → VER-4。
**工作量**：VER-1 小；VER-2 中（要建 complex join）；VER-3/4/5 中。

---

## 3. 主题 B — 聚类分析（更少 bias）

<!-- PLACEHOLDER-CLUSTERING：聚类方法学调查 agent 仍在跑，回来后补全 -->
*（此节待聚类方法学调查完成后补全：现有 4 算法 + k=64→9 手动合并的细节、DR-DL 点云形状与网格边界建议、k-sweep/silhouette 是否已有、以及对「k 直接设 7–9」「DR-DL 网格分格 A1/A2/B1」「连续数据二维聚类是否合适、有无更好方法」三个想法的具体评估。）*

---

## 4. 主题 C — Paralog 系统化分析

**你的问题**：paralog 分析能不能更系统，比如分析所有有 paralog 的基因？

**现状**：
- 特征矩阵有 `paralog_count` 列（5126 基因中 **3738 个 ≥1 paralog，2535 个恰好 1 个**，最多 24）。
- 原始 paralog 对：`resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv`（15721 行，含 `homology_type` = `within_species_paralog` 与 %identity）。第二来源：`PMID20473289/paralogs_from_PMID20473289.csv`。
- 已有一个**旧的、孤立的** notebook `Feature_organization/notebooks/paralog_analysis.ipynb`：按 `Paralog` 列建 group，画每组 Gompertz 曲线并按相对 mRNA 丰度着色。但它用的是旧拟合文件（`fitting_v7_36_params.csv`），**没算 paralog 间 DR-DL 距离，也没做冗余分类**。

### 任务清单

**PAR-1｜构建 paralog 对/组 ↔ 当前 DR-DL 的整合表（前置）**
- 把 ensembl paralog 边表 join 到当前 `18_gene_level_clustering` 的 DR/DL/`RevisedDeletion_essentiality`；保留 %identity 作为协变量。
- 产出：`results/paralog/HD_DIT_HAP/paralog_pairs_with_phenotype.tsv`（每对 paralog：两端 DR/DL/essentiality + %identity + paralog group size）。

**PAR-2｜冗余信号：paralog 对是否倾向「都非必需」**
- 对每对/每组：分类 (i) 都非必需 vs 一个/都必需（essentiality + DR 阈值）；(ii) paralog 间 DR-DL 距离 vs 随机背景（复用主题 D 的 permutation 方法）；(iii) 用 %identity 分层看冗余强度是否随序列相似度增加。
- 科学假设：高相似 paralog 对 → 功能冗余 → 双非必需比例显著高于随机。
- 产出：冗余分类汇总 + %identity vs 冗余的趋势图。

**PAR-3｜paralog 补偿解释验证失败/差异基因**
- 直接服务主题 A 的 VER-4：某基因单敲无表型，但其 paralog 高表达 → 代偿。用相对 mRNA 丰度（`mean_EMM_Proliferating_Cell_RNA_Abundance`）做代偿方向判断。
- 产出：一个「paralog 补偿」候选基因表，回填到验证 outcome 解读里。

**优先级**：PAR-1 前置 → PAR-2（系统化主体）→ PAR-3（与主题 A 联动）。**工作量**：中（数据分散在 3 个 repo，整合是主要成本）。

---

## 5. 主题 D — Coherence 分析（扩展 + 归因）

**数据**：`complex_coherence_metrics.csv`（173 行 complex × 49 列）。
- 术语源 = **仅 GO cellular_component macrocomplex**（PomBase `macromolecular_complex_annotation.tsv`，476 term）。不是 CORUM。
- 过滤：`DR(um)>0.3` 去掉非耗竭基因 → 每个 term 保留 `3≤n≤300` 成员 → 173 个 complex。
- 每个 dispersion 指标（mean/median distance-to-centroid、mean/median/max pairwise、mean kNN）都转成 **1000 次 permutation 的 z-score + p-value**（负 z = 比随机更紧 = coherent）。DIT-HAP 和 gRNA 各算一套。
- **最佳「incoherence」指标** = `DIT_HAP_median_pairwise_distance_zscore`（notebook 主轴；coherent 定义 = median & mean pairwise z 都 < −1）。
- 现状极值：最 coherent = small-subunit processome (z=−5.76)、线粒体大亚基核糖体 (−4.01)、CCT/chaperonin、MCM；最 incoherent = eIF3e subcomplex (z=+3.41)、43S 前起始复合体 (+2.54)、CLRC、eIF3、SAGA、Ino80、TORC1、Swr1。

### 任务清单

**COH-1｜扩展到更多 term（对应问题①：更多 term、符合/不符合预期的）**
- notebook 已经加载了带传播（`part_of`/`is_a`）的完整 BP/MF/CC GAF（`gene2go`/`go2genes`）——把 coherence 主循环从只喂 `macrocomplex` 改成也喂 GO BP / MF term group，即可扩展，无需新数据。
- KEGG pathway / CORUM 需要新注释文件（KEGG brite 有，CORUM 无）——列为可选。
- 谨慎下调 `n≥3` 下限（n<4 时 permutation z 不稳）。
- 产出：扩展后的 `coherence_metrics_all_namespaces.tsv`；标注「非常符合预期」（强 coherent，如已知稳定复合体）与「非常不符合预期」（本该 coherent 却 incoherent）的 term 清单。

**COH-2｜incoherence 的成因归因（对应问题②）**
现状**完全没有**归因逻辑，只算了 dispersion。要新建两类诊断：
- **(a) major/minor 亚基之分**：在 DR-DL 空间对某 incoherent complex 的成员做子聚类（2-component GMM / silhouette），看是否分成「核心紧簇 + 离群少数」。核心=化学计量主亚基，离群=调控/亚化学计量亚基。用特征矩阵的 mRNA/蛋白丰度、half-life、`evolutionary_rate` 佐证（旧 notebook cells 58–64 已原型化）。找公认好例子（如核糖体 core vs 松散结合因子）。
- **(b) 多复合体共享基因**：一个基因出现在多个 `GO_term_name` group 里 → 直接可数（现在没做）。共享基因会把两个复合体「拉」到中间，造成表观 incoherence。产出共享基因清单 + 具体例子（如 eIF3 vs eIF3e subcomplex 共享成员）。
- 产出：每个 incoherent complex 一行归因（major/minor 分裂 / 共享基因 / 数据缺失 / 真实生物学异质），带 2–3 个可写进论文的具体例子。

**优先级**：COH-1 中（改造现有循环）→ COH-2（科学价值最高，是「进一步分析」的核心）。**工作量**：COH-1 小–中；COH-2 中。

---

## 6. 主题 E — ML 分析（提升 R² + 可解释性）

**你的目标**：用 non-phenotypic 特征预测 DR/DL，当前 R²≈0.1，想再提高并做到可解释。

**现状**（`workflow/scripts/ml/{prepare_ml_data,train_automl}.py`）：
- 目标：DR、DL（`A` 带着但没建模）。
- **关键瓶颈**：`prepare_ml_data.py` 用 `dr_filter=0.3` 只保留 `DR>0.3` 的基因——4518 拟合基因里只剩 **1561**（→ ~1249 train / 312 test）。**近零 DR / WT 那一大团被整个丢掉，而不是建模**，这是 R² 天花板的主因，也截断了目标值域。
- 70 个特征列（正确排除了所有 phenotypic 标签，避免泄漏）；train-only Yeo-Johnson PowerTransformer 同时作用于 X 和 y；mljar 内部 impute。
- mljar `AutoML(mode=Explain|Perform, explain_level=2, total_time_limit=14400)`。旧 leaderboard：DR 保留集 R²≈0.09–0.10，DL≈0.16–0.30。确认「~0.1」。
- **可解释性工具已有**：mljar 写 per-model `learner_fold_*_shap_importance.csv` + permutation importance；`aggregate_feature_importance()` 求均值。旧 SHAP top 预测子：DR ← `mRNA_synthesis_rate_per_minute`（主导）、`evolutionary_rate`、氮饥饿 RNA 丰度、`GI_degree`、`mean.phylop`；DL ← `evolutionary_rate`（主导）。
- **注意**：`resources/curated/final_clusters.tsv` 在新 repo 里还没生成（只有 placeholder），ML pipeline 在新 repo 尚未跑过。所以本主题起步可先直接读旧 `18_gene_level_clustering` 的 DR/DL + 特征矩阵，不必等 finalize_gene_clusters。

### 任务清单（按性价比排序）

**ML-0｜建专用 conda env（前置）**
- 新建 env 装 `mljar-supervised`（版本 pin，见 [[phase2-core-chain-decisions]] 决策 3）+ `shap` + sklearn/pandas。所有其他 env 都没有这两个包。

**ML-1｜重新处理 `DR>0.3` 截断（最大杠杆，cheap win）**
- 方案 A：建全值域模型（不丢近零基因）。
- 方案 B：two-stage hurdle——先分类 essential/near-zero，再对正值域回归。
- 对比两者与现状的 CV R²。这一步预计对 R² 影响最大。

**ML-2｜报告 CV R² 而非单次 holdout（cheap win）**
- N≈1561 + 70 特征时，312 行 holdout 的 R² 方差很大。用 Perform 模式 5-fold，报 mean±SD。

**ML-3｜特征清理（cheap win，也提升可解释性）**
- 丢弃或加缺失指示：`copies_per_cell_*`（~49% 缺失）、`protein_half_life_minutes`（59% 缺失）——重 impute 可能加噪。
- 共线性剪枝：大量 `aa_percent_*`、`Mass/Residues/Peptide_length`、`ENC/GC3/tAIg/CAI` 冗余。做相关/VIF 过滤或分组，稳定正则线性模型。

**ML-4｜可解释基线：ElasticNet + SHAP（medium）**
- 在 Yeo-Johnson 特征上跑 ElasticNet，给出**有符号的加性系数**（效应方向），补 tree SHAP 只有重要性没方向的短板。

**ML-5｜GBM（CatBoost/LightGBM）+ 完整 SHAP 导出（medium）**
- 已知最好的学习器。导出 beeswarm / dependence / interaction values，显式呈现每个特征的效应方向与交互。

**ML-6｜目标工程 & 异方差（medium–high）**
- DL 零膨胀（中位数 0）→ Tweedie / log-hurdle / two-part 更合适；DR 右偏 → log/rank 目标。
- 残差扇形发散 → quantile GBM 给预测区间、处理非常数方差。

**ML-7｜分 regime 建模（high）**
- 复用 config 里已定义但闲置的 `dr_split_threshold=0.35` / `wt_cluster=9` 切分，对慢/快耗竭分别建模。

**优先级**：ML-0 前置 → ML-1/2/3（cheap wins，先做）→ ML-4/5（可解释性主体）→ ML-6/7（进阶）。

---

## 7. 主题 F — 非编码基因分析

**现状**：
- `non_coding_RNA_analysis.ipynb` 加载 `non_coding_rna.bed`（7484 ncRNA：lncRNA 7166 / tRNA 196 / snoRNA 66 / rRNA 49 / snRNA 7），merge `19_insertion_in_non_coding_genes/Non_coding_genes_Gene_level_statistics_fitted.tsv`（3135 行，3133 Success）。拟合参数含 `um`(DR)/`A`(DL)/`lam`。
- **目前只实际分析了 tRNA**（cells 6–12）：从 Systematic ID 解析氨基酸、从 GtRNAdb 名解析反密码子，`tRNA_copy_number` = 按 (氨基酸, 反密码子) 分组计数，`single_copy = copy_number==1`。加了 telomere/centromere Location、Marguerat-2012 mRNA 丰度。产出 `all_nuclear_tRNAs.xlsx`。cell 12 画单拷贝 tRNA 上下游 200bp 序列组成热图。
- **ncRNA 没有 essentiality「判定」**，只有连续 `um`；`FYPOviability` 基本 unknown，`DeletionLibrary_essentiality` 全 Not_determined。
- lncRNA/snoRNA/rRNA/snRNA **加载了但没处理**。拟合成功的：lncRNA 3022 / tRNA 48 / rRNA 35 / snoRNA 27 / snRNA 1。

### 任务清单

**NC-1｜tRNA 必需性的其他必要条件（对应问题①）**
- 现有 `all_nuclear_tRNAs.xlsx` 已有 copy_number、反密码子、Location、telo/centro 距离、mRNA 丰度。
- 新增维度检验：tRNA 表达量（需 tRNA-seq，缺）、**每个 tRNA 的 TTAA 位点数**（缺，需从基因组/`12_concatenated` 算，见 NC-5）、反密码子对应密码子的使用频率（codon usage，可从特征矩阵推）。
- 假设：single-copy 是必要非充分；叠加「高使用密码子 + 有可插入 TTAA」才更可能必需。
- 产出：tRNA 必需性多因素表 + single vs multi-copy 的 DR 比较（仅 48/196 有拟合，先做能做的）。

**NC-2｜lncRNA 还能发现什么（对应问题②）**
- 低垂：先把 3022 个已拟合 lncRNA 按 `um` 排序，挑高 `um` 的有表型候选（现成数据，无需新数据）。
- de-novo 新 lncRNA 发现需要原始插入坐标 vs 基因间区 + 表达证据，数据不全——列为后续。

**NC-3｜有表型的 ncRNA 的具体例子（对应问题③）**
- 在 `Non_coding_genes_Gene_level_statistics_fitted_annotated.tsv` 里按 Feature class 分别阈值 `um`（如 >0.5），挑 lncRNA/snoRNA/rRNA 里的高表型个例，逐个查文献注释。**现成数据，低垂果实**。

**NC-4｜基因上游区 与 丰度 的相关性**
- `upstream_and_downstream_analysis.ipynb` cells 33–36 **已做了 90%**：5UTR `um_ratio` vs log10 mRNA 丰度回归（slope/R²/p）。
- 扩展：延伸到 3UTR / N 端；分条件（增殖 vs 氮饥饿）；用蛋白丰度再做一遍。

**NC-5｜TTAA 分布 & tRNA 的 TTAA 情况**
- 需要 TTAA 位点坐标（从 `12_concatenated` 或基因组扫），与 171 个 tRNA 位点求交。
- 产出：每个 tRNA 的可插入 TTAA 数；解释「某 tRNA 无表型」是不是因为**根本没有可插入位点**（假必需/假非必需的技术性解释）。回填 NC-1。

**优先级**：NC-3 → NC-2（低垂果实，先出结果）→ NC-1 → NC-5（依赖 TTAA 计算）→ NC-4（扩展已有）。

---

## 8. 主题 G — 上游/下游 & 位置相关分析

**现状**（`upstream_and_downstream_analysis.ipynb`，39 cells，2.4MB）：
- `assign_UTR_type`（cells 6–7，`distance_threshold=400bp`）把基因间插入分类到 5UTR/3UTR，计算 **`um_ratio` = insertion_um / gene_um**（**这就是你说的「liwen 5UTR」概念** = UTR 插入 DR 对亲本基因 DR 归一化）和 `A_ratio`。
- cell 13 `plot_UTR_gene_boundary_insertions`：Upstream / N 端 / C 端 / Downstream 四联箱线图，按到边界距离分箱——**这直接对应你说的「中间表型弱、两边强」**（限 `um_gene>0.5` 的生长贡献基因）。forward/reverse 链分开（cells 18–21），有 `mug126`/`rib1` 等单基因图（cells 24–26）。
- cells 28–30：5UTR 聚类（`um_ratio` 在 −400→0bp 上插值 → clustermap）。
- **「piggyBac 下游 3x reporter」分析不存在**——字符串只出现在 base64 图片输出里，代码里没有。属于待做/规划。

### 任务清单

**UD-1｜「中间弱、两边强」的基因系统化（对应笔记「看基因中间表型弱，两边强」）**
- 复用 cell 13 的边界分箱逻辑，但改成**逐基因**打分：定义一个「N端+C端 DR vs 中段 DR」的比值/差值指标，扫全基因组，挑出中段插入耐受、两端敏感的基因。
- 生物学含义：中段可能是可容忍的 linker/loop，两端含关键结构域或调控元件。
- 产出：候选基因排序表 + 代表性 per-gene 插入-表型剖面图。

**UD-2｜piggyBac 下游 3x reporter 分析（新建，对应笔记）**
- 这是全新分析，代码不存在。**需要你先说明实验设计**：3x reporter 指什么读出、数据在哪、要回答什么问题（下游插入对报告基因表达的极性效应？）。列为**待澄清**任务。

**UD-3｜基因上游 与 丰度相关性（与 NC-4 合并）**
- 见 NC-4，同一套 `um_ratio` vs 丰度回归，扩展维度。

**UD-4｜低平台基因（对应笔记「看看关于低平台的能不能再拿出来」）**
- **待澄清**：「低平台」推测指 depletion 曲线 plateau 低（`A`/`DL` 小但 `DR` 不低？）的一类基因。需要你确认定义，再决定是从 `A`/`DL`/`y_inflection` 哪个参数切分、想回答什么。列为**待澄清**任务。

**UD-5｜找一两个「点相关」的例子（对应笔记「最好能找一两个点相关的」）**
- **待澄清**：「点相关」推测指单个插入位点层面的关联（如某特定 TTAA 位点的表型），而非 gene-level 汇总。需你确认，再从 `15_insertion_level_curve_fitting` 的位点级数据切入。列为**待澄清**任务。

**优先级**：UD-1 中（有现成逻辑可复用）；UD-2/UD-4/UD-5 待澄清后再排。

---

## 9. 需要你澄清的点（阻塞部分任务）

1. **piggyBac 下游 3x reporter（UD-2）**：实验读出是什么？数据文件在哪？要回答什么科学问题？
2. **「低平台」基因（UD-4）**：具体指哪个参数低（`A` / `DL` / plateau）？想从中「拿出来」什么？
3. **「点相关」（UD-5）**：是指 insertion-level 单点表型关联，还是别的？
4. **Notion 验证数据表**（你给的链接）：里面是否有本地 xlsx 之外的额外列（如已整理的 pathway/complex 归属）？如果有，VER-2 可以省去自己 join complex 的步骤——需要你导出或确认。
5. **优先做哪几个主题**：6 个主题全铺工作量很大，建议先确认 top-2/3。我的推荐见 §10。

---

## 10. 建议的执行顺序（我的推荐）

按「科学价值 × 数据就绪度 × 工作量」排：

1. **主题 A（验证）VER-1→VER-2/VER-5**——数据完全就绪，直接服务论文结论，且 VER-1 是很多分析的基础表。
2. **主题 D（coherence）COH-2**——「进一步分析」的核心诉求（成因归因），科学新意最高。
3. **主题 E（ML）ML-0→ML-1/2/3**——cheap wins 就能明显动 R²，且 ML-1（去掉 DR>0.3 截断）是明确的天花板。
4. **主题 C（paralog）+ 主题 F（ncRNA 低垂果实 NC-3/NC-2）** 并行——都能快速出个例结果。
5. 其余（COH-1 扩展、NC-1/5、UD-1）随后。
6. 待澄清任务（UD-2/4/5、Notion）等你回复。

每个主题落地为 `notebooks/<theme>/*.ipynb`（探索）+ 稳定部分下沉 `workflow/scripts/{stage}/`。所有新 stage 的 results 写到语义化目录 `results/{verification,paralog,coherence,ml,noncoding,utr}/HD_DIT_HAP/`。

---
---

# 深化设计：主题 A（验证深挖）+ 主题 D（coherence 归因）

> 以下 §A–§D 是 2026-07-19 头脑风暴确认后的可落地设计，优先于其他主题实施。

## §A 整体架构与数据流

两个主题共享一条数据地基，在 complex 层面交汇：

```
              ┌─ 上游只读数据 (DIT_HAP_pipeline 旧结果 + release/) ─┐
              │  验证表(411)  DR-DL拟合  GO-CC/BP注释  KEGG brite   │
              └──────────────────────┬────────────────────────────┘
                                     │
      [S0] 构建共享地基表 (notebooks/verification_complex/00_foundation.ipynb)
           · gene → {GO-CC complex, GO-BP term, KEGG pathway} 长表映射
           · gene → {DR, DL, 文库essentiality, 验证essentiality, outcome} 主表
                                     │
            ┌────────────────────────┴────────────────────────┐
            ▼                                                  ▼
   主题A: 验证深挖 (01_verification_deepdive)         主题D: coherence归因 (02_coherence_attribution)
   [A1] outcome标签(VER-1)                            [D1] 扩展term (GO-CC/BP/KEGG)
   [A2] 正确基因落哪些分组(VER-2)                      [D2] incoherence成因诊断
   [A3] 新颖基因(VER-3)                                     · major/minor子聚类(+curated先验)
   [A4] 失败归因(VER-4)                                     · 多complex共享基因
   [A5] 差异解释(VER-5)
            └────────────────────┬───────────────────────────┘
                                 ▼
              [X] A×D 闭环 (10_verification_x_coherence.ipynb)
                 · 58个含flip基因的complex
                 · 用「验证后essentiality」重算coherence
                 · 找「验证前incoherent → 纠正后coherent」的例子
```

- **[S0] 地基表是唯一前置**，A 和 D 都从它出发，保证同一套基因分组 + 同一套 DR-DL。
- 全程先在 `notebooks/verification_complex/` 探索；`[S0]` 与 `[D1]` 的确定性计算稳定后选性下沉。
- pombase 版本统一 `2026-06-01`（与特征矩阵一致）。

## §B 地基表 [S0] + 主题 A

**[S0] 两张地基表**（`00_foundation.ipynb`，产出到 `results/verification/HD_DIT_HAP/`）：
1. `gene_annotation_long.tsv`：gene ×（分组类型, 分组ID, 分组名）长表。三源——GO-CC macrocomplex（`resources/external/pombase/2026-06-01/ontologies_and_associations/macromolecular_complex_annotation.tsv`，477 term/1957 基因）；GO-BP（GAF 带 `part_of`/`is_a` 传播）；KEGG pathway（`DIT_HAP_pipeline/resources/KEGG/pombe_kegg_brite.xlsx` 的 gene→KO→pathway）。附 `n_groups_per_gene`（为 D2 共享基因检测预埋）。
2. `verification_master.tsv`：411 行 + DR/DL/R²/RMSE + 文库essentiality + 验证essentiality + 质控字段（`filtering_fraction`/`total_colonies`/colony-area day3–6）。

**主题 A 五分析**（`01_verification_deepdive.ipynb`）：
- **[A1] outcome 标签（前置）**：`outcome ∈ {concordant_E, concordant_V, flip_V→E(80), flip_E→V(17), ambiguous}`。`ambiguous` 触发规则 = `comments` 中文关键词（如"部分可见"）∪ 18 个 `noised_genes` ∪ 低 `total_colonies` ∪ 低 `R²`。产出 `verification_outcome_table.tsv`。
- **[A2] 正确基因落哪些分组**：concordant + flip 基因按三套分组聚合，输出 `分组 × outcome` 汇总，标注「整组被改判」的分组。
- **[A3] 新颖基因**：`gene_name` 仍为 systematic-ID 形式的未表征基因 ∩ 高 R² ∩ outcome 明确，优先在 80 个 flip_V→E 里找。
- **[A4] 失败归因**：ambiguous/矛盾基因按「数据质量（`filtering_fraction`/`total_colonies`/`R²`）vs 生物学（paralog 补偿 / FYPO 条件依赖）」分类。
- **[A5] 差异解释**：97 个 flip 交叉 `FYPOviability`（41 条件依赖）+ colony-area 曲线 + DIT-HAP/gRNA 双平台确认。

## §C 主题 D（coherence 扩展 + 归因）

产出到 `results/coherence/HD_DIT_HAP/`（`02_coherence_attribution.ipynb`）。

**[D1] 扩展 term**（改造现有循环，无需新数据）：
- 现状只跑 GO-CC 173 个 complex。把 coherence 主循环改成也喂 GO-BP term group + KEGG pathway（来自 [S0]）。
- 沿用量化：geometric-median centroid + 1000 次 permutation z-score，`median_pairwise_distance_zscore` 为主轴；`n≥3` 下限保留（n<4 时 z 不稳）。
- 产出 `coherence_metrics_all_namespaces.tsv`（加 `namespace ∈ {CC,BP,KEGG}`）+ 两张清单：「非常符合预期」（强 coherent）与「非常不符合预期」（本该紧却 incoherent）。→ 回答问题①。

**[D2] incoherence 成因诊断**（全新）：对每个 incoherent complex（z > 阈值）跑两条诊断线：
- **(a) major/minor 亚基分裂**：DR-DL 空间 2-component GMM + silhouette 判「核心紧簇 + 离群少数」；用特征矩阵 mRNA/蛋白丰度、half-life、`evolutionary_rate` 佐证核心=化学计量主亚基。**结合文献先验**：手工 curated `resources/curated/complex_subunit_roles.xlsx`（核心 vs 边缘/调控亚基名单，标注教科书/综述出处），GMM 结果与之对照。公认好例子：核糖体 core vs 松散因子、proteasome core vs regulatory。
- **(b) 多 complex 共享基因**：用 [S0] `n_groups_per_gene` 检测出现在多分组的基因——把复合体「拉」向别处造成表观 incoherence。产出共享基因清单 + 例子（如 eIF3 vs eIF3e subcomplex）。
- 每个 incoherent complex 输出归因标签 `{major/minor分裂 / 共享基因 / 数据缺失 / 真实异质}`。→ 回答问题②。

## §D A×D 闭环 + 固化与测试

**[X] 验证 × coherence 闭环**（`10_verification_x_coherence.ipynb`）：
- 聚焦 58 个含 flip 基因的 complex（交集已实测充足）。
- 每个：用「验证后 essentiality」替换文库判定，重算 coherence z-score，对比验证前后。
- 找 **「验证前 incoherent → 纠正后 coherent」** 的例子（抓手：90S preribosome / SNARE / small-subunit processome）。
- 反向：某 complex 因个别成员被文库误判而「显得矛盾」，验证澄清后认识如何改变。→ 回答问题①后半句。

**固化边界（先探索后固化）**：
- 确定性、稳定 → 下沉脚本：[S0] 三源 join（`workflow/scripts/verification/build_foundation.py`）、[D1] 扩展 coherence（`workflow/scripts/coherence/compute_coherence.py`，含 permutation，pin `random_state=42`）。各配 `.smk` 规则。
- 含人工判断 → 留 notebook + curated：ambiguous 判定、major/minor 角色表、闭环例子文字解读。
- **测试**：pytest 加入现有 `tests/`——地基表行数/唯一基因数断言；coherence 在已知 coherent complex（small-subunit processome，z 应显著为负）上的回归测试。

**执行前置**：
- 用 `data_analysis` conda env（有 openpyxl/pandas 2.2/sklearn）。
- 别经 worktree 的 snakemake 跑（Snakefile 硬编码 workdir 指向主树）；直接跑脚本。
- worktree 已软链 `resources/external`、`resources/literature`、`results` 到主树。

**交付物**：notebooks（`00_foundation`/`01_verification_deepdive`/`02_coherence_attribution`/`10_verification_x_coherence`）；results（`results/{verification,coherence}/HD_DIT_HAP/`）；curated（`complex_subunit_roles.xlsx`）；本计划文档。



