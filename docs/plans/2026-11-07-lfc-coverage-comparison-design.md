# LFC 基因覆盖度对比设计（2a_coverage 扩展）

**日期**：2026-11-07  
**状态**：设计定稿，待实现  
**范围**：`workflow/rules/2a_coverage.smk`、`workflow/src/coverage/`、`workflow/scripts/coverage/`、`tests/test_coverage.py`

## 1. 背景与目标

项目里有三套可比的插入文库：

- `HD_DIT_HAP`：高密度、二倍体构建、3 个生物学重复；
- `LD_DIT_HAP`：低密度、二倍体构建、3 个生物学重复；
- `LD_haploid`：低密度、单倍体、1 个样本。

目标：按 `FYPOviability` 分层，比较三套数据在初始时间点（t0）的基因覆盖度，突出 HD/LD 在 essential 基因上的优势。

结论先行：

| FYPO 类 | LD_haploid | HD_DIT_HAP | LD_DIT_HAP | 分母 |
|---|---:|---:|---:|---:|
| inviable（essential） | 767 (63.8%) | 1082 (89.9%) | 1044 (86.8%) | 1203 |
| viable | 3363 (93.9%) | 3159 (88.2%) | 3050 (85.1%) | 3582 |
| depends_on_conditions | 142 (81.1%) | 159 (90.9%) | 154 (88.0%) | 175 |
| unknown | 118 (71.1%) | 111 (66.9%) | 92 (55.4%) | 166 |
| total | 4390 (85.6%) | 4511 (88.0%) | 4340 (84.7%) | 5126 |

因此主张限定为 **essential 覆盖优势**，并明确 viable 的取舍：

- essential：HD +26.1 个百分点，LD +23.0 个百分点；
- viable：HD −5.7 个百分点，LD −8.8 个百分点；
- total：HD +2.4 个百分点，LD −0.9 个百分点。

不能写“DIT-HAP 覆盖的基因总数更多”——LD 的总数低于 haploid。

## 2. 定义与口径

### 2.1 覆盖定义

一个基因被算作“覆盖”，当且仅当它的注释区间内有至少一个插入位点出现在该数据集 release 的 `insertion_level/LFC.tsv` 中，且该插入通过 `IN_GENE_FILTER`：

```
Type != 'Intergenic region' and Distance_to_stop_codon > 4
```

LFC 表是上游 hard filter 和生物学重复合并后的产物：

- hard filter：每个样本在 `initial_time_point` 的 reads ≥ 8（`hard_filtering_cutoff: 8`）；
- 重复合并：HD/LD 走 DESeq2 分支，插入必须在该数据集所有样本中保留；`LD_haploid` 单样本，无合并。

因此 LFC 覆盖 = “有可用于下游分析的插入证据”的基因，而不是“注释表里出现过任意插入”的基因。

### 2.2 为什么用 LFC 而不是 DR

- DR 覆盖 = 基因水平曲线拟合成功（`gene_level/fitting_results.tsv` 的 `DR` 非空）。
- LFC 覆盖 = 插入层面有 LFC 记录。
- 两者在本项目只差 1–2 个基因（见 §3），因为拟合成功率接近 100%。
- LFC 更上游：不依赖曲线拟合；DR 覆盖把“文库覆盖”和“拟合成功”混在一起。
- 选择 LFC 作为主口径；DR 作为并列表列，供对照。

### 2.3 FYPO 分类与分母

`FYPOviability` 四个取值全列：

- `viable`
- `inviable`
- `depends_on_conditions`
- `unknown`

另加 `total` 行。分母来自 1c 注释参考的蛋白编码基因全集：**5,126 个有 FYPO 注释的基因**。17 个无 FYPO 注释的基因不入分母，在日志与附表单列。

`inviable` 是本文的 essential 类；`depends_on_conditions` 单独成类，不并入 essential。

### 2.4 t0 口径

- `LD_haploid`：`initial_time_point = 0h`，`0h` 与 `YES0` 未合并；
- `HD_DIT_HAP` / `LD_DIT_HAP`：`0h` 与 `YES0` 合并为 `YES0`，release 已 drop `0h`。

LFC 表使用各项目自己的 `initial_time_point` 做 hard filter；覆盖比较不另取时间点列。

### 2.5 对照集

主比较三套：`HD_DIT_HAP`、`LD_DIT_HAP`、`LD_haploid`。

不纳入：

- `HD_diploid`：registry 有，但 `release/` 不在磁盘上，无法构建；
- `Spore2YES6_1328`：流程对照，实验口径不同；
- `Spikein`：无基因水平数据。

### 2.6 统计

不写 p 值，也不给置信区间。LFC 表已经是 hard filter + 生物学重复合并后的产物：
HD/LD 要求插入在所有样本中保留，`LD_haploid` 只有单样本。按基因数计的 Wilson
区间只反映基因之间的二项波动，不代表重复或建库不确定性；三个库又对同一批基因
做覆盖，观测也不独立。覆盖率按点估计陈述，差异用百分点表示。

## 3. 预期数字

### 3.1 LFC 覆盖与 DR 覆盖

| FYPO 类 | LD_haploid LFC / DR | HD_DIT_HAP LFC / DR | LD_DIT_HAP LFC / DR | 分母 |
|---|---:|---:|---:|---:|
| inviable | 767 / 767 | 1082 / 1082 | 1044 / 1044 | 1203 |
| viable | 3363 / 3364 | 3159 / 3161 | 3050 / 3052 | 3582 |
| depends_on_conditions | 142 / 142 | 159 / 159 | 154 / 154 | 175 |
| unknown | 118 / 118 | 111 / 111 | 92 / 92 | 166 |
| total | 4390 / 4391 | 4511 / 4513 | 4340 / 4342 | 5126 |

两个定义只差 1–2 个基因，全部落在 viable 类。

### 3.2 覆盖率

| FYPO 类 | LD_haploid | HD_DIT_HAP | LD_DIT_HAP |
|---|---:|---:|---:|
| inviable | 63.8% | 89.9% | 86.8% |
| viable | 93.9% | 88.2% | 85.1% |
| depends_on_conditions | 81.1% | 90.9% | 88.0% |
| unknown | 71.1% | 66.9% | 55.4% |
| total | 85.6% | 88.0% | 84.7% |

数字由当前 release 直接计算；实现后以脚本输出为准，差 1–2 个基因属预期。

## 4. 架构与 DAG

不改 stage 编号，在 `2a_coverage.smk` 追加两条 rule：

```
prepare_coverage_data (per dataset, 已有)
    -> compute_lfc_coverage (per dataset)
        inputs:  release/insertion_level/LFC.tsv
                 release/insertion_level/annotations.tsv.gz
                 results/2a_coverage/{dataset}/_work/gene_result.parquet
        outputs: results/2a_coverage/{dataset}/lfc_coverage_stats.tsv
                 results/2a_coverage/{dataset}/_work/lfc_gene_coverage.parquet
    -> compare_lfc_coverage (cross-dataset)
        inputs:  三个数据集的 lfc_coverage_stats.tsv
        outputs: results/2a_coverage/cross_dataset/lfc_coverage_comparison.tsv
                 results/2a_coverage/cross_dataset/lfc_coverage_by_fypo.pdf
                 results/2a_coverage/cross_dataset/lfc_coverage_by_fypo_horizontal.pdf
                 results/2a_coverage/cross_dataset/lfc_coverage_by_fypo_two_libraries.pdf
                 （各 + .review.png）
```

### 4.1 为什么是两条 rule

`compute_lfc_coverage` 复用 `_work/gene_result.parquet`，拿到：

- 5,126 个蛋白基因全集；
- `FYPOviability` 注释；
- 每个基因的 `DR`（用于 `n_dr` 对照列）。

不需要再读 1c 注释参考。每条 rule 的输出只被下一条消费，没有像 `annotations.parquet` 那样的多点 fan-out，所以不需要再拆 prepare / compute。

### 4.2 数据集列表

`config/analysis.yaml` 新增：

```yaml
# --- Gene insertion coverage (2a_coverage.smk) ---
coverage:
  # 跨数据集 LFC 覆盖对比的成员；顺序决定图的柱序。
  lfc_comparison_datasets:
    - LD_haploid
    - HD_DIT_HAP
    - LD_DIT_HAP
```

`compare_lfc_coverage` 按这个列表展开输入；缺一个数据集时 Snakemake 在路径解析阶段失败。

### 4.3 `rule all`

先不加，按输出路径单独构建，与 6a 的现状一致：

```bash
snakemake --use-conda --cores 8 \
  results/2a_coverage/cross_dataset/lfc_coverage_comparison.tsv
```

确认稳定后再决定是否加入 `rule all`。

## 5. 输出契约

### 5.1 每数据集 stats TSV

`results/2a_coverage/{dataset}/lfc_coverage_stats.tsv`

| 列 | 类型 | 说明 |
|---|---|---|
| `dataset` | str | 数据集名 |
| `fypo_class` | str | `inviable` / `viable` / `depends_on_conditions` / `unknown` / `total` |
| `n_total` | int | 该类基因数（分母） |
| `n_dr` | int | DR 覆盖数（参考列） |
| `n_lfc` | int | LFC 覆盖数（主口径） |
| `pct_dr` | float | `n_dr / n_total`，0–1 |
| `pct_lfc` | float | `n_lfc / n_total`，0–1 |

行序固定：`inviable`、`viable`、`depends_on_conditions`、`unknown`、`total`。
`pct_*` 以 0–1 存储；图的百分号在渲染层乘 100。

### 5.2 跨数据集 comparison TSV

`results/2a_coverage/cross_dataset/lfc_coverage_comparison.tsv`

列与 §5.1 相同，三个数据集纵向拼接；`dataset` 列区分。行序按 `coverage.lfc_comparison_datasets`，组内按固定 FYPO 顺序。

### 5.3 图

三张图共用同一份 comparison TSV。

**① `lfc_coverage_by_fypo.pdf`（+ `.review.png`，竖向）**

- 单面板，5 组：`inviable`、`viable`、`depends_on_conditions`、`unknown`、`total`；
- x = FYPO 类；每组 3 根柱（`LD_haploid`、`HD_DIT_HAP`、`LD_DIT_HAP`）；
- y = 覆盖率 0–100%，不截断；
- 柱顶标签 = `n_lfc`（绝对数）；
- `total` 组的 `LD_DIT_HAP` 柱加星号，脚注：

  ```
  LD_DIT_HAP 4340 < LD_haploid 4390
  ```

- 标题中立：`Gene coverage by FYPO viability class`；
- 颜色沿用 cnsplots house style；`inviable` 组用强调色；
- 用 `save_dual()` 输出 PDF + PNG，不硬编码尺寸。

**② `lfc_coverage_by_fypo_horizontal.pdf`（+ `.review.png`，x/y 互换）**

- 与 ① 同数据、同柱数，改为横向：y = FYPO 类，x = 覆盖率 0–100%；
- 柱顶标签、`inviable` 强调、`total` 星号与脚注均保留。

**③ `lfc_coverage_by_fypo_two_libraries.pdf`（+ `.review.png`，两库横置）**

- 同样横置，只画 `LD_haploid` 与 `HD_DIT_HAP`；
- 图例显示名：`LD_haploid` → `Traditional Tn-Seq`，`HD_DIT_HAP` → `DIT_HAP`；
- 两库版本不含 `LD_DIT_HAP`，因此不加星号/脚注。

## 6. 代码组织

### 6.1 `workflow/src/coverage/core.py`

新增纯函数：

- `build_lfc_gene_table(lfc_path, annotations_path, gene_result: pd.DataFrame) -> pd.DataFrame`
  - 复用 `load_insertion_level(lfc_path, annotations_path)`：LFC 表作为第一个参数，函数返回的 annotations 已按 LFC index 重排并去重；
  - 对 annotations 应用 `IN_GENE_FILTER`，取 `Systematic ID` 唯一集合；
  - 以 `gene_result` 的 `Systematic ID` 为全集，左连接得到 `lfc_covered`；
  - `dr_covered = DR.notna()`（`gene_result` 已带 DR）；
  - 保留 `Systematic ID`、`FYPOviability`、`lfc_covered`、`dr_covered`。
- `summarize_lfc_coverage(gene_table: pd.DataFrame) -> pd.DataFrame`
  - 按 `FYPOviability` 分组，加 `total` 行；
  - 输出 §5.1 的列（不含 `dataset`，由脚本补）；
  - `n_total` 只计非空 `FYPOviability`；17 个 NaN 由脚本日志单列。

### 6.2 `workflow/scripts/coverage/compute_lfc_coverage.py`

CLI：

```
--lfc                  {release}/insertion_level/LFC.tsv
--annotations          {release}/insertion_level/annotations.tsv.gz
--gene-result          results/2a_coverage/{dataset}/_work/gene_result.parquet
--dataset              {dataset}
--output-stats         results/2a_coverage/{dataset}/lfc_coverage_stats.tsv
--output-gene-table    results/2a_coverage/{dataset}/_work/lfc_gene_coverage.parquet
```

遵循 `python-script-conventions`：argparse + loguru，`main() -> int`，边界 fail-fast，core 函数不写防御分支。

### 6.3 `workflow/scripts/coverage/compare_lfc_coverage.py`

CLI：

```
--stats                        三个 lfc_coverage_stats.tsv（nargs="+"）
--output-tsv                   results/2a_coverage/cross_dataset/lfc_coverage_comparison.tsv
--output-figure                results/2a_coverage/cross_dataset/lfc_coverage_by_fypo.pdf
--output-figure-horizontal     results/2a_coverage/cross_dataset/lfc_coverage_by_fypo_horizontal.pdf
--output-figure-two-libraries  results/2a_coverage/cross_dataset/lfc_coverage_by_fypo_two_libraries.pdf
```

读入三份 stats，纵向拼接，写 TSV，再渲染 §5.3 的三张图。图读的就是写出的 TSV 内容，
保证图/表一致；后两个参数可选，省略时只画竖向图。

### 6.4 `workflow/rules/2a_coverage.smk`

两条 rule：

- `compute_lfc_coverage`
  - `input`: `lfc`（release LFC）、`annotations`（release 注释）、`gene_result`（`_COVWORK/gene_result.parquet`）
  - `output`: `stats`、`gene_table`
  - `conda`: `../envs/statistics_and_figure_plotting.yml`
  - `log`: `logs/2a_coverage/compute_lfc_coverage_{dataset}.log`
- `compare_lfc_coverage`
  - `input`: 按 `config["coverage"]["lfc_comparison_datasets"]` 展开的 stats 列表
  - `output`: `comparison_tsv`、`figure_pdf`/`figure_review`、`figure_horizontal_pdf`/`figure_horizontal_review`、`figure_two_libraries_pdf`/`figure_two_libraries_review`
  - `conda`: `../envs/cnsplots.yml`
  - `log`: `logs/2a_coverage/compare_lfc_coverage.log`

### 6.5 `config/analysis.yaml`

新增 `coverage.lfc_comparison_datasets`（见 §4.2）。

### 6.6 测试

`tests/test_coverage.py` 新增：

- `build_lfc_gene_table`：小型合成表，覆盖 in-gene / intergenic、重复注释行、`lfc_covered` / `dr_covered` 的交叉；
- `summarize_lfc_coverage`：五行的列与顺序；`total.n_total == 各非 total 行 n_total 之和`；
- 跨数据集拼接：三份合成 stats → 15 行，列齐全。

测试不跑 Snakemake，只测 `workflow/src/coverage/`。

## 7. 验证

1. 构建三个每库 stats：

   ```bash
   snakemake --use-conda --cores 8 \
     results/2a_coverage/LD_haploid/lfc_coverage_stats.tsv \
     results/2a_coverage/HD_DIT_HAP/lfc_coverage_stats.tsv \
     results/2a_coverage/LD_DIT_HAP/lfc_coverage_stats.tsv
   ```

2. 构建跨库输出：

   ```bash
   snakemake --use-conda --cores 8 \
     results/2a_coverage/cross_dataset/lfc_coverage_comparison.tsv
   ```

3. 对照 §3：

   - `LD_haploid` total 4,390；`HD_DIT_HAP` 4,511；`LD_DIT_HAP` 4,340；
   - inviable 767 / 1,082 / 1,044；
   - viable 3,363 / 3,159 / 3,050。

4. 检查 DR 与 LFC 的差：−1 / −2 / −2。

5. 检查图：竖向与横置版本各 5 组 × 3 柱 = 15 根柱，`total` 组星号与脚注；两库版本 5 组 × 2 柱，图例为 `Traditional Tn-Seq` / `DIT_HAP`。

6. 质量门：

   ```bash
   ruff check workflow/scripts/coverage/compute_lfc_coverage.py \
              workflow/scripts/coverage/compare_lfc_coverage.py \
              workflow/src/coverage/core.py
   ruff format --check <同上>
   ty check <同上>
   pytest tests/test_coverage.py
   ```

## 8. 论文写法与限制

可以写：

- “DIT-HAP 的 essential 基因 LFC 覆盖为 89.9% / 86.8%，haploid 文库为 63.8%。”
- “DIT-HAP 以 viable 覆盖的小幅下降换取 essential 覆盖的大幅提升。”

不能写：

- “DIT-HAP 覆盖的基因总数更多”（LD 4,340 < 4,390）；
- 任何 p 值或“显著”字样（未做独立性检验）。

限制：

- HD/LD 的 LFC 表经过 DESeq2 重复合并，要求插入在所有样本中保留；`LD_haploid` 单样本。两边的重复门槛不同。
- essential 优势在该更严格口径下仍成立（+26 / +23 个百分点），因此结论稳健。
- viable 的差距可能部分来自门槛差异；正文按观测值陈述，不做过强机制解释。
- 覆盖率是点估计，不给 CI 或 p 值（理由见 §2.6）。

## 9. 不做的事

- 不建 6b 阶段；
- 不替换 2a 现有的 DR 覆盖输出；
- 不做交互页；
- 不补 `HD_diploid`；
- 不改 upstream release；
- 不改 2a 现有 `coverage_stats.tsv`、`detailed_genes.xlsx` 与既有图。

## 10. 变更清单

| 文件 | 变更 |
|---|---|
| `config/analysis.yaml` | 新增 `coverage.lfc_comparison_datasets` |
| `workflow/src/coverage/core.py` | 新增 `build_lfc_gene_table`、`summarize_lfc_coverage`、`combine_lfc_coverage_stats`；`load_insertion_level` 增加注释表必需列校验 |
| `workflow/scripts/coverage/compute_lfc_coverage.py` | 新建 |
| `workflow/scripts/coverage/compare_lfc_coverage.py` | 新建 |
| `workflow/rules/2a_coverage.smk` | 新增 `compute_lfc_coverage`、`compare_lfc_coverage` |
| `tests/test_coverage.py` | 新增 LFC 覆盖与跨数据集拼接测试 |
| `docs/plans/2026-11-07-lfc-coverage-comparison-design.md` | 本文件 |
