# 论文图版说明（7a_figures.smk）

**用途**：本文件是最终正文图 / 附图的唯一说明文档。每个 figure 一节，每个 panel 一小节，按"输入 → 变换 → 画法 → 版式数字 → 产物"五要素写清。

**唯一真相**：脚本常量是**执行真相**，本文件是**描述**。两者不一致时以脚本为准，并在改图时同步更新本文件（`workflow/scripts/figures/figureN.py` → 本文对应小节）。

**Stage 约定**：

- 规则：`workflow/rules/7a_figures.smk`，一个 figure 一条 rule（`plot_figure1`、`plot_figure_s1` …）。
- 脚本：`workflow/scripts/figures/`，每个 figure 一个独立脚本。
- 产物：`results/7a_figures/`（**真实目录**，不是软链；worktree 里写这里不会污染主仓库）。
- 环境：`workflow/envs/final_figures.yml`（cnsplots 0.8.0 + pypdf + poppler）。
- 脚本**不 import `workflow/src`、不调用任何上游 renderer**，只读上游已产出的表。改图不会触发其它 stage 的重算。

**重建单张图**（worktree 里必须带 `--directory`）：

```bash
env -i HOME="$HOME" USER="$USER" TMPDIR=/tmp \
  CONDA_PKGS_DIRS=/data/c/yangyusheng_optimized/DIT_HAP_analysis/.snakemake/conda-pkgs \
  PATH="/data/a/yangyusheng/miniforge3/condabin:/data/a/yangyusheng/miniforge3/bin:/usr/local/bin:/usr/bin:/bin" \
  /data/a/yangyusheng/.local/share/mamba/envs/snakemake/bin/snakemake \
  --directory "$(pwd)" --use-conda --cores 2 \
  --conda-prefix /data/c/yangyusheng_optimized/DIT_HAP_analysis/.snakemake/conda \
  -R plot_figure1 -- results/7a_figures/figure1.pdf
```

`-R` 是必须的：脚本不在 rule 的 `input` 里，改脚本不会触发重跑；`--` 之前不能跟目标（会被 `--forcerun` 吃掉）。其它图把 `plot_figure1` 与目标路径换成对应的（如 `plot_figure_s1` / `results/7a_figures/figure_s1.pdf`）。

---

## Figure 1 — DIT-HAP 流程、覆盖度与 spike-in 线性

**画布**：493 × 273.217 pt（Cell 双栏宽）。高度由 A 面板的比例推出：

```
10（A 字母带） + 123.217（= 493 × 113.4 / 453.72） + 6（行距）
+ 10（B/C 字母带） + 80（面板高） + 42（下方装饰） + 2（底边）
```

**产物**：

| 文件 | 说明 |
|---|---|
| `results/7a_figures/figure1.pdf` | 最终矢量图（单页） |
| `results/7a_figures/figure1.review.png` | 360 dpi 评审位图，由最终 PDF 栅格化（含 A 面板） |
| `results/7a_figures/_work/figure1_canvas.pdf` | 叠加前的画布，A 区为空白；排查版式时看这个 |

**为什么 A 面板用 pypdf 叠加**：matplotlib / cnsplots 没有把 SVG/PDF 当矢量嵌入的 API，`ax.imshow` 一定栅格化。脚本先用 matplotlib 画整张画布（A 区留白），再用 `pypdf` 把示意图 PDF 等比缩放、平移进 A 区，最后 `merge_transformed_page` 合成单页 PDF。示意图的 4 个 Arial 子集与 4 张 601 ppi 位图都原样保留。

---

### Panel A — DIT-HAP 工作流示意图（外部绘制）

| 要素 | 内容 |
|---|---|
| 输入 | `resources/figures/figure1a.pdf`（单页，页面框即内容框，453.72 × 113.4 pt） |
| 变换 | 等比缩放到 493 pt 宽（×1.0865），高 123.217 pt；在预留区内水平/垂直居中（本例正好铺满，无偏移） |
| 画法 | matplotlib 只建一个空 axes（`set_axis_off()`）；`pypdf` 的 `PdfWriter(clone_from=canvas)` + `merge_transformed_page(schematic, Transformation().scale(s).translate(x,y))`；坐标原点是 PDF 左下角 |
| 版式 | 预留框 `x=0, y_top=10, w=493, h=123.217`；不裁边、不补白边 |
| 产物 | 并入 `figure1.pdf` |

可编辑母版是 `resources/figures/figure1a.svg`；**管线只读 PDF**。SVG 换新后需重新导出 PDF（当前 PDF 由 PDFKit 生成，文字保持矢量）。

---

### Panel B — Traditional Tn-Seq vs DIT-HAP 的 LFC 覆盖度

| 要素 | 内容 |
|---|---|
| 输入 | `results/2a_coverage/cross_dataset/lfc_coverage_comparison.tsv`（`compare_lfc_coverage` 产出） |
| 必需列 | `dataset`、`fypo_class`、`n_lfc`、`pct_lfc`（缺列在 I/O 边界直接报错） |
| 变换 | 只取 `LD_haploid`（显示名 `Traditional Tn-Seq`）与 `HD_DIT_HAP`（显示名 `DIT-HAP`）；`coverage_pct = pct_lfc × 100` |
| 类顺序 | `inviable → viable → depends_on_conditions → unknown → total`，显示名 `Inviable / Viable / Depends on\nconditions / Unknown / Total` |
| 画法 | `cns.barplot(x="coverage_pct", y="class_label", orient="h", hue="dataset_label", hue_order=[Traditional Tn-Seq, DIT-HAP], palette=SERIES_COLORS[2,3], order=类顺序)`；柱端 `ax.bar_label` 标 `n_lfc` 千分位整数，无星号、无脚注 |
| 高亮 | Inviable 组（索引 0）画 `axhspan(-0.5, 0.5)`，accent 色（palette 0）8% alpha；该组刻度文字同色加粗 |
| 轴 | x `Covered genes (%)`，刻度 0/20/40/60/80/100，上限 112 只给柱端数字留头空；y `FYPO viability class` |
| 版式 | axes `x=130, y_top=149.217, w=120, h=80` pt；图例 `bbox_to_anchor=(0.5, -0.34)`、两列、无边框 |
| 产物 | 并入 `figure1.pdf` |

数据口径：`n_lfc` = 至少一个通过 hard filter + 重复合并（`release/insertion_level/LFC.tsv`）且满足 `IN_GENE_FILTER` 的插入所覆盖的基因数；分母为 `gene_result` 的 5,126 个蛋白编码基因。无置信区间、无 p 值（2026-11-07 用户决定）。

---

### Panel C — Spike-in 稀释线性

| 要素 | 内容 |
|---|---|
| 输入 | `results/1a_pcr_qc/spike_in_stats.tsv` |
| 必需列 | `Sample`、`Name`、`Relative_Dilution_Ratio`、`Relative_Read_Ratio` |
| 变换 | 丢弃 `Spikein0`（零稀释参考，读数为 0，会把拟合锚到 −inf）；`Name` 的 `Spike-in Insertion ` 前缀换成 `#`；按 `Name` 稳定排序固定颜色顺序 |
| 画法 | `cns.regplot(x="Relative_Dilution_Ratio", y="Relative_Read_Ratio", color="Name", add_equation=True)`——`color=` 而不是 `hue=`：每个位点一个图例项，但全局只拟合一条线。cnsplots 自带左上 `r / P` 标注与右下拟合方程 + R² |
| 图例 | `cns.take_legend_out(title="Spike-in")`，右侧竖排 #1–#5 |
| 轴 | `log₂(relative dilution ratio)` / `log₂(relative read ratio)` |
| 版式 | axes `x=300, y_top=149.217, w=110, h=80` pt |
| 产物 | 并入 `figure1.pdf` |

不用 y=x 参考线：稀释比是设定值、读数是测量值，参考线会暗示两者应当逐点相等。

---

### 面板字母

三个字母都是手摆（不用 `mp.panel` 自动标签）：8 pt 粗体 Arial，`fig.text(x, 1 − (y_top − 2)/画布高, va="bottom", ha="left")`，即每个面板左上角上方 2 pt。A 在 `(0, 10)`，B 在 `(130, 149.217)`，C 在 `(300, 149.217)`。

---

### 可调旋钮（全部在 `workflow/scripts/figures/figure1.py` 顶部）

| 常量 | 当前值 | 作用 |
|---|---:|---|
| `CANVAS_WIDTH_PT` | 493 | 画布宽（Cell 双栏） |
| `PANEL_LETTER_ROOM_PT` | 10 | 每行上方给面板字母的净空 |
| `ROW_GAP_PT` | 6 | A 与下排之间 |
| `ROW_PANEL_HEIGHT_PT` | 80 | B/C 的 axes 高度 |
| `ROW_BELOW_DECORATION_PT` | 42 | x 刻度 + xlabel + B 图例所占高度 |
| `PANEL_B_X_PT` / `PANEL_B_WIDTH_PT` | 130 / 120 | B 的横轴位置与轴宽 |
| `PANEL_C_X_PT` / `PANEL_C_WIDTH_PT` | 300 / 110 | C 的横轴位置与轴宽 |
| `LEGEND_BELOW_Y` | −0.34 | B 图例在 axes 坐标系中的纵向锚点 |
| `BAR_COUNT_FONT_PT` | 5 | 柱端数字字号（fit-critical 覆盖，house 允许） |

改高度类常量后画布高度自动重算；改 `PANEL_*_X_PT` 时注意 C 的图例在轴右侧约 55 pt，`PANEL_C_X_PT + PANEL_C_WIDTH_PT` 不宜超过约 420。

---

## Figure S1 — 文库制备重复性：PBL/PBR、技术重复、生物学重复

**画布**：493 × 162 pt（Cell 双栏宽）：

```
10（字母带） + 120（面板高） + 30（x 刻度 + xlabel） + 2（底边）
```

三个 120 pt 方形面板并排，间距 40 pt（`PANEL_PITCH_PT = 160`）。没有外部矢量面板，因此不需要 pypdf 叠加。

**产物**：

| 文件 | 说明 |
|---|---|
| `results/7a_figures/figure_s1.pdf` | 最终矢量图（单页） |
| `results/7a_figures/figure_s1.review.png` | 360 dpi 评审位图 |

**为什么复制而不是复用**：这三张图原本由 `workflow/scripts/pcr_qc/plot_pcr_qc.py` 画在 1a 的 2×2 图里。本脚本自带一份密度着色、对数轴与 n/r/P 标注的实现（KDE 部分对应 `workflow/src/figure_render/_point_density.py`，其余对应 `_scatter.py` 的密度路径），不 import 任何 stage 模块——改 S1 的版式不会动 1a 的 QC 图。

**输入为什么是 `_work/` 下的 parquet**：1a_pcr_qc 对 A–C 三个面板只发布了一张 2×2 PDF，没有发布表；`pbl_pbr / tech / bio.parquet` 就是它准备的表，rule 把它们声明为 input，1a 重算时 S1 会自动跟着重建。

---

### Panel A — PBL vs PBR（同一文库两个引物）

| 要素 | 内容 |
|---|---|
| 输入 | `results/1a_pcr_qc/_work/pbl_pbr.parquet` |
| 必需列 | `PBL`、`PBR`（缺列在 I/O 边界直接报错） |
| 变换 | 只保留两列都 > 0 的行（对数轴画不出 0）：164,820 → 108,717 行 |
| 画法 | 在 log10 空间做 2D 核密度（`gaussian_kde`，>5000 点时固定种子抽 5,000 点、128×128 网格插值），逐点按密度排序后 `ax.scatter(s=1.5, alpha=0.6, linewidths=0, rasterized=True)`，viridis 着色并归一化到 [0,1]；灰色虚线 x=y；左上 `n / r / P`（log10 空间 Pearson，`P < 1e-9` 显示为 `P < 1E-9`，否则 `P = mantissa E±exp`）；右下 inset 密度色条（0→low、1→high，标题 Density）；双对数轴（本图中 matplotlib 只给 y 轴次刻度，按 house 规则两侧一并关闭，与 1a 的 QC 图一致） |
| 版式 | axes `x=30, y_top=10, w=120, h=120` pt；字母 A 在 `(30, 10)` |
| 产物 | 并入 `figure_s1.pdf` |

### Panel B — 技术重复

| 要素 | 内容 |
|---|---|
| 输入 | `results/1a_pcr_qc/_work/tech.parquet`（同一份样品在两个上游项目各跑一次，按 (Chr, Coordinate, Strand) 合并） |
| 必需列 | `Reads_1`、`Reads_2` |
| 变换 | 两列都 > 0：152,276 → 117,297 行 |
| 画法 | 与 A 相同；轴标签 `Reads of Technical Replicate 1 / 2` |
| 版式 | axes `x=190`，其余同 A；字母 B 在 `(190, 10)` |
| 产物 | 并入 `figure_s1.pdf` |

### Panel C — 生物学重复

| 要素 | 内容 |
|---|---|
| 输入 | `results/1a_pcr_qc/_work/bio.parquet`（同一项目内两个样品） |
| 必需列 | `Reads_1`、`Reads_2` |
| 变换 | 两列都 > 0：125,566 → 92,704 行 |
| 画法 | 与 A 相同；轴标签 `Reads of Biological Replicate 1 / 2` |
| 版式 | axes `x=350`，其余同 A；字母 C 在 `(350, 10)` |
| 产物 | 并入 `figure_s1.pdf` |

---

### 可调旋钮（全部在 `workflow/scripts/figures/figure_s1.py` 顶部）

| 常量 | 当前值 | 作用 |
|---|---:|---|
| `CANVAS_WIDTH_PT` | 493 | 画布宽（Cell 双栏） |
| `PANEL_SIZE_PT` | 120 | 每个方形面板的边长 |
| `PANEL_PITCH_PT` | 160 | 面板起点间距（= 边长 + 40 pt 放刻度与 ylabel） |
| `PANEL_FIRST_X_PT` | 30 | 第一个面板的横轴位置（给 y 刻度留左边距） |
| `BELOW_DECORATION_PT` | 30 | x 刻度 + xlabel 所占高度 |
| `PANEL_LETTER_ROOM_PT` | 10 | 面板上方给字母的净空 |
| `DENSITY_SCATTER_S` / `DENSITY_SCATTER_ALPHA` | 1.5 / 0.6 | 点大小与透明度（密云例外） |
| `DENSITY_CBAR_BOUNDS` | (0.82, 0.02, 0.035, 0.30) | 密度色条在 axes 分数坐标中的位置与大小 |

改 `PANEL_SIZE_PT` 后画布高度自动重算；改 `PANEL_PITCH_PT` 时注意 `PANEL_FIRST_X_PT + 2 × PITCH + SIZE` 必须 < 493。

---

## Supplementary figures

（其余 S 图待加入：每个 S 图同样一条 rule、一节说明。）
