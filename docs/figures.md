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

## Figure 1 — DIT-HAP 流程、覆盖度、spike-in、LFC 轨迹与 FACS 验证

**画布**：493 × 462.42 pt（Cell 双栏宽）。整页由 `cns.multipanel(max_width=493, unit="pt")` 流式排版，两个 `mp.newline()` 隔开三行；高度由各行面板尺寸推出（实测值）：

```
A 行：     10（字母带 = 7.2 标签高 + 2.8 净空） + 121.718（= 487 × 113.4 / 453.72） + 6（行距）      = 137.72
B/C/D 行： 7.2（字母带） + 17.3（图例，实测） + 100（面板高） + 26（下方装饰）                        = 150.50
E 行：     6（行距） + 7.2（字母带） + 161（容器高） + 0（末行不留下边距）                            = 174.20
```

B/C/D 用同一个 100 × 100 pt 方形 axes，三块绘图区大小一致、长宽相等（house 面板契约 80–150 pt 内）。E 是一个 487 pt 宽的空白容器 axes，内部再排 4 列 × 2 行的 mini axes（每个 70 × 55 pt，字号沿用 house 的 8/7 pt），右侧空出的竖带放图例。

**版式机制**：`mp.panel()` 的 `width/height` 是 axes 尺寸，面板的外框在渲染时再加上**实测**的 y 轴装饰宽度（刻度标签 + ylabel）与右边距 10 pt。所以 B/C/D 从左往右流式排列、行首贴页面左边；ylabel 变长只会把它自己这一列往右推，不会压到前一个面板。面板字母由 multipanel 放在每列外框的左缘（8 pt 粗体），A、B、E 的字母因此都在 x = 0。E 的 mini axes 用 `container.inset_axes()` 建，位置是容器 axes 的分数坐标，渲染时自动跟随容器。

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
| 变换 | 等比缩放（×1.0733，几乎满框）到预留的 axes 框；框的宽高比按示意图页面比例设定，因此正好铺满，无偏移 |
| 画法 | multipanel 建一个空 axes（`set_xticks([]) + set_yticks([]) + set_axis_off()`，清刻度是为了让 multipanel 量到 0 宽的左装饰）；`pypdf` 的 `PdfWriter(clone_from=canvas)` + `merge_transformed_page(schematic, Transformation().scale(s).translate(x,y))`；坐标原点是 PDF 左下角 |
| 版式 | 预留框 `x=5.775, y_top=10, w=487, h=121.718`（487 = 493 − 6 pt，6 pt 是字母 "A" 占的那一列）；不裁边、不补白边 |
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
| 版式 | axes `100 × 100` pt，由 flow 排到行首（x 含实测 ylabel 宽度）；图例两列、无边框，用 `bbox_to_anchor=(0, 1)` 挂在 axes 上方，由 multipanel 实测为顶部装饰 |
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
| 版式 | axes `100 × 100` pt，紧随 B；图例 `cns.take_legend_out` 移到 axes 上方右侧竖排 |
| 产物 | 并入 `figure1.pdf` |

不用 y=x 参考线：稀释比是设定值、读数是测量值，参考线会暗示两者应当逐点相等。

---

### Panel D — 基因级 LFC 时间轨迹

| 要素 | 内容 |
|---|---|
| 输入 | `<HD_DIT_HAP release>/gene_level/LFC.tsv` |
| 必需列 | `Systematic ID`、`FYPOviability`、`YES0`…`YES4`（缺列在 I/O 边界直接报错） |
| 变换 | 只留 `inviable` 与 `viable` 两类；宽表 melt 成每基因每世代一行；`YES0…YES4` 映射到世代 `0 / 2.352 / 5.588 / 9.104 / 12.48`（上游 config 的 time_points） |
| 画法 | 每个类别把全部基因 pivot 成"一列一基因"，一次 `ax.plot` 画完整类：线宽 0.5 pt、alpha 0.05，密集处靠叠加变深。inviable 用 accent（palette 0），viable 用 palette 1；inviable 的 zorder 更高，重叠处保留高亮色。`axhline(0)` 灰色虚线 |
| 图例 | 每类一个不透明色块（线本身不带 label，否则上千条曲线各成一个图例项），标签带样本量 `Inviable (n = 1,082)` / `Viable (n = 3,161)`；单列、无边框，横跨 axes 上边缘（下半在图内），`pad_top` 把下半高度还给行基线 |
| 轴 | x 刻度 `0 / 2.4 / 5.6 / 9.1 / 12.5`（世代数取三位小数会在 100 pt 内打架，标签只做显示舍入）；x `Generations`，y `Gene-level LFC (log₂)` |
| 版式 | axes `100 × 100` pt，紧随 C；x 刻度 + xlabel 由 `ROW_BELOW_DECORATION_PT` 预留 |
| 产物 | 并入 `figure1.pdf` |

---

### Panel E — FACS 验证基因的耗竭曲线 vs DIT-HAP

| 要素 | 内容 |
|---|---|
| 输入 | `resources/curated/FACS_depletion_curves_from_XH.xlsx`（XH 提供的 FACS 耗竭曲线）+ 上一个 panel 同源的 `gene_level/LFC.tsv` |
| 必需列 | FACS：`Name`、`Strain_alias`、`Generations`、`LFC`；release：`Name` 与 `YES0`…`YES4`（缺列在 I/O 边界直接报错） |
| 数据形状 | workbook 是 tidy 长表（每行 = 基因 × 菌株 × 时间点），10 个基因、每株 9 个时间点（Spore + YES0…YES7），共 171 行；panel 去掉 `orc6`、`srb4` 后画 8 个基因（`FACS_EXCLUDED_GENES`）；`fal1` 只有 `Strain#1`，没有 `Strain#2` |
| 变换 | 每个基因把 FACS 两条独立菌株曲线与 release 的 5 点 LFC 拼成一张长表：`gene / source / series / generation / lfc`。**符号对齐**：FACS 的 LFC 以耗竭为正，release 以耗竭为负，故 FACS 乘 −1 后与 panel D 同向（耗竭基因向下）。共享 y 范围取两源合并后的最小值/最大值向外取整 |
| 画法 | 一个 487 pt 容器 axes（`set_axis_off()`）内嵌 8 个 `inset_axes`，4 列 × 2 行、按基因名排序（每行 4 个）。每格画 3 条线：DIT-HAP（palette 4 绿）、FACS `Strain#1`（palette 5 鲑红）、FACS `Strain#2`（palette 6 紫）——两条是独立菌株而非技术重复，所以都画，各用一色、不用线型区分；三条数据线线宽 2 pt（用户指定，覆盖 house 默认 1.5）、`alpha=0.8`，图例色块保持不透明；零线灰色虚线 0.5 pt；基因名斜体加粗作 title（house 8 pt） |
| 轴 | 全格共享 x `0–12.8`（刻度 0/4/8/12）与 y（刻度间隔 2、向外取整）；只有第一列显示 y 刻度标签、只有底行显示 x 刻度标签；两个共享轴标题用容器 `text()` 居中在整个网格上，不用某个 mini axes 的 label |
| 图例 | 竖排在网格右侧的留白带里（左边缘距网格 12 pt），在网格高度上垂直居中：`DIT-HAP`、`FACS Strain#1`、`FACS Strain#2`，单列、无边框 |
| 版式 | mini axes `70 × 55` pt；列距 10 pt、行距 16 pt（给下一行 8 pt title 留位）；左侧 34 pt（y 刻度 + 旋转 y 标题）、底部 22 pt（x 刻度 + x 标题）、顶部 13 pt（首行 title）；右侧 12 pt 后是图例带；容器高 = 13 + 2×55 + 16 + 22 = 161 pt；E 是最后一行，`margin_bottom=0`（否则 multipanel 默认再留 10 pt 页脚） |
| 产物 | 并入 `figure1.pdf` |

数据口径：FACS 表的 `Generations` 每株各自记录（同名的 YES 时间点不等于同一世代值），因此不能共用 x 向量；`Spore` 恒为 0，有三个 block 的 YES0 略小于 0（如 `arg6/Strain#1` 的 −0.18），x 下限留到 −0.3。

---

### 面板字母

五个字母由 `mp.panel(label=...)` 自动放置：8 pt 粗体 Arial，每列外框左缘、面板顶端上方（A 用 `pad_top` 留 2.8 pt，D 用 `pad_top` 抵消图例下半，其余为 0）。`multipanel` 把它放在"实测左装饰外侧"，而列的外框左缘就是页面的列起点，所以 A、B、E 的字母都在 x = 0，不随 ylabel 宽度漂移。**不要**再手摆 `fig.text` 字母——会和 multipanel 的自动标签重影。

---

### 可调旋钮（全部在 `workflow/scripts/figures/figure1.py` 顶部）

| 常量 | 当前值 | 作用 |
|---|---:|---|
| `CANVAS_WIDTH_PT` | 493 | 画布宽（Cell 双栏），即 `multipanel(max_width=...)` |
| `PANEL_LABEL_COLUMN_PT` | 6 | 给每列字母留的那一列；`PAGE_CONTENT_WIDTH_PT = 493 − 6`（A 与 E 用它做宽度） |
| `PANEL_LETTER_PAD_TOP_PX` | 14 | A 上方字母净空（显示像素，≈2.8 pt @360 dpi） |
| `ROW_GAP_PT` | 6 | A 行与 B/C/D 行之间（A 的 `margin_bottom`）；E 的行距（`margin_top`） |
| `ROW_PANEL_SIZE_PT` | 100 | B/C/D 共用的方形 axes 边长（x 由 flow 决定） |
| `ROW_BELOW_DECORATION_PT` | 26 | x 刻度 + xlabel 所占高度（`margin_bottom`） |
| `BAR_COUNT_FONT_PT` | 5 | 柱端数字字号（fit-critical 覆盖，house 允许） |
| `LFC_LEGEND_HALF_HEIGHT_PT` | 8.4 | D 的图例跨过 axes 顶边时，`pad_top` 补偿的下半高度 |
| `FACS_GRID_COLUMNS` | 4 | E 的列数；行数由 `ceil(基因数 / 4)` 自动推出 |
| `FACS_MINI_AXES_WIDTH_PT` / `FACS_MINI_AXES_HEIGHT_PT` | 70 / 55 | E 每个基因 mini axes 的尺寸（70×4 + 10×3 + 34 = 344，右侧留给图例） |
| `FACS_EXCLUDED_GENES` | `orc6`、`srb4` | 不画的基因（用户 2026-11-12 决定），在 frame 构建时剔除 |
| `FACS_LINE_WIDTH_PT` | 2 | E 的三条数据线线宽（用户指定） |
| `FACS_COLUMN_GAP_PT` / `FACS_ROW_GAP_PT` | 10 / 16 | E 的列距/行距；行距要给下一行基因名 title 留位 |
| `FACS_TITLE_BAND_PT` | 13 | E 首行基因名 title 的顶部留白 |
| `FACS_LEGEND_GAP_PT` | 12 | E 网格右缘到右侧图例左缘的距离 |
| `FACS_BOTTOM_DECORATION_PT` | 22 | E 的 x 刻度标签 + 共享 x 标题 |
| `FACS_LEFT_DECORATION_PT` | 34 | E 的 y 刻度标签 + 共享 y 标题 |
| `FACS_LFC_SIGN` | −1 | FACS LFC 符号翻转，与 release 的"耗竭为负"对齐 |

改 axes 尺寸后画布高度自动重算；B/C/D 不会再被"摆到中间"，它们总是从行首左对齐。若想改变面板之间的空隙，调 `mp.panel(..., margin_right=...)`（默认 10 pt）。

---

## Figure S1 — 文库制备重复性 + 每基因插入条数分布

**画布**：493 × 304.4 pt（Cell 双栏宽），由 `cns.multipanel(max_width=493, unit="pt")` 流式排版：

```
A–C 行：7.2（字母带） + 120（面板高） + 30（x 刻度 + xlabel）
D 行：  10（行距）+ 7.2（字母带） + 100（面板高） + 30（x 刻度 + xlabel）
```

A–C 是三个 120 pt 方形面板并排（面板外框 = axes 120 + 实测 y 轴装饰 ≈29 + 右边距 10 ≈ 159 pt，三列共 477 pt，占满 493 pt 页宽）；D 在第二行由 `mp.newline()` 从行首开始，用 house 的 SQUARE 尺寸 100 × 100 pt（单 panel 必须落在 80–150 pt 内，不能通栏），右侧留白。没有外部矢量面板，因此不需要 pypdf 叠加。

**版式机制**：`mp.panel(width, height)` 只给 axes 尺寸；每列的横向宽度和字母位置在渲染时按实测的 y 轴装饰宽度算出，所以 A/B/C 的字母都落在各自列的外框左缘，A 与 D 的字母因此都在 x = 0（左对齐）。ylabel 变长只会把它自己这一列往右推，不会压到前一个面板。

**产物**：

| 文件 | 说明 |
|---|---|
| `results/7a_figures/figure_s1.pdf` | 最终矢量图（单页） |
| `results/7a_figures/figure_s1.review.png` | 360 dpi 评审位图 |

**为什么复制而不是复用**：这三张图原本由 `workflow/scripts/pcr_qc/plot_pcr_qc.py` 画在 1a 的 2×2 图里。本脚本自带一份密度着色、对数轴与 n/r/P 标注的实现（KDE 部分对应 `workflow/src/figure_render/_point_density.py`，其余对应 `_scatter.py` 的密度路径），不 import 任何 stage 模块——改 S1 的版式不会动 1a 的 QC 图。

**输入为什么是 `_work/` 下的 parquet**：1a_pcr_qc 对 A–C 三个面板只发布了一张 2×2 PDF，没有发布表；`pbl_pbr / tech / bio.parquet` 就是它准备的表，rule 把它们声明为 input，1a 重算时 S1 会自动跟着重建。

**D 面板的输入**：HD_DIT_HAP release 的 `insertion_level/annotations.tsv.gz`，只读 `Type`、`Systematic ID` 两列；路径在 `7a_figures.smk` 里由 `DATASETS` 拼出。没有基因内插入的 56 个基因不出现在注释表里，因此分布覆盖的是有 ≥1 条插入的 5,070 个基因。（把 1b 的 5,126 基因全集也接进来会让 8 个特征收集任务重跑，为 1.1% 的零值不值得。）

---

### Panel A — PBL vs PBR（同一文库两个引物）

| 要素 | 内容 |
|---|---|
| 输入 | `results/1a_pcr_qc/_work/pbl_pbr.parquet` |
| 必需列 | `PBL`、`PBR`（缺列在 I/O 边界直接报错） |
| 变换 | 只保留两列都 > 0 的行（对数轴画不出 0）：164,820 → 108,717 行 |
| 画法 | 在 log10 空间做 2D 核密度（`gaussian_kde`，>5000 点时固定种子抽 5,000 点、128×128 网格插值），逐点按密度排序后 `ax.scatter(s=1.5, alpha=0.6, linewidths=0, rasterized=True)`，viridis 着色并归一化到 [0,1]；灰色虚线 x=y；左上 `n / r / P`（log10 空间 Pearson，`P < 1e-9` 显示为 `P < 1E-9`，否则 `P = mantissa E±exp`）；右下 inset 密度色条（0→low、1→high，标题 Density）；双对数轴（本图中 matplotlib 只给 y 轴次刻度，按 house 规则两侧一并关闭，与 1a 的 QC 图一致） |
| 版式 | axes `w=120, h=120` pt；由 flow 排到 `x=28.8, y_top=7.2`（x 含实测 ylabel 宽度）；字母 A 由 multipanel 放在列首，x = 0 |
| 产物 | 并入 `figure_s1.pdf` |

### Panel B — 技术重复

| 要素 | 内容 |
|---|---|
| 输入 | `results/1a_pcr_qc/_work/tech.parquet`（同一份样品在两个上游项目各跑一次，按 (Chr, Coordinate, Strand) 合并） |
| 必需列 | `Reads_1`、`Reads_2` |
| 变换 | 两列都 > 0：152,276 → 117,297 行 |
| 画法 | 与 A 相同；轴标签 `Reads of Technical Replicate 1 / 2` |
| 版式 | axes `w=120, h=120`，紧随 A 排到 `x=187.6`（比 A 宽 159 pt，因为 ylabel 更长）；字母 B 在列首，x ≈ 158.8 |
| 产物 | 并入 `figure_s1.pdf` |

### Panel C — 生物学重复

| 要素 | 内容 |
|---|---|
| 输入 | `results/1a_pcr_qc/_work/bio.parquet`（同一项目内两个样品） |
| 必需列 | `Reads_1`、`Reads_2` |
| 变换 | 两列都 > 0：125,566 → 92,704 行 |
| 画法 | 与 A 相同；轴标签 `Reads of Biological Replicate 1 / 2` |
| 版式 | axes `w=120, h=120`，排到 `x=346.5`，右缘 466.5 pt（页宽 493，行未占满）；字母 C 在列首，x ≈ 317.6 |
| 产物 | 并入 `figure_s1.pdf` |

### Panel D — 每个基因的插入条数分布

| 要素 | 内容 |
|---|---|
| 输入 | `…/HD_DIT_HAP/release/insertion_level/annotations.tsv.gz`（release 的逐插入注释表） |
| 必需列 | `Type`、`Systematic ID`（缺列在 I/O 边界直接报错） |
| 变换 | 只数 `Type == "Coding gene"` 的插入行；按 `Systematic ID` 计数得到 5,070 个基因的条数（中位 62，四分位 36–98，最大 1,254）；再按 25 条一档分箱，0–299 共 12 档，≥300 合并为最后一档（46 个基因）。第一档的刻度是 bin 左边界 `0`，实际装的是 1–24 条（零插入基因不在表里） |
| 画法 | `cns.barplot(x="bin_label", y="n_genes", color=palette[2])`，柱色为 house palette 第 2 位（灰度可分）；13 个刻度在 100 pt 里放不下，`ax.set_xticks(...[::3])` 隔两档显示（0 / 75 / 150 / 225 / ≥300）；x `Insertions per gene`、y `Number of genes` |
| 版式 | axes `w=100, h=100` pt（house SQUARE），排在 `x=33.6, y_top=174.4`；字母 D 在列首，x = 0（与 A 左对齐） |
| 产物 | 并入 `figure_s1.pdf` |

**为什么是普通直方图而不是累计曲线**：横轴跨 0–1254，中位 62、99% 的基因 ≤ 290；线性横轴上的累计曲线会把 99% 的曲线挤进左侧 24% 的宽度。普通直方图能直接读出众数（25–75）与长尾，最后一档标成 `≥300`，避免被误读成又一个 25 宽的 bin。

---

### 可调旋钮（全部在 `workflow/scripts/figures/figure_s1.py` 顶部）

| 常量 | 当前值 | 作用 |
|---|---:|---|
| `CANVAS_WIDTH_PT` | 493 | 画布宽（Cell 双栏），即 `multipanel(max_width=...)` |
| `PANEL_SIZE_PT` | 120 | A–C 每个方形面板的 axes 边长 |
| `BELOW_DECORATION_PT` | 30 | x 刻度 + xlabel 所占高度（`margin_bottom`） |
| `ROW_GAP_PT` | 10 | A–C 行与 D 行之间（D 的 `margin_top`） |
| `PANEL_D_SIZE_PT` | 100 | D 的 axes 边长（80–150 pt 契约） |
| `COUNT_BIN_WIDTH` / `COUNT_TAIL_START` | 25 / 300 | 直方图档宽与最后一档的起点 |
| `DENSITY_SCATTER_S` / `DENSITY_SCATTER_ALPHA` | 1.5 / 0.6 | 点大小与透明度（密云例外） |
| `DENSITY_CBAR_BOUNDS` | (0.82, 0.02, 0.035, 0.30) | 密度色条在 axes 分数坐标中的位置与大小 |

改 `PANEL_SIZE_PT` 后列宽和画布高度自动重算。三列总宽 = Σ(axes + 实测装饰 + `margin_right` 10)，当前 3 × 158.8 ≈ 476 pt；若把它加到 ~130 pt，C 会被 `max_width` 换行到第二行，改前先跑一次看渲染。面板间距用 `mp.panel(..., margin_right=...)` 调。

---

## Supplementary figures

（其余 S 图待加入：每个 S 图同样一条 rule、一节说明。）
