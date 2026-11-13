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

## Figure 1 — DIT-HAP 流程、覆盖度、spike-in、LFC 轨迹、FACS 验证与聚合算法

**画布**：A4 整页 595.28 × 841.89 pt（210 × 297 mm）。`cns.multipanel` 只决定宽度，所以内容画布宽 515.28 pt（= 595.28 − 2 × 40 pt 页边距）、高由四行面板推出（实测 658.1 pt）；组装时 pypdf 把内容页贴到 A4 页 (40, 40) pt 处（距左、距上各 40 pt，右侧也留 40 pt），页面其余部分留白。内容画布本身由 `cns.multipanel(max_width=515.28, unit="pt")` 流式排版，三个 `mp.newline()` 隔开四行：

```
A 行：     10（字母带 = 7.2 标签高 + 2.8 净空） + 141.6（示意图页面实际高度） + 6（行距）            = 157.6
B/C/D 行： 7.2（字母带） + 17.3（图例，实测） + 100（面板高） + 26（下方装饰）                        = 150.50
E 行：     6（行距） + 7.2（字母带） + 6.6（图例，实测） + 171（容器高） + 0（行距由 F 的 margin_top 给） = 190.80
F 行：     6（行距） + 7.2（字母带） + 120（容器高） + 26（下方装饰：x 刻度 + 标题）                  = 159.20
            内容高合计 658.1，贴到 841.89 高的 A4 页上距顶 40 pt 处
```

B/C/D 用同一个 100 × 100 pt 方形 axes，三块绘图区大小一致、长宽相等（house 面板契约 80–150 pt 内）；A4 宽度下这一行右端会留出空白。E 是一个 509.28 pt 宽的空白容器 axes，内部再排 4 列 × 2 行的 mini axes（每个 80 × 60 pt，字号沿用 house 的 8/7 pt）；图例一行放在容器上方。F 是一个 **400 × 120 pt** 的空白容器 axes，左对齐、右侧留白；内部按点坐标排五个 inset axes（左侧镜像剖面、中间过渡、右侧聚合曲线）。F 的 x 刻度标签与 x 标题放在容器**下方**的 26 pt 装饰带里（与 B/C/D 的 `margin_bottom` 同理），所以 120 pt 容器内只放绘图区本身。

**版式机制**：`mp.panel()` 的 `width/height` 是 axes 尺寸，面板的外框在渲染时再加上**实测**的 y 轴装饰宽度（刻度标签 + ylabel）与右边距 10 pt。所以 B/C/D 从左往右流式排列、行首贴页面左边；ylabel 变长只会把它自己这一列往右推，不会压到前一个面板。面板字母由 multipanel 放在每列外框的左缘（8 pt 粗体），A、B、E、F 的字母因此都在 x = 0。E、F 的 inset axes 用 `container.inset_axes()` 建，位置是容器 axes 的分数坐标，渲染时自动跟随容器。

**产物**：

| 文件 | 说明 |
|---|---|
| `results/7a_figures/figure1.pdf` | 最终矢量图：单页 A4，内容距左、上、右各 40 pt |
| `results/7a_figures/figure1.review.png` | 360 dpi 评审位图，由最终 PDF 栅格化（含 A 面板） |
| `results/7a_figures/_work/figure1_canvas.pdf` | 叠加示意图前、贴 A4 前的内容画布（A 区空白，515.28 × 470.7 pt）；排查版式时看这个 |

**为什么 A 面板用 pypdf 叠加**：matplotlib / cnsplots 没有把 SVG/PDF 当矢量嵌入的 API，`ax.imshow` 一定栅格化。脚本先用 matplotlib 画整张画布（A 区留白），再用 `pypdf` 把示意图 PDF 等比缩放、平移进 A 区，最后 `merge_transformed_page` 合成单页 PDF。示意图的 4 个 Arial 子集与 4 张 601 ppi 位图都原样保留。

---

### Panel A — DIT-HAP 工作流示意图（外部绘制）

| 要素 | 内容 |
|---|---|
| 输入 | `resources/figures/figure1a.pdf`（单页，页面框即内容框，510 × 141.6 pt） |
| 变换 | 不缩放（×1.0000）：预留框就是示意图页面框 510 × 141.6 pt，pypdf 原样平移到该框，不加边、不缩放 |
| 画法 | multipanel 建一个空 axes（`set_xticks([]) + set_yticks([]) + set_axis_off()`，清刻度是为了让 multipanel 量到 0 宽的左装饰）；`pypdf` 的 `PdfWriter(clone_from=canvas)` + `merge_transformed_page(schematic, Transformation().scale(s).translate(x,y))`；坐标原点是 PDF 左下角 |
| 版式 | 预留框 `x=5.775, y_top=10, w=510, h=141.6`（示意图页面框原尺寸；右缘比 515.28 pt 的内容画布多约 0.5 pt，pypdf 合成时会被页面框裁掉这一条）；不裁边、不补白边 |
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
| 画法 | 一个 509.28 pt 容器 axes（`set_axis_off()`）内嵌 8 个 `inset_axes`，4 列 × 2 行、按基因名排序（每行 4 个）。每格画 3 条线：DIT-HAP（palette 4 绿）、FACS `Strain#1`（palette 5 鲑红）、FACS `Strain#2`（palette 6 紫）——两条是独立菌株而非技术重复，所以都画，各用一色、不用线型区分；三条数据线线宽 2 pt（用户指定，覆盖 house 默认 1.5）、`alpha=0.8`，图例色块保持不透明；零线灰色虚线 0.5 pt；基因名斜体加粗作 title（house 8 pt） |
| 轴 | 全格共享 x `0–12.8`（刻度 0/4/8/12）与 y（刻度间隔 2、向外取整）；只有第一列显示 y 刻度标签、只有底行显示 x 刻度标签；两个共享轴标题用容器 `text()` 居中在整个网格上，不用某个 mini axes 的 label |
| 图例 | 一行三列、无边框，挂在容器 axes 上方（`loc="lower center"`，锚点 x = 网格水平中心 = 209/509.28 ≈ 0.411），水平居中于整个网格：`DIT-HAP`、`FACS Strain#1`、`FACS Strain#2`；multipanel 把它量作 E 的顶部装饰（实测高 6.6 pt、宽 152.6 pt） |
| 版式 | mini axes `80 × 60` pt（4 列共 350 pt，加左侧 34 pt 后仍留白）；列距 10 pt、行距 16 pt（给下一行 8 pt title 留位）；左侧 34 pt（y 刻度 + 旋转 y 标题）、底部 22 pt（x 刻度 + x 标题）、顶部 13 pt（首行 title）；容器高 = 13 + 2×60 + 16 + 22 = 171 pt；图例在容器外上方另计；E 的 `margin_bottom=0`（到 F 的行距由 F 的 `margin_top` 给） |
| 产物 | 并入 `figure1.pdf` |

数据口径：FACS 表的 `Generations` 每株各自记录（同名的 YES 时间点不等于同一世代值），因此不能共用 x 向量；`Spore` 恒为 0，有三个 block 的 YES0 略小于 0（如 `arg6/Strain#1` 的 −0.18），x 下限留到 −0.3。

---

### Panel F — 从 insertion level 到 gene level 的加权聚合（cdc7 实例）

| 要素 | 内容 |
|---|---|
| 输入 | `<HD_DIT_HAP release>/insertion_level/LFC.tsv`（每插入方向一行 × YES0…YES4）、`insertion_level/padj.tsv`（同形状，画棒棒糖权重）、`insertion_level/insertion_weights.tsv`（长表：插入 × 时间点的基因内归一化权重，同基因-时间点求和为 1）、`gene_level/LFC.tsv`（右侧珊瑚色基因级曲线）、`resources/external/pombase/<version>/genome_region/coding_gene_primary_transcripts.bed`（primary transcript 的 CDS/intron 分段） |
| 必需列 | release 三表：`Chr`、`Coordinate`、`Strand`、`Target` + `YES0`…`YES4`，weights 另加 `Timepoint`、`Weight`、`Systematic ID`；BED：`#Chr`、`Systematic ID`、`Feature`、`Start`、`End`、`Strand`、`ParentalRegion_start`、`ParentalRegion_end`（缺列在 I/O 边界直接报错） |
| 示例基因 | 固定在 `AGGREGATION_GENE_ID = SPBC21.06c`（cdc7）；选它而不是逐次挑：在 inviable 基因里，它是唯一同时有 ≥10 个测得 TTAA 位点、跨过 ORF 10/10 个 10% 分箱、且有 ≥3 个末端位点在 12.5 代仍停在 0 附近（LFC > −0.5）又带内含子的基因。换基因只改这一个常量，但聚合核对（见下）会要求该基因在 gene-level 表里存在 |
| 变换 | ① 三条 release 表（LFC / padj / weights）按 4 级插入键 + `Timepoint` 拼成一张长表；每个位点的两个 TTAA 方向是独立测量，但共享同一物理坐标，故按 `Coordinate × Timepoint` 合并：`site_weight = ΣWeight`、`site_lfc = Σ(Weight·LFC)/site_weight`（权重加权平均，合并前后基因级总和严格不变）。② 每个位点的 `raw_weight = −log10(clip(padj, 1e-6, 1−1e-6))`，取 YES4（12.5 代）并对其两个方向取均值，作为棒棒糖高度；`padj < 0.05`（即 `raw_weight > 1.301`）画高权重色，其余画灰。③ 坐标换到 ORF 坐标：plus 链 `coord − ParentalRegion_start`、minus 链 `ParentalRegion_end − coord`，所以 5′ 在左、3′ 在右且保留真实 bp 间距。④ 右侧灰色轨迹用 YES0…YES4 全五个时间点，珊瑚色曲线直接取 release 的 gene-level 行，不重算发布的数字 |
| 聚合核对 | 画之前把 `Σ(site_weight · site_lfc)` 与 release 的 gene-level 行逐时间点比较；偏差 > `AGGREGATION_TOLERANCE = 1e-2` 直接报错（实测 3.63e-4，覆盖 release 的三位小数舍入）。这一步只防图内 join 错位，不改变任何发布值 |
| 画法 | **左**：上下两层共享 x、无垂直间距；上层 `vlines + scatter` 棒棒糖（高权重 palette 1 青 `#2f7e8f`、低权重 house 灰 `#A3A3A3`，竖线 **1.0 pt**、圆点 **20 pt²**），y 0–6.4（权重上限 6 是上游 clip 的结果，多留 0.4 给点）；下层四条空间剖面折线（每个世代一条，按位点坐标顺序连），**同一珊瑚色相加白从浅（2.4 代，35% 色）到深（12.5 代，100% 色）**，线宽 **2 pt**，不画点；最下方基因结构条 = 一条 0.5 pt 横线（intron）+ CDS `barh` 实心块（**palette 3 蓝 `#4e5a8a`**，块高 0.5 的 axes 高度 ≈ 7.5 pt），5′/3′ 标在 bar 的左右两端、与 bar 共用 `GENE_MODEL_BASELINE_Y` 并 `va="center"`（严格纵向居中）；基因名 `cdc7` 斜体加粗、在左侧 120 pt 块上**水平居中**（`set_title(loc="center")`）。**中**：黑色右箭头 + 上方 `Weighted aggregation` + 下方单行公式 `LFC_gene,t = Σ(w_i·LFC_i,t)/Σw_i`（mathtext 7 pt，斜杠式而非堆叠分式，好让 120 pt 高的中间块与上下图例各留约 10 pt 等距）。**右**：每个位点一条灰轨迹（**1.0 pt**、alpha 0.3，含 YES0），上面一条珊瑚粗线（palette 0，**3.0 pt**）画 gene-level；零线灰虚线 |
| 图例 | 三组、全部无边框；中间块自上而下：① `High weight`/`Low weight`（weight axes 右侧，**一列两行**，图例上缘与 weight y 轴顶端齐平）→ ② `Weighted aggregation` + 箭头 + 公式（公式块与上下两个图例各留约 **10 pt** 等距）→ ③ `Generations`（过渡 axes **左下、一列**，标题在上、四个色块竖排，左缘与 ① 左对齐）；`Insertion-level LFC`/`Gene-level LFC` 在 aggregate axes **顶部、一列两行**，两行中线正好落在 y 轴顶端（图例跨出 axes 上缘），水平锚点 **0.7**（偏右） |
| 轴 | 左上层 `Weight (−log10 padj)`、y 刻度 0/2/4/6；左下层 `LFC (log₂)`、y 刻度 0/−2/−4/−6/−8；两层在 y = 0 处共用一条**可见的 x 轴线**（weight axes 的 bottom spine，0.5 pt），profile 不再画灰色虚线零线；x 刻度短线全部去掉，只有刻度标签（0/1,000/2,000/3,000，千分位）出现在基因条**下方**的 26 pt 装饰带里，下面接 `Position along gene (bp)`；右图 x `Generations`（0/2.4/5.6/9.1/12.5）、y `Gene-level LFC (log₂)`，刻度仍为 0…−8，但 y 轴上限留到 **+1.5**，让 y = 0 虚线不贴顶、也给跨出上缘的两行图例留出空间（0 以上无刻度） |
| 版式 | 容器 **400 × 120** pt；左侧装饰 30 pt、左 axes 宽 **120** pt；上层 **120 × 30**、下层 **120 × 75**、基因条 **120 × 15**（三层相加正好 120）；中间过渡区 **100 × 120**、右 axes **120 × 120**（左缘 x = 280）；三块之间的横向间距为 0（靠 inset 内部留白分开）；x 刻度标签 + 标题放在容器下方 **26 pt** 的 `margin_bottom` 里；所有几何常量在脚本顶部，改一个轴宽/高，右侧起点与画布高度自动重算 |
| 产物 | 并入 `figure1.pdf` |

数据口径：权重来自上游 `-log10(padj)` 再按基因-时间点归一化；棒棒糖画的是**归一化前**的 `−log10(padj)`（YES4），所以高度上限恰好是 6——这是 `padj` 被 clip 到 `1e-6` 的结果，不是坐标轴人为截断。两个插入方向的 `Weight` 相加即该位点的归一化权重；`site_lfc` 是这两个方向的权重加权均值，因此 `Σ(site_weight · site_lfc)` 与逐方向求和完全等价，右侧珊瑚曲线可由灰色轨迹严格复现。

---

### 面板字母

六个字母由 `mp.panel(label=...)` 自动放置：8 pt 粗体 Arial，每列外框左缘、面板顶端上方（A 用 `pad_top` 留 2.8 pt，D 用 `pad_top` 抵消图例下半，其余为 0）。`multipanel` 把它放在"实测左装饰外侧"，而列的外框左缘就是页面的列起点，所以 A、B、E、F 的字母都在 x = 0，不随 ylabel 宽度漂移。**不要**再手摆 `fig.text` 字母——会和 multipanel 的自动标签重影。

---

### 可调旋钮（全部在 `workflow/scripts/figures/figure1.py` 顶部）

| 常量 | 当前值 | 作用 |
|---|---:|---|
| `A4_WIDTH_PT` / `A4_HEIGHT_PT` | 595.28 / 841.89 | 最终页面尺寸（210 × 297 mm） |
| `PAGE_MARGIN_PT` | 40 | 内容页距 A4 页左、上、右的留白 |
| `CANVAS_WIDTH_PT` | 515.28 | 内容画布宽（= A4 宽 − 2 × 页边距），即 `multipanel(max_width=...)` |
| `PANEL_LABEL_COLUMN_PT` | 6 | 给每列字母留的那一列；`PAGE_CONTENT_WIDTH_PT = 515.28 − 6`（E、F 用它做宽度；A 用示意图自身的页面框） |
| `PANEL_LETTER_PAD_TOP_PX` | 14 | A 上方字母净空（显示像素，≈2.8 pt @360 dpi） |
| `ROW_GAP_PT` | 6 | A 行与 B/C/D 行之间（A 的 `margin_bottom`）；E、F 的行距（各自的 `margin_top`） |
| `ROW_PANEL_SIZE_PT` | 100 | B/C/D 共用的方形 axes 边长（x 由 flow 决定） |
| `ROW_BELOW_DECORATION_PT` | 26 | x 刻度 + xlabel 所占高度（`margin_bottom`） |
| `BAR_COUNT_FONT_PT` | 5 | 柱端数字字号（fit-critical 覆盖，house 允许） |
| `LFC_LEGEND_HALF_HEIGHT_PT` | 8.4 | D 的图例跨过 axes 顶边时，`pad_top` 补偿的下半高度 |
| `FACS_GRID_COLUMNS` | 4 | E 的列数；行数由 `ceil(基因数 / 4)` 自动推出 |
| `FACS_MINI_AXES_WIDTH_PT` / `FACS_MINI_AXES_HEIGHT_PT` | 80 / 60 | E 每个基因 mini axes 的尺寸（80×4 + 10×3 + 34 = 384，图例在上方另占一行） |
| `FACS_EXCLUDED_GENES` | `orc6`、`srb4` | 不画的基因（用户 2026-11-12 决定），在 frame 构建时剔除 |
| `FACS_LINE_WIDTH_PT` | 2 | E 的三条数据线线宽（用户指定） |
| `FACS_COLUMN_GAP_PT` / `FACS_ROW_GAP_PT` | 10 / 16 | E 的列距/行距；行距要给下一行基因名 title 留位 |
| `FACS_TITLE_BAND_PT` | 13 | E 首行基因名 title 的顶部留白 |
| `FACS_BOTTOM_DECORATION_PT` | 22 | E 的 x 刻度标签 + 共享 x 标题 |
| `FACS_LEFT_DECORATION_PT` | 34 | E 的 y 刻度标签 + 共享 y 标题 |
| `FACS_LFC_SIGN` | −1 | FACS LFC 符号翻转，与 release 的"耗竭为负"对齐 |
| `AGGREGATION_GENE_ID` | `SPBC21.06c` | F 的示例基因（cdc7，固定不做逐次选择） |
| `AGG_PANEL_WIDTH_PT` / `AGG_PANEL_HEIGHT_PT` | 400 / 120 | F 容器尺寸（左侧三层相加正好 120） |
| `AGG_LEFT_DECORATION_PT` / `AGG_LEFT_AXES_WIDTH_PT` | 30 / 120 | F 左侧 y 装饰宽、镜像剖面的 axes 宽 |
| `AGG_WEIGHT_HEIGHT_PT` / `AGG_LFC_HEIGHT_PT` / `AGG_GENE_AXES_HEIGHT_PT` | 30 / 75 / 15 | 权重区 / LFC 区 / 基因条高度（30 : 75 = 1 : 2.5，合计 120） |
| `AGG_TRANSITION_WIDTH_PT` / `AGG_TRANSITION_HEIGHT_PT` | 100 / 120 | 中间过渡区尺寸 |
| `AGG_RIGHT_DECORATION_PT` / `AGG_RIGHT_X_PT` | 30 / 280 | 右侧 y 装饰宽、aggregate axes 左缘 |
| `AGG_BELOW_DECORATION_PT` | 26 | F 的 x 刻度标签 + x 标题装饰带（`margin_bottom`） |
| `AGGREGATION_WEIGHT_LINE_WIDTH_PT` / `AGGREGATION_WEIGHT_MARKER_SIZE` | 1.0 / 20 | 棒棒糖竖线宽、圆点面积（pt²） |
| `AGGREGATION_PROFILE_LINE_WIDTH_PT` / `AGGREGATION_INSERTION_LINE_WIDTH_PT` / `AGGREGATION_GENE_LINE_WIDTH_PT` | 2.0 / 1.0 / 3.0 | 左下层剖面、右侧灰轨迹、右侧 gene-level 曲线线宽 |
| `GENERATION_COLOR_LIGHT_MIX` / `GENERATION_COLOR_DARK_MIX` | 0.35 / 1.0 | 世代单色珊瑚加白比例（浅 → 深） |
| `GENE_MODEL_COLOR` | palette 3 `#4e5a8a` | 基因结构条的蓝色（用户指定） |
| `AGGREGATE_HEADROOM` | 1.5 | 右图 y = 0 上方的空带（给两行图例；刻度仍到 0） |
| `AGGREGATION_TOLERANCE` | 1e-2 | 图内加权聚合与 release gene-level 行的最大允许偏差 |

改 axes 尺寸后画布高度自动重算；B/C/D 不会再被"摆到中间"，它们总是从行首左对齐。若想改变面板之间的空隙，调 `mp.panel(..., margin_right=...)`（默认 10 pt）。

---

## Figure 2 — 耗竭动力学：一条拟合曲线如何定义 A、DR、DL

**画布**：A4 整页 595.28 × 841.89 pt（210 × 297 mm）。内容画布宽 515.28 pt（= 595.28 − 2 × 40 pt 页边距）、高由两行面板推出（实测 382.4 pt = A/C 行 7.2 + 100 + 行距 10 + B 行 7.2 + 258）；组装时 pypdf 把内容页贴到 A4 页 (40, 40) pt 处（距左、距上、距右各 40 pt），页面其余部分留白。

**当前是 A、B、C 三个面板**（2026-11-13 起逐步搭建，其余面板待加入）。第一行左侧 A 是"模型定义"面板：数据 → 公式 → 参数的分析流程，右侧用 mtf1 的真实拟合把三个参数画成几何元素；A 右侧的 C 是 DR 与 DL 的分布（一行两列直方图）；B 在 A 正下方、等宽，用 2×2 的四个基因展示同一模型在真实数据上的拟合。

**产物**：

| 文件 | 说明 |
|---|---|
| `results/7a_figures/figure2.pdf` | 最终矢量图：单页 A4，内容距左、上、右各 40 pt |
| `results/7a_figures/figure2.review.png` | 360 dpi 评审位图，由最终 PDF 栅格化 |
| `results/7a_figures/_work/figure2_canvas.pdf` | 贴 A4 前的内容画布（515.28 × 382.4 pt）；排查版式时看这个 |

---

### Panel A — 模型定义：数据 → Gompertz 公式 → 三个参数（mtf1 示例）

**版式**：一个 `236 × 100` pt 的容器 axes（字母 A 由 multipanel 放在列首）。容器内部分三块：

```
x 2–80    左侧流程列（78 pt 宽）：数据框 / 公式框 / 参数框，竖直箭头连接
x 80–101  水平箭头（约 20 pt），指向右侧曲线
x 101–126 解剖图的 y 轴装饰（y 刻度 + 旋转 ylabel）
x 126–226 解剖图 axes（100 × 76 pt）
x 226–236 右边距（给最右 x 刻度标签留位）
```

| 要素 | 内容 |
|---|---|
| 输入 | `<HD_DIT_HAP release>/gene_level/gene_level_fitting_statistics.tsv`（release 的基因级拟合表，一行一个基因） |
| 必需列 | `Name`、`time_points`、`A`、`DR`、`DL`、`R2`、`t_inflection`、`y_inflection`、`YES0`…`YES4`、`YES0_fitted`…`YES4_fitted`（缺列在 I/O 边界直接报错） |
| 变换 | 按 `Name == "mtf1"` 取唯一一行（缺失或重复都报错）；`time_points` 的逗号串解析为世代向量 `0 / 2.352 / 5.588 / 9.104 / 12.48`；平滑曲线由发布的 `A/DR/DL` 代入模型求值（300 点）；`*_fitted` 五列只用来核对模型一致（逐点偏差 ≤ 5e-3，覆盖参数舍入的传播），不用于摆线。mtf1：A = −5.509、DR = −0.656、DL = 4.606、R² = 0.999、t_inflection = 7.698、y_inflection = −2.027 |
| 流程框 | 圆角框 `78 × 17 / 33.5 / 19.5` pt，浅灰底 `#f7f7f7`、灰色 0.6 pt 描边；框间竖箭头（灰色 0.8 pt，`mutation_scale=6`）；顶部框文字 `Time-series LFC data` 6.5 pt；底部框 `Extracted parameters` 6.5 pt + 彩色参数名 `DR / DL / A` 7 pt 加粗，横向居中 |
| 公式 | 放在中间框，两行 6.5 pt、行距 9 pt、各自水平居中：① `LFC(t) = A·exp(−exp[`；② `(DR·e/A)(DL−t) + 1])`。参数名 `A / DR / DL / A` 带同色 25% alpha 圆角底（`FancyBboxPatch`）；公式下方灰色 `R² = 0.999` 6.5 pt。**实现注意**：matplotlib 的 mathtext 不支持 `\colorbox`，所以公式被拆成若干 `FormulaPiece`，用 Agg renderer 逐段量宽后在框内左→右排，参数段再按实测 bbox 画圆角底；全程矢量、可进 Illustrator 编辑 |
| 参数配色 | `A` = Cell 槽 2 赭黄 `#e1a22e`；`DR` = 槽 0 砖红 `#c84c3a`；`DL` = 槽 1 青绿 `#2f7e8f`。数据与拟合线用中性色：曲线 = 槽 3 蓝灰 `#4e5a8a`（1.5 pt），实测点 = 槽 9 灰 `#6b6b6b`（`s=12`、白边 0.4 pt）。颜色只承担"参数含义"，不与数据语义冲突 |
| 解剖图 | 实测 5 点 + 拟合曲线；灰色点线 `y=0` 基线。三条 guide 全部用参数色：① **A**：`y=A` 横虚线（4,2），右端圆点，标签 `A (plateau)` 在线的左端上方；② **DR**：拐点切线 `y=DR·(x−DL)`（点划线 5,1.4,0.8,1.4），在交点处画圆点 `(t_inflection, y_inflection)`，标签 `DR (slope)` 在切线上方；③ **DL**：`x=DL` 竖点线（1,1.6）从 `y=0` 到轴底，两端各一个圆点，标签 `DL (intercept)` 在顶端右侧。所有 guide 线宽 1 pt、圆点 `s=16` |
| 轴 | x `Generations`，刻度为五个世代值（一位小数，`2.352 → 2.4`、`12.48 → 12.5`，末尾 `.0` 去掉），6.5 pt；y `LFC (log₂)`，7 pt，刻度只取数据范围内的 2 单位网格（mtf1 为 −4/−2/0），6.5 pt |
| 比例 | y 轴数据下限**包含 A**（即使曲线在观测窗口内还没到达平台，否则平台线会贴到轴上）；下留白 28%、上留白 31%（上方给切线和标签）；x 两侧各 3% |
| 字体 | 流程/公式/刻度 6.5 pt；参数名与轴标题 7 pt；全部 Arial（TrueType 嵌入） |
| 产物 | 并入 `figure2.pdf` |

数据口径：面板里**不写基因名**（这是模型定义面板，不是基因故事），示例基因 mtf1 的记录在本文档与图注；实测点保留，用来连接左侧"Time-series LFC data"这一步。`A` 是平台值（负值 = 耗竭到底）、`DR` 是拐点处切线的斜率（最大耗竭速率）、`DL` 是切线与 `y=0` 的交点（延迟）。mtf1 的曲线在 12.48 代时约 −4.45，尚未到达 A = −5.51，平台线因此画在数据末端之下。

---

### Panel B — 四个基因的示例拟合（2 × 2）

**版式**：一个 `236 × 258` pt 的容器（与 A 等宽），内部 2 列 × 2 行 mini axes，每个 `100 × 76` pt（与 A 右侧 plot 同一 footprint），下方 58 pt 是参数表带。行主序：`tna1`、`pst2` 在上排，`rpl3801`、`nrs1` 在下排。左装饰 25 pt（y 刻度 + 共享 ylabel）、列距 6 pt、行距 14 pt（给下排基因名 title 留位）、顶部 title 带 12 pt、底部 x 装饰 22 pt；网格下方 10 pt 净空 + 46 pt 表格（表头 10 + 4 行 × 9）+ 2 pt 底边。右侧留 5 pt 给最右 x 刻度标签。

| 要素 | 内容 |
|---|---|
| 输入 | 与 A 同一张 release 表，按 `Name` 取四行 |
| 基因 | `tna1`、`pst2`、`rpl3801`、`nrs1`（`PANEL_B_GENES` 的顺序就是行主序） |
| 变换 | 每个基因各自 `build_gene_fit`：解析世代串、用发布的 A/DR/DL 求 300 点曲线、`*_fitted` 五列核对（≤ 5e-3）。四行都必须存在，缺一行直接报错 |
| 轴 | **共享** x 与 y。x 用 B 自己的留白：左 9%、右 3% → −1.12…12.85（左侧多留的 0.75 代是为了让 `DL = 0` 的竖线离左轴脊足够远，三个基因的 DL 都是 0）；y 数据范围取四基因的 min（含各自 A）到 max，四个基因合并后为 −6.677…0.03，下留白 8%、上留白 12% → 约 −7.21…0.83；y 刻度 −6/−4/−2/0，2 单位网格；x 刻度为五个世代值。只有左列显示 y 刻度标签、只有底行显示 x 刻度标签；共享 `LFC (log₂)` 在网格左侧居中、共享 `Generations` 在网格底部居中 |
| 画法 | 每格画拟合曲线（槽 3 蓝灰、1.5 pt）、实测点（槽 9 灰、`s=12`、白边 0.4 pt）、灰色点线 `y=0` 基线；**每格都带 A 的三条 guide**（A 平台线 + 右端圆点、DR 拐点切线 + 拐点圆点、DL 竖点线 + 两端圆点），但**不重复 guide 的文字标签**——标签只在 A 出现，四格靠颜色 + 下方参数表对应。B 的 guide 比 A 轻：线宽 0.8 pt、圆点 `s=8`（A 为 1 pt / `s=16`） |
| 参数表 | 网格下方一张表：行 = 四个基因（斜体黑字），列 = `Gene`、`A`、`DR`、`DL`、`R²`。表头与数值按参数配色（A 赭黄 / DR 砖红 / DL 青绿、R² 灰），数值右对齐（A/DR/DL 两位小数、R² 三位小数、`-0.0` 显示为 `0.00`），表头下一条 0.4 pt 灰横线；列锚点 x = 42 / 111 / 144 / 177 / 214 pt |
| 标题 | 每个 mini axes 用基因名作 title（house 8 pt，斜体加粗） |
| 产物 | 并入 `figure2.pdf` |

数据口径（四个基因的发布拟合）：`tna1` A = 0.025、DR = 0.003、DL = 0.000、R² = −0.225——这是一条退化的近平线（实测点都在 0 ± 0.05），在面板里就是"不耗竭"的例子，曲线贴着零线；`pst2` A = −2.901、DR = −0.200、DL = 0.830、R² = 0.993；`rpl3801` A = −3.963、DR = −0.454、DL = −0.000、R² = 0.979；`nrs1` A = −6.677、DR = −1.119、DL = −0.000、R² = 0.968。共享 y 轴让四者的耗竭深度可直接比较（tna1 ≈ 0 → nrs1 ≈ −6.7）。

---

### Panel C — DR 与 DL 的分布（一行两列）

**版式**：一个 `252 × 100` pt 的容器，放在 A 同一行的右侧（A 236 + 列距/字母留白后仍 < 515.28）。容器内两个 `100 × 76` pt 的 mini axes，左 DR、右 DL，间距 10 pt；左装饰 34 pt（四个计数刻度 + 共享 ylabel）、右侧留 8 pt。两个面板**共享一个 count y 轴**（同一批 4,513 个基因，DL 的 0 峰 3,113 与 DR 的峰 1,263 可直接比较），只有左图显示 y 刻度标签。

| 要素 | 内容 |
|---|---|
| 输入 | 与 A/B 同一张 release 表，只读 `Status`、`DR`、`DL` |
| 筛选 | `Status == "Success"`；HD_DIT_HAP 的 4,513 行全部成功，无缺失 |
| 分箱 | DR：`−2.2…0.3`，宽 0.1（25 箱）；DL：`0…10`，宽 0.5（20 箱）。DL 的 0 峰包含 2,620 个精确 0 与 (0, 0.5) 的 493 个基因，合计 3,113 |
| 画法 | `cns.histplot(x=values, bins=edges, color=参数色, legend=False)`；DR 用砖红（槽 0）、DL 用青绿（槽 1），与 A/B 的参数配色一致；柱间无描边。两图都清掉 histplot 默认的 `Count` ylabel，只保留容器左侧的 `Number of genes` |
| 轴 | x `Depletion rate (DR)` / `Depletion lag (DL)`（7 pt）；y 共享 `0…3,362`（最大值 × 1.08），左图刻度 `0/1,000/2,000/3,000`（6.5 pt）；右图保留刻度线但不显示标签 |
| 产物 | 并入 `figure2.pdf` |

数据口径：DR 中位 −0.065、四分位 −0.637…0.001、范围 −2.091…0.168；DL 中位 0、75% 0.92、95% 4.89、最大 9.511。DL 的零膨胀是拟合下界 + L1 惩罚的结果，直方图第一根柱即该边界；DR 在 0 附近的双峰（[−0.1,0) 1,235、[0,0.1) 1,263）与 −1 附近的肩部是 DR 分布的真实形状。

---

### 可调旋钮（全部在 `workflow/scripts/figures/figure2.py` 顶部）

| 常量 | 当前值 | 作用 |
|---|---:|---|
| `PANEL_A_GENE` | `mtf1` | A 面板示例基因（固定） |
| `PANEL_B_GENES` | `tna1, pst2, rpl3801, nrs1` | B 面板四基因（行主序） |
| `PANEL_WIDTH_PT` / `PANEL_HEIGHT_PT` | 236 / 100 | A 面板容器尺寸 |
| `PANEL_B_WIDTH_PT` / `PANEL_B_HEIGHT_PT` | 236 / 258 | B 面板容器尺寸（含 58 pt 参数表带） |
| `FLOW_X_PT` / `FLOW_WIDTH_PT` | 2 / 78 | A 左侧流程列的位置与宽度 |
| `DATA_BOX_Y_PT` / `FORMULA_BOX_Y_PT` / `PARAMETER_BOX_Y_PT` | 80 / 40 / 14 | A 三个流程框的底边 |
| `PLOT_X_PT` / `PLOT_Y_PT` / `PLOT_WIDTH_PT` / `PLOT_HEIGHT_PT` | 126 / 22 / 100 / 76 | A 解剖图的 inset 位置与大小（B 的 mini axes 复用 100 × 76） |
| `B_LEFT_DECORATION_PT` / `B_COLUMN_GAP_PT` / `B_ROW_GAP_PT` | 25 / 6 / 14 | B 网格的左侧装饰、列距、行距 |
| `B_TOP_TITLE_BAND_PT` / `B_BOTTOM_DECORATION_PT` | 12 / 22 | B 顶部基因名带、底部 x 装饰 |
| `B_Y_PADDING_BOTTOM` / `B_Y_PADDING_TOP` | 0.08 / 0.12 | B 共享 y 轴上下留白 |
| `B_X_PADDING_LEFT` / `B_X_PADDING_RIGHT` | 0.09 / 0.03 | B 共享 x 轴左右留白（左多留是为了 DL = 0 的竖线） |
| `B_GUIDE_LINE_WIDTH_PT` / `B_GUIDE_DOT_SIZE` | 0.8 / 8 | B 每格 guide 的线宽与圆点大小（A 为 1 / 16） |
| `B_TABLE_BAND_PT` / `B_TABLE_BOTTOM_PT` | 58 / 2 | B 参数表带高、表格底边 |
| `B_TABLE_HEADER_HEIGHT_PT` / `B_TABLE_ROW_HEIGHT_PT` | 10 / 9 | 表头高、每行高 |
| `B_TABLE_GENE_X_PT` … `B_TABLE_R2_RIGHT_PT` | 42 / 111 / 144 / 177 / 214 | 表列锚点（Gene 左对齐，其余右对齐） |
| `B_TABLE_FONT_PT` / `B_TABLE_RULE_Y_PT` | 6.5 / 38 | 表字号、表头下横线高度 |
| `PANEL_C_WIDTH_PT` / `PANEL_C_HEIGHT_PT` | 252 / 100 | C 面板容器尺寸 |
| `C_LEFT_DECORATION_PT` / `C_COLUMN_GAP_PT` | 34 / 10 | C 的左侧装饰（y 刻度 + 共享 ylabel）与两图间距 |
| `C_PLOT_Y_PT` | 22 | C 两个直方图的底边（底部 22 pt 给 x 刻度 + xlabel） |
| `C_Y_LABEL_CENTER_X_PT` | 10 | C 共享 `Number of genes` 的 x 位置 |
| `C_Y_HEADROOM` | 1.08 | C 共享 y 轴上限 = 最大 bin × 1.08 |
| `DR_HIST_RANGE` / `DR_BIN_WIDTH` | −2.2…0.3 / 0.1 | DR 直方图范围与箱宽 |
| `DL_HIST_RANGE` / `DL_BIN_WIDTH` | 0…10 / 0.5 | DL 直方图范围与箱宽 |
| `ROW_GAP_PT` | 10 | A 行与 B 行之间 |
| `CHIP_ALPHA` | 0.25 | 公式参数底色透明度 |
| `Y_PADDING_BOTTOM` / `Y_PADDING_TOP` | 0.28 / 0.31 | A 面板 y 轴上下留白占数据跨度的比例；下限包含 A |
| `Y_TICK_STEP` | 2 | y 刻度间隔（只在数据范围内取刻度） |
| `X_PADDING` | 0.03 | x 两侧留白占跨度比例 |
| `FLOW_TEXT_FONT_PT` / `EQUATION_FONT_PT` | 6.5 | 流程文字与公式字号 |
| `PLOT_TICK_FONT_PT` / `PLOT_LABEL_FONT_PT` | 6.5 / 7 | 解剖图刻度与轴标题字号 |
| `EQUATION_LINE_GAP_PT` / `CHIP_PAD_PT` | 0.7 / 0.8 | 公式分段的水平间距、圆角底外扩 |
| `OBSERVED_MARKER_SIZE` / `GUIDE_DOT_SIZE` | 12 / 16 | 实测点与 guide 圆点大小 |

颜色、基因、公式分段与三条 guide 的标签都在同一批常量里；A 换基因只改 `PANEL_A_GENE`，轴范围与曲线跟着该行重算（但 guide 标签的绝对位置是按 mtf1 调的，换基因需复查）。B 换基因改 `PANEL_B_GENES`，共享 y 范围与网格跟着四行重算。

---

## Figure S1 — 文库制备重复性 + 每基因插入条数分布

**画布**：A4 整页 595.28 × 841.89 pt（210 × 297 mm）。内容画布宽 515.28 pt（= 595.28 − 2 × 40 pt 页边距）、高由两行面板推出（实测 304.4 pt）；组装时 pypdf 把内容页贴到 A4 页 (40, 40) pt 处（距左、距上、距右各 40 pt），其余留白。内容画布由 `cns.multipanel(max_width=515.28, unit="pt")` 流式排版：

```
A–C 行：7.2（字母带） + 120（面板高） + 30（x 刻度 + xlabel）
D 行：  10（行距）+ 7.2（字母带） + 100（面板高） + 30（x 刻度 + xlabel）
```

A–C 是三个 120 pt 方形面板并排（面板外框 = axes 120 + 实测 y 轴装饰 ≈29 + 右边距 10 ≈ 159 pt，三列共 477 pt，A4 宽度下右侧留白）；D 在第二行由 `mp.newline()` 从行首开始，用 house 的 SQUARE 尺寸 100 × 100 pt（单 panel 必须落在 80–150 pt 内，不能通栏），右侧留白。没有外部矢量面板，因此不需要示意图叠加；pypdf 只负责把内容页贴到 A4 页左上角内缩 40 pt 处。

**版式机制**：`mp.panel(width, height)` 只给 axes 尺寸；每列的横向宽度和字母位置在渲染时按实测的 y 轴装饰宽度算出，所以 A/B/C 的字母都落在各自列的外框左缘，A 与 D 的字母因此都在 x = 0（左对齐）。ylabel 变长只会把它自己这一列往右推，不会压到前一个面板。

**产物**：

| 文件 | 说明 |
|---|---|
| `results/7a_figures/figure_s1.pdf` | 最终矢量图：单页 A4，内容距左、上、右各 40 pt |
| `results/7a_figures/figure_s1.review.png` | 360 dpi 评审位图 |
| `results/7a_figures/_work/figure_s1_canvas.pdf` | 贴 A4 前的内容画布（515.28 × 304.4 pt）；排查版式时看这个 |

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
| 版式 | axes `w=120, h=120`，排到 `x=346.5`，右缘 466.5 pt（内容宽 515.28，行未占满）；字母 C 在列首，x ≈ 317.6 |
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
| `A4_WIDTH_PT` / `A4_HEIGHT_PT` | 595.28 / 841.89 | 最终页面尺寸（210 × 297 mm） |
| `PAGE_MARGIN_PT` | 40 | 内容页距 A4 页左、上、右的留白 |
| `CANVAS_WIDTH_PT` | 515.28 | 内容画布宽（= A4 宽 − 2 × 页边距），即 `multipanel(max_width=...)` |
| `PANEL_SIZE_PT` | 120 | A–C 每个方形面板的 axes 边长 |
| `BELOW_DECORATION_PT` | 30 | x 刻度 + xlabel 所占高度（`margin_bottom`） |
| `ROW_GAP_PT` | 10 | A–C 行与 D 行之间（D 的 `margin_top`） |
| `PANEL_D_SIZE_PT` | 100 | D 的 axes 边长（80–150 pt 契约） |
| `COUNT_BIN_WIDTH` / `COUNT_TAIL_START` | 25 / 300 | 直方图档宽与最后一档的起点 |
| `DENSITY_SCATTER_S` / `DENSITY_SCATTER_ALPHA` | 1.5 / 0.6 | 点大小与透明度（密云例外） |
| `DENSITY_CBAR_BOUNDS` | (0.82, 0.02, 0.035, 0.30) | 密度色条在 axes 分数坐标中的位置与大小 |

改 `PANEL_SIZE_PT` 后列宽和画布高度自动重算。三列总宽 = Σ(axes + 实测装饰 + `margin_right` 10)，当前 3 × 158.8 ≈ 476 pt，内容宽 515.28 pt 下还有富余；加到 ~130 pt 仍能并排，再大就会被 `max_width` 换行。面板间距用 `mp.panel(..., margin_right=...)` 调。

---

## Supplementary figures

（其余 S 图待加入：每个 S 图同样一条 rule、一节说明。）
