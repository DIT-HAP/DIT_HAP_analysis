# CONTEXT.md

## 项目术语表

### 计算与绘图解耦（Computation-Plotting Decoupling）
架构模式，将数据计算和图表绘制分离为独立的脚本和环境。

- **计算脚本**：负责统计分析、聚类、富集计算等，输出 Parquet/TSV 格式的结构化数据。运行在 `statistics_and_figure_plotting.yml` 环境（pandas 3.x, matplotlib 3.11+）。
- **绘图脚本**：负责读取计算结果并渲染图表，输出 PDF/PNG。运行在 `cnsplots.yml` 环境（pandas 2.x, matplotlib 3.10, cnsplots 0.7.0）。

**动机**：cnsplots 包要求 matplotlib<3.11，与计算环境的新版本依赖冲突。通过环境隔离避免依赖冲突。

**命名约定**：
- 规则名：`compute_<功能>` / `plot_<功能>`（前缀模式，对齐上游项目）
- 产物路径：直接描述内容，不加类型后缀（如 `coherence.parquet`, `coherence.pdf`）

---

### 中间数据格式（Intermediate Data Format）
Snakemake 规则之间传递的数据格式。

- **优先格式**：**Parquet**（精度无损，跨语言可读，节省空间）
- **必要时使用**：Pickle（scikit-learn 模型对象、复杂 Python 数据结构）
- **人类可读输出**：TSV（最终的、供人工审查的结果表）

**历史背景**：项目曾约定"中间文件用 pickle"（见 ADR-XXXX），现改为"Parquet 优先"以提升可读性和工具兼容性。

---

### 样式系统（Style System）
控制 matplotlib 图表视觉风格的机制。

- **当前实现**：`workflow/src/figures.py`（Python 代码，调用 `cnsplots.setup_matplotlib()` + 自定义参数）
- **已废弃**：`config/DIT_HAP.mplstyle`（403 行 matplotlib 样式文件，将在迁移完成后删除）

**核心函数**：
- `apply_house_style()`：设置全局调色板和布局参数
- `save_dual(output_stem)`：保存 PDF（期刊质量）+ PNG（评审预览）

---

### cnsplots
外部 PyPI 包（版本 0.7.0），提供 50+ 高级绘图函数（箱线图、热图、火山图等）。

- **环境**：`workflow/envs/cnsplots.yml`
- **依赖约束**：pandas 2.2-2.4, matplotlib 3.10-3.11（与计算环境隔离）
- **使用场景**：所有绘图脚本统一通过 cnsplots API 生成图表

---

### figure.smk
Snakemake 规则文件，集中管理所有绘图规则。

- **职责**：包含所有输出 PDF/PNG 的规则（`plot_variant_clusters`, `plot_coherence` 等）
- **对比**：计算规则（生成数据）保留在各自的领域 .smk 文件中（`clustering.smk`, `enrichment.smk` 等）

**动机**：清晰的职责分离，便于统一管理绘图环境的切换。

---

### 上游项目（Upstream Project）
指 DIT_HAP_snakemake 项目（`../DIT_HAP_snakemake/`），本项目的数据来源和架构参考。

- **关系**：本项目消费上游的 `release/` 输出，进行下游分析（聚类、富集、ML）
- **对齐点**：
  - 解耦架构（脚本级分离 + 环境隔离）
  - cnsplots 包的使用（版本 0.7.0）
  - 样式系统（`figures.py` 的 API 设计）
  - 命名约定（`compute_`/`plot_` 前缀）
