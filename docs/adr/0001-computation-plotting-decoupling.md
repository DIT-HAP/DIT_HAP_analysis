# ADR-0001: 计算与绘图解耦架构

**状态**: 已接受  
**日期**: 2026-09-02  
**决策者**: Yusheng Yang

---

## 背景

### 问题
本项目 (DIT_HAP_analysis) 当前使用单一的 `statistics_and_figure_plotting.yml` 环境处理计算和绘图，包含：
- 高版本计算库：pandas≥3.0.3, numpy≥2.5.0, matplotlib≥3.11
- 绘图需求：54 个 scripts 中 11 个涉及绘图，其中 4 个是计算+绘图混合脚本

### 约束
上游项目 DIT_HAP_snakemake 已实现计算与绘图解耦，使用 cnsplots 0.7.0 包提升绘图质量。但 cnsplots 0.7.0 要求：
- `matplotlib>=3.10,<3.11`（与计算环境的 3.11+ 冲突）
- `pandas>=2.2,<2.4`（与计算环境的 3.x 冲突）

### 目标
1. 对齐上游项目的架构模式（完全对齐）
2. 引入 cnsplots 包提升绘图质量
3. 避免依赖版本冲突
4. 支持计算任务在无图形环境运行

---

## 决策

### 核心架构：脚本级解耦 + 环境级隔离

#### 1. 环境分离
创建独立的绘图环境 `workflow/envs/cnsplots.yml`：
```yaml
dependencies:
  - python=3.12
  - pandas>=2.2,<2.4          # 显式固定，降级以兼容 cnsplots
  - matplotlib>=3.10,<3.11    # 显式固定
  - numpy>=2.4,<2.5
  - seaborn>=0.13.2,<0.14
  - pyarrow                   # Parquet 支持
  - pytest>=9.0
  - loguru
  - pip:
    - cnsplots==0.7.0
```

计算环境 `statistics_and_figure_plotting.yml` 保持高版本：
- pandas≥3.0.3, matplotlib≥3.11

#### 2. 脚本拆分
将 4 个混合脚本（`coherence/compute_coherence.py` 等）拆分为：
- **计算脚本**：`compute_<功能>.py` → 输出 Parquet 数据
- **绘图脚本**：`plot_<功能>.py` → 读取 Parquet，输出 PDF

**示例**（coherence 模块）：
```python
# 计算规则
rule compute_coherence:
    conda: "../envs/statistics_and_figure_plotting.yml"
    output: "results/coherence/{dataset}/coherence.parquet"
    script: "../scripts/coherence/compute_coherence.py"

# 绘图规则
rule plot_coherence:
    conda: "../envs/cnsplots.yml"
    input: "results/coherence/{dataset}/coherence.parquet"
    output: "results/coherence/{dataset}/coherence.pdf"
    script: "../scripts/coherence/plot_coherence.py"
```

#### 3. 样式系统迁移
- **废弃**：`config/DIT_HAP.mplstyle`（403 行 matplotlib 样式文件）
- **采用**：`workflow/src/figures.py`（复制上游 409 行，通过 Python 代码控制样式）
  - `apply_house_style()`：调用 `cnsplots.setup_matplotlib()` + 自定义参数
  - `save_dual(output_stem)`：保存 PDF + PNG 双输出

#### 4. 中间数据格式约定
- **优先**：Parquet（精度无损，跨语言可读，pandas 原生支持）
- **必要时**：Pickle（scikit-learn 模型对象、复杂 Python 数据结构）
- **最终输出**：TSV（人类可读，供审查）

#### 5. 规则组织
创建 `workflow/rules/figure.smk`，集中管理所有绘图规则：
- 将 11 个绘图相关规则从各领域 .smk 移至 `figure.smk`
- 清晰的职责分离：计算规则（生成数据）vs 绘图规则（渲染图表）

---

## 命名约定

### 规则命名（前缀模式，对齐上游）
- 计算规则：`compute_<功能>`（如 `compute_coherence`）
- 绘图规则：`plot_<功能>`（如 `plot_coherence`）
- 数据生成：功能性名词短语（如 `coherence_metrics`）

### 产物路径（不加类型后缀）
- `results/coherence/{dataset}/coherence.parquet`（不是 `coherence_metrics.parquet`）
- `results/coherence/{dataset}/coherence.pdf`（不是 `coherence_figure.pdf`）

---

## 实施策略

### Phase 1: 基础设施
1. 创建 `workflow/envs/cnsplots.yml`
2. 复制上游 `workflow/src/figures.py`（409 行）
3. 创建空的 `workflow/rules/figure.smk`

### Phase 2: 纯绘图脚本迁移（验证环境）
**严格顺序执行**：先完成所有 7 个纯绘图脚本迁移，验证环境可用性，再进入 Phase 3。

迁移步骤（每个脚本）：
1. 修改导入：`from workflow.src.plotting` → `from workflow.src import figures`
2. 修改 Snakemake 规则的 `conda:` 指向 `cnsplots.yml`
3. 将规则从原 `.smk` 移动到 `figure.smk`
4. 生成新图表，人工目视对比旧图（验证质量）
5. 在旧脚本文件名中标记为 deprecated

### Phase 3: 样式调整
- 对比新旧图表的视觉效果
- 如需要，调整 `figures.py` 中的样式参数（颜色、字体等）

### Phase 4: 混合脚本拆分
拆分 4 个混合脚本（`coherence/`, 等），生成对应的 `compute_` 和 `plot_` 规则。

### Phase 5: 清理
- 标记 `workflow/src/plotting/` 为 deprecated（在文件名中标记）
- 删除 `config/DIT_HAP.mplstyle`
- 迁移完成后删除旧模块

---

## 后果

### 优点
1. **依赖隔离**：计算环境和绘图环境的版本冲突完全避免
2. **架构对齐**：与上游项目保持一致，便于代码复用和维护
3. **绘图质量**：cnsplots 包提供更专业的图表样式和 API
4. **环境最小化**：计算任务可在无图形环境运行（如 HPC 计算节点）
5. **可读性提升**：Parquet 格式比 pickle 更通用，支持跨语言工具链
6. **职责清晰**：计算和绘图的边界明确，规则组织更清晰

### 缺点
1. **环境数量增加**：从 3 个 conda 环境增至 4 个（增加维护负担）
2. **依赖冗余**：两个环境都包含 pandas/matplotlib，但版本不同
3. **迁移成本**：11 个脚本需要改写，4 个混合脚本需要拆分
4. **规则数量增加**：拆分混合脚本后，Snakemake 规则数量增加
5. **学习曲线**：团队需要学习 cnsplots API（替代直接使用 matplotlib）

### 风险与缓解
| 风险 | 影响 | 缓解措施 |
|------|------|----------|
| cnsplots 包停止维护 | 长期可维护性 | 保留 `workflow/src/plotting/` 的核心逻辑（特别是 `gene_level.py`），确保可以回退到纯 matplotlib |
| 图表样式与现有不一致 | 论文图表需要重新审查 | Phase 3 提供样式调整窗口，必要时在 `figures.py` 中覆盖 cnsplots 默认设置 |
| Parquet 与旧代码不兼容 | 现有分析脚本可能失败 | 渐进式迁移，保留旧 pickle 文件直到验证完成 |
| 环境构建时间增加 | CI/CD 性能下降 | 使用 conda/mamba 缓存，固定版本确保可复现 |

---

## 替代方案（已拒绝）

### 方案 A: 环境级解耦（保持脚本不变）
**描述**：创建独立绘图环境，但不拆分混合脚本。同一个脚本根据需要在不同环境运行。

**拒绝理由**：
- 混合脚本仍然包含计算和绘图逻辑，无法在无图形环境运行计算部分
- 与上游项目的脚本级分离架构不一致

### 方案 B: 不引入 cnsplots，保持现有 plotting 模块
**描述**：只创建环境隔离，但继续使用 `workflow/src/plotting/` 模块，不引入 cnsplots。

**拒绝理由**：
- 无法获得 cnsplots 的绘图质量提升
- 与上游项目的 API 不对齐，日后合并代码困难

### 方案 C: 升级 cnsplots 版本，避免依赖冲突
**描述**：等待或推动 cnsplots 包支持 matplotlib 3.11+。

**拒绝理由**：
- cnsplots 不是本项目可控的外部依赖，升级时间不确定
- 即使升级，pandas 2.x → 3.x 仍有 breaking changes，兼容性不保证

---

## 参考资料

- 上游项目：`../DIT_HAP_snakemake/`
  - `workflow/envs/cnsplots.yml`
  - `workflow/src/figures.py`
  - `workflow/rules/figure.smk`
- cnsplots 文档：https://pypi.org/project/cnsplots/
- 本项目现状调查：见 grilling session 2026-09-02

---

## 审查记录

| 日期 | 审查者 | 决定 | 备注 |
|------|--------|------|------|
| 2026-09-02 | Yusheng Yang | 接受 | 通过 grilling session 确认所有设计细节 |
