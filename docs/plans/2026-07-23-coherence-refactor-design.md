# complex.smk → coherence.smk 重构设计

**日期**: 2026-07-23
**worktree/分支**: `.worktrees/optimize-complex-smk` / `optimize-complex-smk`
**范围**: 拆分 `complex.smk`，从 clustering finalize-variant 解耦，改为按 **source** 扇出的
「整理分组数据 → 计算 coherence → 绘图」三段式，方便后续加入新的分组数据库（GO / KEGG / …）。
**前置背景**: 本设计把 [2026-07-19-followup-analysis-expansion.md](2026-07-19-followup-analysis-expansion.md)
§B [S0] + §C [D1] 的 notebook 探索方案落进 Snakemake DAG。attribution（§C [D2]）本次不接。

---

## 决策摘要（2026-07-23 头脑风暴）

经与用户逐项确认：

- **解耦输入**: coherence 分析不再挂 `final_clusters.tsv`（clustering finalize-variant 系统）。
  实测两个旧脚本只读 `Systematic ID / DR / DL`，**从不读 `cluster` 列**，而 final 表里的
  DR/DL 就是从 `release/gene_level/fitting_results.tsv` 原样带过来的。改为直接读
  `fitting_results.tsv`，结果不变，去掉一层多余耦合。
- **架构**: 分离 + source 扇出（镜像 clustering 的 `method` 扇出、enrichment 的 `ontology` 扇出）。
  加新库 = 加一个 source adapter + 一行 config，`compute`/`plot` 逻辑一行不改。
- **stage 重命名**: `complex` → `coherence`（`complex.smk`→`coherence.smk`、`results/complex/`→
  `results/coherence/{dataset}/{source}/`、`scripts/complex/`→`scripts/coherence/`）。理由:
  它现在跨 GO-CC 复合体 + GO-BP 过程两类分组，"complex" 名不副实；与设计文档命名一致。
- **本次实现 3 个 source / 2 个 loader**: `go_macrocomplex`（flat TSV，原逻辑）、`go_cc`、`go_bp`
  （后两者复用 enrichment 的 GAF + goatools 传播，按 namespace 参数化）。
- **module viz 改造**: 旧的硬编码五模块图 → 通用绘图模块：给定 term/complex（config 名单驱动），
  画其成员在 DR-DL 空间散点（叠背景云）+ term 信息 + coherence 结果标注。
- **attribution 本次不接**: `workflow/src/coherence/attribution.py`（major/minor GMM + 共享基因）
  留作下一个 PR；long table 的 `n_groups_per_gene` 已为它铺路。
- **不碰 byte-faithful 科学核心**: Weiszfeld geometric median、MPD permutation、DR/DL 归一化、
  `random_state=42`。`tests/test_complex_coherence.py` 的 24 个回归断言必须保持通过。
- **git mv + 改造**保留 blame 历史；`max_term_genes` 对所有 source 一致生效（含 go_macrocomplex）。

---

## §1 架构总览与数据流

```
release/gene_level/fitting_results.tsv  (DR/DL/A, dataset-specific)
   │        ┌─ (per SOURCE: go_macrocomplex / go_cc / go_bp) ───────┐
   │        │  PomBase 分组注释:                                      │
   │        │   · go_macrocomplex → macromolecular_complex_annotation.tsv (flat) │
   │        │   · go_cc  → GO GAF cellular_component NS + 传播 (part_of/is_a)     │
   │        │   · go_bp  → GO GAF biological_process NS + 传播                    │
   ▼        ▼                                                        │
[1] prepare_coherence_annotation  (per dataset × source)            │
      → results/coherence/{dataset}/{source}/group_annotation_long.tsv
        统一契约: source, group_id, group_name, Systematic ID, Name, n_group_genes
              │
              ▼
[2] compute_coherence  (per dataset × source, 确定性 seeded permutation)
      → results/coherence/{dataset}/{source}/coherence_metrics.tsv
      → results/coherence/{dataset}/{source}/coherence_analysis.pdf
              │
              ▼
[3] plot_group_scatter  (config 名单驱动的通用绘图)
      → results/coherence/{dataset}/{source}/group_scatter.pdf
```

**目录契约**: `results/coherence/{dataset}/{source}/...` — Snakemake 产出、可删可重建、git-ignored。
`{source}` 是新路径维度。`[2]/[3]` 对 source 完全无感，只吃 `[1]` 的统一契约 long table。

---

## §2 prepare 层 —— source-adapter 接口

**统一 long-table 契约**（所有 source 产出同一 schema）:

| 列 | 含义 |
|---|---|
| `source` | `go_macrocomplex` / `go_cc` / `go_bp` / …（冗余带上便于合并比较） |
| `group_id` | 分组稳定 ID（GO term ID，如 `GO:0005840`） |
| `group_name` | 人类可读名（`GO_term_name`，如 `cytosolic ribosome`） |
| `Systematic ID` | 成员基因 |
| `Name` | 基因 symbol；**缺失时回填该基因的 Systematic ID**（不留空） |
| `n_group_genes` | 该 term 传播后**全部**注释基因数（与 DR 无关；`max_term_genes` 过滤依据 + 为 D2 共享基因铺路） |

**source adapter 接口**（`workflow/src/coherence/sources.py`，新建）:

```python
def load_macrocomplex(pombase_dir: Path) -> pd.DataFrame: ...          # flat TSV，原 go_cc 逻辑
def load_gaf_namespace(pombase_dir: Path, namespace: str) -> pd.DataFrame:  # obo+gaf，按 NS 过滤+传播
SOURCE_LOADERS = {
    "go_macrocomplex": load_macrocomplex,
    "go_cc": lambda d: load_gaf_namespace(d, "CC"),
    "go_bp": lambda d: load_gaf_namespace(d, "BP"),
}   # 加库 = 加一项；GAF 三 namespace 共用一函数，加 go_mf 几乎零成本
```

- `load_gaf_namespace` 复用 `workflow/src/enrichment/ontology.py` 加载 obo + gaf 的 goatools 逻辑，
  按 namespace 过滤 term，用 `part_of`/`is_a` 传播展开成 term×gene 长表。
- prepare 脚本 `workflow/scripts/coherence/prepare_annotation.py`: `--source` + `--pombase-dir`
  → 查 `SOURCE_LOADERS[source]` → 写统一 long table。
- **env**: prepare 规则挂含 goatools 的 env（enrichment 那套 / `biopython.yml`），因为 go_cc/go_bp
  要用 goatools；go_macrocomplex 只是 flat TSV 改名，但同规则统一 env 无妨。

**三层过滤**（config 驱动，落在 compute 层比较，prepare 层备好 `n_group_genes`）:
- `min_group_size: 3` — DR>阈值成员下限（permutation 需要）。
- `max_group_size: 300` — DR>阈值成员上限（统计有效性 + go_macrocomplex byte-faithful `3≤n≤300`）。
- `max_term_genes: 500` — term 总注释基因数上限（砍 GO-BP/CC 宽泛父节点）。**对所有 source 一致生效**
  （含 go_macrocomplex，实际是 no-op 但代码路径统一）。500 为起步值，跑出 GO-BP 存活 term 数后可调。

---

## §3 compute 层 + 绘图

**`[2] compute_coherence`**（`compute_complex_coherence.py` → `compute_coherence.py`，git mv + 改造）:
- 输入: `--fitting-results`（release gene_level）+ `--annotation`（`[1]` long table）+ `--source`
  + （可选）`--features`（本 repo features 表，面板 B/D 用）。不再读 final_clusters、不再自解析 GO 关键词。
- **byte-faithful 科学核心一字不动**: 背景点云、Weiszfeld geometric median、seeded MPD permutation
  （`random_state=42`）。分组来源从"内部 groupby GO_term_name"改成"读 long table 的 group_id/group_name"。
- 三层过滤在此落地。输出 `coherence_metrics.tsv`: 新增 `source`/`group_id`/`n_group_genes`，
  `complex`→`group_name`；其余 schema 不变。

**`coherence_analysis.pdf`（确定性主图，全部进 PDF）**:
- 直方图: complex size 分布 + z-score 分布。
- 质心位置图: x=典型 DR, y=典型 DL（geometric-median centroid）, 色=z-score, 点大小=成员数。
  一眼看到每个模块在耗竭空间的位置 + 紧密度。
- **coherence×X 面板**（y 轴恒为 z-score，x 轴换独立生物学指标，解耦共线性）:
  - 面板 A｜× **共享亚基比例**: 从 long table `n_groups_per_gene` 算，零外部数据，**始终产出**。
    含义: 共享/moonlighting 成员多 → 质心被拉开 → 表观 incoherent。对应设计文档 D2(b)。
  - 面板 B｜× **蛋白丰度均匀度**（features 表 abundance CV）: `features_panels` 开关。
    含义: 化学计量主亚基丰度均匀=稳定核心→该 coherent。对应 D2(a)。
  - 面板 D｜× **保守性均匀度**（features 表 `evolutionary_rate` CV）: 同开关。
  - features 缺列 → 跳过该面板 + warn。features 是本 repo 确定性 Snakemake 产物，可安全作输入。
- **面板 C（内部 STRING/PPI 连接度）本次不做**（STRING 是 web API，非确定性，会污染确定性主链）。

**`[3] plot_group_scatter`**（新 `workflow/scripts/coherence/plot_group_scatter.py`）:
- 输入: `--fitting-results` + `--annotation` + `--metrics`（取 z/p/size）+ `--groups`（config 名单）。
- 每个指定 term/complex: 背景云（全 DR>阈值基因灰点）+ 成员高亮散点（复用
  `plotting/gene_level.py::plot_given_genes_on_feature_space`），标题/角注写
  `group_name (group_id)`、n_members、z-score、p-value。
- config `coherence.scatter_groups: {go_macrocomplex: [...], go_cc: [...], go_bp: [...]}`，
  按 source 分；名单项支持 group_name 或 group_id；缺某 source 跳过。

---

## §4 config / Snakefile / 测试 / 边界

**`analysis.yaml`: `complex:` → `coherence:`**

```yaml
coherence:
  sources: [go_macrocomplex, go_cc, go_bp]   # 扇出维度；加库在这里加一行
  min_group_size: 3          # DR>阈值成员下限 (permutation 需要)
  max_group_size: 300        # DR>阈值成员上限 (统计有效性 + go_macrocomplex byte-faithful)
  max_term_genes: 500        # term 总注释基因数上限 (砍宽泛父节点); 对所有 source 一致生效
  dr_threshold: 0.3
  n_permutations: 1000
  random_state: 42
  features_panels: true      # 开启 coherence×蛋白丰度/保守性面板 (读 features 表)
  scatter_groups:            # plot_group_scatter 手动名单, 按 source 分; 支持 group_name 或 group_id
    go_macrocomplex: [kinetochore, "mitochondrial large ribosomal subunit"]
    go_cc: []
    go_bp: []
```

**Snakefile**: `include: complex.smk` → `coherence.smk`；`rule all` 里注释掉的 `results/complex/...`
更新为 `results/coherence/{dataset}/{source}/coherence_metrics.tsv`。

**产物路径**（`{source}` 新维度）:
```
results/coherence/{dataset}/{source}/group_annotation_long.tsv
results/coherence/{dataset}/{source}/coherence_metrics.tsv
results/coherence/{dataset}/{source}/coherence_analysis.pdf      # 直方图+质心图+面板A/B/D
results/coherence/{dataset}/{source}/group_scatter.pdf           # config 名单驱动
```

**测试**（`tests/test_complex_coherence.py` → `test_coherence.py`，扩充）:
- 保留现有 24 个 byte-faithful 回归（geometric median、permutation z、known-coherent complex z<0）——不变。
- 新增: `sources.py` 每个 loader 契约测试（列=统一 schema、Name 缺失回填 Systematic ID、`n_group_genes` 正确）；
  `max_term_genes` 过滤逻辑；三层 size 过滤边界；go_bp GAF 传播用小型 fixture（不依赖真实全库）。

**落地边界**:
- 下沉脚本: `prepare_annotation.py` / `compute_coherence.py`（git mv 改造）/ `plot_group_scatter.py`。
- 本次不接 attribution（[D2] 下个 PR，`n_groups_per_gene` 已铺路）。
- 不碰 byte-faithful 科学核心。

**执行环境提醒**（见 [[complex-smk-optimization-worktree]]）:
- 测试跑在 `statistics_and_figure_plotting` env（需 `pip install pytest`）；go_bp loader 的 goatools
  需在含 goatools 的 env（enrichment/`biopython.yml`）里验证。
- 基线: 125 passed + 8 个与本 stage 无关的 collection error（缺 gffutils/Bio/requests/mljar），非回归。

