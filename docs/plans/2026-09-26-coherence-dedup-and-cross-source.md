# coherence：dedup 层的表/图/归因 + 跨 source 对比图

- 日期：2026-09-26
- 相关代码：`workflow/rules/3a_coherence.smk`、`workflow/scripts/coherence/`、
  `workflow/src/coherence/{attribution,palette}.py`、`Snakefile`

## 背景

`3a_coherence` 对 5 个 source 各自产出 `coherence.pdf` + `coherence_metrics.parquet`，
合并成 `coherence_metrics_combined.parquet`（q 在并集上重新 BH），再去冗余得到
`coherence_terms_deduplicated.tsv` / `coherence_terms_representatives.tsv`。
去冗余后的代表集合是整个 stage 的最终产物，但它此前只有 TSV：没有图，也没有
incoherence 归因；5 个 source 的 coherence 图也彼此独立，无法一眼横向比较。

同时发现一个潜在 bug：`group_id` **跨 source 不唯一**（实测 173 个，例如 `GO:0032040`
同时是 `go_cc` 的 term 和 `go_macrocomplex` 的 complex），而 attribution 的成员查找
只按 `group_id` 过滤。此前所有调用都是 per-source 所以没暴露，一旦把 pooled 表喂进去
就会静默取到两个 source 成员集的并集。

## 改动

### 1. dedup 层成为一等产物

| 新产物 | 规则 |
|---|---|
| `coherence_dedup.pdf` | `plot_coherence_dedup` |
| `dedup_incoherence_attribution.tsv` + `dedup_incoherence_split_points.parquet` | `compute_coherence_attribution_dedup` |
| `dedup_incoherence_attribution.pdf` | `plot_coherence_attribution_dedup` |

**表沿用已有的 `coherence_terms_representatives.tsv`**，不新建中间产物：它是
`coherence_terms_deduplicated.tsv` 的 `is_representative` 严格行子集，已经带
`source` / `group_id` / pooled `q_value` / `redundancy_cluster`。

q 值**原样带过**，不在代表集合上重新 BH。理由写在 `combine_metrics.py` 的 FDR 注释里：
代表是按「簇内 min q」选出来的，一族 per-cluster 最小值在原假设下不是均匀分布，
在代表集合上重跑 BH 是**反保守**而不只是冗余。

dedup 归因复用同一个 `compute_incoherence_attribution.py`，只是把 5 个 source 的
`group_annotation_long.tsv` 一起传进去（`--annotation` 改为 `nargs="+"`）：
per-source 调用传一张、pooled 调用传 5 张，和 `combine_metrics.py --metrics` 的写法一致。

### 2. `(source, group_id)` 键（根因修复）

`workflow/src/coherence/attribution.py`：

- `shared_subunit_fractions` 的返回键在表带 `source` 列时改为 `(source, group_id)`
  （原先只有 `group_id`，冲突的 source 会被后者静默覆盖）；不带 `source` 列时保持旧行为，
  与 `member_pairs` 自己的列判断保持一致。
- `shared_subunits(long_table, group_id, source=None)` 增加显式 `source`：给了就直接用，
  不给则维持「从 group_id 反推」的旧逻辑（单 source 表下等价，pooled 表下会取到任意一个）。

`compute_incoherence_attribution.py` 把 metrics 行的 `source` 一路透传，
split-points parquet 增加 `source` 列（pooled 图按 pair join 需要）。

**回归保证**：单 source 调用的输出逐字节不变（go_macrocomplex 的 attribution TSV
与改动前 `diff` 为空）；pooled 调用逐行等于对应的 per-source 行（只有 `q_value` 不同，
因为族不同）。这两条都有测试兜底。

### 3. 跨 source 对比图

`plot_coherence.py` 新增两个选项，`--color-by none`（默认）时行为**逐像素不变**：

- `--color-by source`：每个 panel 按 source 画一遍并按 source 上色。panel A/B 改为
  共享分箱的 step 轮廓（6 条实心柱会互相盖住），图例只挂在 A 上；panel C 让出 z 的连续配色
  （z 本来就是其余每个 panel 的 y 轴），去掉 colorbar，只留组大小图例。
- `--dedup-series PATH`：把代表集合以 `source="dedup"` 追加成第 6 个系列，
  和 5 个完整集合画在同一套坐标上，直接看出「去冗余后还剩什么」。

配色与顺序集中在 `workflow/src/coherence/palette.py`：固定顺序
`[go_macrocomplex, go_cc, go_bp, kegg_brite, kegg_pathway, dedup]`，颜色是手工挑的
Cell palette index，**不是 `house_colors(range(6))`** —— palette 的 0/1 两个位置
（红、青）灰度亮度只差 0.003，印出来是同一个灰。前三个 source 的颜色沿用
`plot_redundancy_network.py` 的既有分配，同一个 source 在两张图里同色。
`palette.py` 不放在 `sources.py`：后者被 compute 阶段 import，颜色 helper 会把
matplotlib 拖进纯计算环境；取色也延迟到调用时（`house_colors` 读的是
`apply_house_style()` 装好的 palette）。

### 4. pooled 归因图的读法

`plot_incoherence_attribution.py` 在三处按 `source` 区分（都以「表里有多于一个 source」
为条件，per-source 图完全不受影响）：panel 标题多一行 source、split-points 按
`(source, group_id)` join、label-frequency panel 改为按 source 堆叠。

top-N 仍是**全局 z 排序**，所以 16 个 panel 可能全部来自同一个 source ——
标题里的 source 让这件事是可见的，而不是被误解成「各 source 均匀分布」。

## 验证

```bash
PY=.snakemake/conda/<cnsplots 或 stats env>/bin/python
pytest -q --continue-on-collection-errors       # 与改动前失败集完全一致
pytest -q tests/test_coherence_figure_layout.py # 新增 color-by-source 布局与配色用例

# per-source 图必须逐像素不变
python - <<'EOF'
from PIL import Image; import numpy as np
a=np.asarray(Image.open('results/3a_coherence/HD_DIT_HAP/go_macrocomplex/coherence.review.png'))
b=np.asarray(Image.open('/tmp/coh_check/gm.review.png'))
assert (a==b).all()
EOF

snakemake --use-conda --cores 8     # 本机配方见 memory: snakemake-conda-env-recipe
```

## 已知遗留

- `attribution_top_n_plot` 对 pooled 表仍是全局排序，没有 per-source 配额。
  想保证每个 source 都露面需要加一个 `--top-n-per-source`。
- `plot_redundancy_network.py` 里的 `_SOURCE_ORDER/_SOURCE_COLORS` 只有 3 个 source，
  已经过期，没有随本轮合并到 `palette.py`（避免扩大改动面）。
