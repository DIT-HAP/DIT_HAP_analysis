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

---

# 追加（2026-09-28）：kegg_brite 跟随上游新 schema

上游 kegg_parser 现在把每一层的标签和 id 分成两列（`Level_X` / `Level_X_ID`），
并且 **id 永远非空**。`sources.py::load_kegg_brite` 相应重写（`load_kegg_pathway` 不动，
`pathway_gene_mapping.tsv` schema 没变）。

## 关键：`Level_D_ID` 非空 ≠ `Level_D_ID` 是 id

kegg_parser 自己的契约（`src/kegg_parser/brite.py:16-21`）写着：

> bracket ids are kept bare (`[PATH:spo00010]` -> `spo00010`, joinable with
> `Pathway_ID`), leading 5-digit codes become `map:00566` / `class:09100`,
> leading EC numbers become `EC:1.1.1.1`, **and labels carrying no id at all fall
> back to the label text itself so ids stay non-empty and text-safe**.

也就是说，KEGG 没给 id 的节点，`Level_D_ID` 就是**标签本身**。而标签恰恰是不唯一的那个东西
（核糖体树里四条分支都以 `Large subunit` 结尾）。所以「非空」不等于「可用作节点主键」。

实测（13,141 行）：

| | 行数 | 跨 tree 碰撞 | 一个 id 对多个节点 |
|---|---|---|---|
| `Level_D_ID != Level_D`（真 id） | 8,642 | 0 | 0 |
| `Level_D_ID == Level_D`（标签回退） | 4,499 | 17 | 68 |

两个 regime 分得干干净净，所以 `Level_D_ID == Level_D` 就是「没有真 id」的判据。

## 新的 group_id 规则（`_brite_group_ids`）

```
有真 id（Level_D_ID != Level_D） →  直接用该 id      （spo00010 / EC:1.1.1.1 / map:00566）
否则                             →  BRITE_ID + ":" + Level_A > Level_B > Level_C > Level_D
```

兜底必须是**完整路径**而不是叶子标签：`Level_D` 在一棵树内就不唯一，按它建 id 会把
核糖体树的四条分支并成一个组。实测 `label_path -> id_path` 是单射、两者都是 1,976 个节点，
而裸 `Level_D_ID` 只有 1,843 个、`BRITE_ID:Level_D_ID` 只有 1,869 个 —— 都少了。

## 效果

| | 旧 | 新 |
|---|---|---|
| kegg_brite 组数（annotation） | 1,869 | 1,937 |
| kegg_brite 组数（过 metrics 过滤） | 369 | 388 |
| 其中真 id / 路径回退 | — | 147 / 241 |
| `Large subunit` | 1 个组（119 基因，四分支并集） | 4 个组（87 / 50 / 32 / 21） |
| combined | 2,568 | 2,587 |
| dedup 代表 | 623 | 626 |

另外 4 个 source 的 per-source 表逐字节不变（`equals()` 验证）。

## 遗留

- **`Level_D` 名字仍会重复**：拆开后 133 个组的 `Level_D` 与别的组重名（`Large subunit` 之类）。
  本轮没动 `group_name`（仍是 `Level_D`），因为改它会把所有 kegg_brite 图/表的标签换一遍。
  要区分可以改成路径里最后两个不同的层（如 `Eukaryotes > Large subunit`）。
- 真 id 里有 130 个与 `kegg_pathway` 的 `Pathway_ID` 重叠（`spo00010` 这类）—— 那是同一批
  KEGG pathway 的两种视图，pooled 表靠 `(source, group_id)` 区分，dedup 也会把它们并到一起，
  这是正确的。
- dedup 的 single-linkage 传递闭包问题见另一节。

---

# dedup 的 single-linkage 链式合并

## 机制

`deduplicate_terms.py::build_clusters` 做的是**单链（single-linkage）连通分量**：
两个 term 之间的边 = 成员集 Jaccard ≥ `dedup_jaccard_threshold`（config 里 0.5）。
连通分量具有传递性 —— A~B 且 B~C 会把 A、B、C 放进同一簇，**即使 A~C 相似度为 0**。
所以一个簇的紧密度只由它最弱的那条链决定，不由阈值决定。

`merge_dag_lineage` 默认关（config 注释里记了原因：开着会让 2,046/2,097 个 GO term
塌成一个簇），所以这里没有第二条边规则。

## 证据（2,587 个 term → 626 个簇）

| 指标 | 值 |
|---|---|
| 簇内 pair 总数 | 36,999 |
| 其中成员**完全不重叠**（J=0） | 13,888（**37.5%**） |
| ≥5 个 term 的簇（121 个）里最弱 pair 的 Jaccard | 中位数 0.25，69 个 < 0.3，18 个 = 0 |
| 最大簇 | **193 个 term，并集只有 191 个基因** |

**最大簇 `all:17`**：193 个 term（go_bp 166 / kegg_brite 11 / kegg_pathway 10 / go_cc 5 /
go_macrocomplex 1），并集 191 个基因，代表是 `GPI anchored protein biosynthesis`。
也就是说表里看到 1 行，另外 192 个 term 被它代表 —— 而这些 term 之间大多是
`organophosphate metabolic process` / `lipid biosynthetic process` 这种互相重叠的宽泛 BP，
度数只有 3~6，靠弱边连成一片。

**核糖体那个簇**（21 个 term，含 `mitochondrial matrix` 133 基因）：胞质和线粒体大亚基
被连在了一起。链是

```
kegg Mito-LSU  --J=0.73-->  kegg Bacteria-LSU
kegg Mito-LSU  --J=0.51-->  cytosolic ribosome
cytosolic LSU  --J=0.52-->  cytosolic ribosome
cytosolic LSU  --J=0.97-->  kegg Eukaryotes-LSU
```

55 条边就把 21 个 term 连成一簇，代表是 `mitochondrial ribosome` ——
于是 `cytosolic large ribosomal subunit`（一个真实的、和线粒体核糖体完全不同的东西）
在 dedup 视图里消失了。

## 阈值敏感性

| threshold | 簇数 | 最大簇 | 非代表 | 单例 |
|---|---|---|---|---|
| 0.5（当前） | 626 | **193** | 1,961 | 252 |
| 0.6 | 920 | 46 | 1,667 | 401 |
| 0.7 | 1,125 | 31 | 1,462 | 528 |
| 0.8 | 1,346 | 27 | 1,241 | 714 |
| 0.9 | 1,579 | 23 | 1,008 | 971 |

193 → 46 只要把阈值从 0.5 提到 0.6：那个巨簇是靠 0.5~0.6 区间的边撑起来的，典型链式特征。

## 可选的处理方向（未实施）

1. **提高阈值** —— 一行 config，0.6 就把最大簇从 193 砍到 46。代价是留下更多冗余。
2. **换连接准则** —— 平均链（average-linkage）或全连接（complete-linkage，簇内所有 pair
   都 ≥ 阈值）。能根治链式，但要改 `build_clusters`，且是 O(n²)。
3. **限簇大小** —— 超过 N 就不再合并，代价是结果依赖合并顺序。
4. **不动算法，只让读者知道** —— `cluster_size` 列已经在表里，`redundancy_network.html`
   就是用来逐个簇审查的。193 个 term 的簇在表里是可见的，只是容易被忽略。

## 这不是"引入的 bug"

是单链聚类的定义本身。repo 里 config 已经把 DAG lineage 关掉（同一类问题的另一种表现），
说明这个折中是明知的。列在这里是因为它的规模（1,961/2,587 个 term 是非代表）值得知道。
