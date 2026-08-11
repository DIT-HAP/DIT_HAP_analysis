# Pombe 基因注释工具（给任意表格加注释列）

**日期：** 2026-08-11
**范围：** 中等。新增一个通用 CLI：给定任意含 pombe 系统名列的表格，追加 24 列注释
（budding yeast 同源基因及其必需性、pombe 自身必需性、功能注释、gRNA 层 DR/DL），
帮助判断基因功能与必需性。注释数据源由一条 Snakemake 规则预先解析成一张参考表缓存，
CLI 只做 join。

**实现状态：** 已完成。`7fea7e9`（注释工具本体，22 列）+ `b2d4bce`（gRNA 层 DR/DL，
→ 24 列）。66 个测试通过；参考表在真实数据上构建为 12,685 基因 × 24 列。下文标注
「实测」的数字来自实现期间的真实数据验证。

## 1. 决策与取舍（用户已定）

| 决策点 | 选定方案 | 理由 |
|---|---|---|
| 使用形态 | 通用 CLI + 预构建参考表 | 既能进分析流程，也能随手注释手头任何表格 |
| SGD 数据获取 | 手动 fetch 脚本，不进 DAG | 与 `resources/external/` 现状一致（pombase/biogrid 全是手工填的，git-ignored） |
| 注释范围 | 同源基因块 + pombe 必需性 + 功能注释 | 定量背景（表达量/网络/paralog）暂不做（YAGNI），已有 features 矩阵可查 |
| 多同源处理 | 每个同源基因单独给值，同序拼接 | 「其中一个必需」与「两个都必需」含义不同，汇总会丢信息 |
| ID 解析 | **不过 `update_sysIDs()`**，直接按系统名 join | 用户表格用的都是最新 ID |
| japonicus 同源 | 不纳入 | 主要用于保守性打分，对判断功能帮助不大 |
| gRNA DR/DL | 只加 `gRNA_DR`/`gRNA_DL` 两列，**不加**拟合参数 | 见 §7；命名加前缀以区别于基因层 DR/DL |

**为什么拆成「参考表 + join 脚本」而不是一个脚本全干：** 注释源解析要读 GAF/OBO
（GO-slim 需要 `GODag`）、20 万行的 SGD phenotype 表、复合体表。这些开销与输入表格无关，
解析一次缓存成 parquet，之后每次注释就是一个 `merge`。也符合仓库「per-stage intermediates
是 parquet、只有最终产物是 TSV」的约定。

## 2. 架构与数据流

```
resources/external/sgd/{sgd_version}/        ← 手动跑一次 fetch_sgd_data.sh
   ├── SGD_features.tab                        (16461 行)
   └── phenotype_data.tab                      (200084 行)
          │
          ▼  规则 build_annotation_reference（按 pombase_version × sgd_version 构建一次）
results/annotation/{pombase_version}/{sgd_version}/gene_annotation_reference.parquet
   宽表，index = pombe systematic id，24 列注释（实测 12,685 基因 × 24 列）
          │
          ▼  通用 CLI（手动跑，不进 DAG）
annotate_pombe_genes.py --input any.tsv --gene-column X --output annotated.tsv
```

另外两个输入不经 SGD，直接来自 `resources/curated/`：
`deletion_library_categories.xlsx`（pombe 缺失文库必需性）和
`260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv`（gRNA 层 DR/DL）。

新增文件：

- `workflow/scripts/annotate/fetch_sgd_data.sh` — `curl` 两个 SGD 文件到
  `resources/external/sgd/{date}/`。手动跑一次，把「怎么拉的」记录下来。下载到 `.part`
  再改名，避免中断留下看起来完整的截断表。**刻意不用 `--continue-at`**：文件已完整时
  SGD 对越过 EOF 的 range 请求返回 HTTP 416，会让每次重跑都失败（实现期实测踩到）。
- `workflow/src/annotation/core.py` — 纯函数库：每个注释块一个 `build_*` 函数，
  加 `assemble_annotation_reference()` 组装、`annotate_table()` 做 join、
  `summarise_match()` 报未匹配。
- `workflow/scripts/annotate/build_annotation_reference.py` — 构建参考表的 CLI。
- `workflow/scripts/annotate/annotate_pombe_genes.py` — 注释 CLI。
- `workflow/rules/annotate.smk` — 只放构建参考表那条规则（`Snakefile` 里 `include:`）。
- `tests/test_annotation.py` — 见 §7。

## 3. 同源基因块（10 列）

`pombe_cerevisiae_orthologs.txt` 的格式（5113 行，pombe id 唯一，1103 个 `NONE`，
795 行含 `|`）有三种形态，解析必须都覆盖：

```
SPBC1711.13    YCL030C                       单一同源
SPAC1002.07c   YEL066W|YPR193C               多个独立同源，| 分隔
SPCC1450.15    YBR265W(N)+YDR302W(C)         融合基因：pombe 一个基因对应 Sc 两个，
                                             (N)/(C) 表示只匹配 N 端 / C 端
SPBC29A3.02b   NONE                          无同源
```

`(N)`/`(C)` 保留不丢——它是真实生物学信息（该 pombe 基因是 Sc 两个基因的融合，功能可能
只对应一半）。`+` 表示同一同源关系的多个片段，`|` 表示互相独立的多个同源基因。两者都拆开
去查名字与必需性，但输出保留原始字符串。

| 列名 | 内容 | 来源 |
|---|---|---|
| `Sc_ortholog_id` | `YCL030C`，多个用 `\|` 拼 | PomBase curated_orthologs |
| `Sc_ortholog_name` | `HIS4`，与上列同序；无常用名回落到系统名 | SGD_features.tab |
| `Sc_ortholog_count` | 独立同源基因数（`\|` 段数），无同源为 0 | 计算 |
| `Sc_ortholog_raw` | 原样字符串，含 `(N)`/`(C)`/`+` | PomBase |
| `Sc_ortholog_qualifier` | `Verified` / `Dubious` / `Uncharacterized`，同序拼接 | SGD_features.tab |
| `Sc_essentiality` | `inviable` / `viable` / `conflicting` / 空，同序拼接 | SGD phenotype_data |
| `Sc_essentiality_evidence` | `inviable:3\|viable:1`（文献条目数） | SGD phenotype_data |
| `Sc_description` | SGD 基因功能描述 | SGD_features.tab (col16) |
| `Hs_ortholog_symbol` | `KDM5A\|KDM5B\|...` | PomBase pombe_human_orthologs |
| `Hs_ortholog_count` | 同上计数 | 计算 |

**基因名回落是常态，不是边缘情况。** SGD 6613 个 ORF 里只有 5312 个有常用名，约 1300 个
只有系统名。

**qualifier 列的必要性。** ORF 分三档：Verified 5783 / Dubious 683 / Uncharacterized 147。
若某 pombe 基因的「同源基因」落在 Dubious ORF 上，这条同源关系的证据强度很弱——SGD 认为
那个 ORF 可能根本不编码蛋白。这直接影响如何看待注释结果。

> **实测后记：这列目前几乎是常数。** PomBase 的 curated ortholog 实际上只指向 Verified
> ORF：被引用的 4219 个 ORF 里 4218 个是 Verified，剩下 1 个（`YER109C`）在 SGD 的 ORF
> 表里已不存在（已退役注释）。保留此列的价值在于将来 PomBase 版本若指向 Dubious ORF 能
> 报警，但今天它信息量很低。

**实测：对齐不变量成立。** 全部 5113 行中，`Sc_ortholog_name` /
`Sc_ortholog_qualifier` / `Sc_essentiality` 的分段数与 `Sc_ortholog_id` 完全一致
（0 处错位）。融合基因如 `SPCC1450.15` → `YBR265W+YDR302W` / `TSC10+GPI11` 正确保持为
一组。同源数分布：0 个的 1103（与 `NONE` 计数吻合）、最多 17 个。

**`conflicting` 的来历（已核实）。** SGD phenotype 表里 `null` 突变体：1245 个基因标
inviable、5223 个标 viable，**其中 202 个基因两种标注都有**（不同文献/背景结论不一致）。
不强行二分，输出 `conflicting` 并由 evidence 列给出各自条目数，让使用者看到证据强度，
而不是被一个可能错的二分标签误导。

判定规则：取 `mutant_type == "null"` 且 `phenotype` 属于 `{inviable, viable}` 的行，
按基因分组计数 → 只有 inviable → `inviable`；只有 viable → `viable`；两者都有 →
`conflicting`；无记录 → 空。只认精确的 `viable`/`inviable`——`viability: decreased`
（321 行）、`viability: increased`（42 行）是分级表型，不是 null 突变体存活判定。

**实测输出：** 6266 个 ORF 有存活判定 —— viable 5021 / inviable 1043 / conflicting 202，
与独立的 awk 基线逐个吻合。典型的 `conflicting`：`YAL035W` 是 `inviable:2|viable:2`
（若强行二分就是掷硬币），`cdc42` 的同源基因是 `inviable:14`（证据很强）。

**实测踩到的坑：SGD 的 phenotype 文件是残缺的。** 20 万行里 39 行多出第 15 个字段、
19 行被内嵌换行截成「13 字段行 + 短续行」。`pd.read_csv` 对此直接抛
`ParserError` 拒绝整个文件（`usecols` 也不行，短行会触发 out-of-bounds）。因此
`read_sgd_phenotype_data()` 自己按行切分并按已知的 14 列补齐/截断。

## 4. pombe 侧列（7 列）

| 列名 | 内容 | 覆盖 |
|---|---|---|
| `gene_name` | `mrx11`；无名回落到系统名 | 全部 |
| `gene_product` | PomBase product 描述 | 全部 |
| `synonyms` | 旧 ID / 别名 | 部分 |
| `Sp_FYPO_viability` | `viable` / `inviable` / `condition-dependent` / `unknown` | 12685 基因，但 7695 是 `unknown` |
| `Sp_deletion_essentiality` | `E` / `V`（Hayles 缺失文库判定） | 4843 基因（E 1267 / V 3576） |
| `Sp_deletion_phenotype` | 如 `misshapen essential` | 4843 |
| `Sp_growth_category` | 如 `microcolonies` / `WT-like` | 4843 |

**两个必需性来源都保留，因为不是一回事：** `gene_viability.tsv` 是 PomBase 汇总的 FYPO
注释（覆盖全基因组但 61% 是 `unknown`）；deletion library 是 Hayles 系统性缺失文库的实测
（只覆盖 4843 个，但每个都有判定）。一起看比压成一列可靠。

来源列名（`resources/curated/deletion_library_categories.xlsx`，4843×24）：
`Gene dispensability. This study` → `Sp_deletion_essentiality`；
`Phenotypic classification used for analysis` → `Sp_deletion_phenotype`；
`Category` → `Sp_growth_category`。

**实测：实际覆盖 4842 而非 4843。** `SPBC8E4.02c` 在 curated 缺失文库里，但已不在
PomBase 2026-06-01 的基因集中（已退役 ID）。这是设计意图——PomBase 基因集定义参考表的
行集合，外部表格多出来的基因不会被引入。

## 5. 功能注释（4 列）

| 列名 | 内容 | 覆盖 |
|---|---|---|
| `GO_slim_BP` | 生物过程 GO-slim 条目名，`\|` 拼接 | 复用 enrichment 的 `load_ontology_data` |
| `GO_slim_CC` | 细胞组分 | 同上 |
| `GO_slim_MF` | 分子功能 | 同上 |
| `complex` | 所属大分子复合体名 | 1957 基因，最多 12 个 |
| `PFAM_domains` | PFAM ID，`\|` 拼接 | 4577 基因 |

GO-slim 走已有的 `workflow/src/enrichment/ontology.py`——`OntologyDataConfig` 已经知道
怎么把三个 slim 表 + OBO + GAF 加载成 `slim_dag`，复用它把基因映射到 slim 条目，不重写。
用 slim 而非完整 GO：完整 GO 一个基因动辄十几条，读不了；slim 是刻意设计的粗粒度可读标签。

`PFAM_domains` 只给 ID 不给名字——PomBase 的 `protein_families_and_domains.tsv` 没有
domain 名称列，要名字得再接一个 Pfam 数据源。暂不加（YAGNI）。

**实测覆盖：** GO-slim BP 4663 / CC 5182 / MF 4348，complex 1957，PFAM 4577。

## 6. gRNA 层 DR/DL（2 列）

| 列名 | 内容 | 覆盖 |
|---|---|---|
| `gRNA_DR` | gRNA 层最大耗竭速率（源表列名 `um`） | 5050 基因 |
| `gRNA_DL` | gRNA 层滞后（源表列名 `lam`） | 5050 基因 |

**为什么加 `gRNA_` 前缀而不直接叫 `DR`/`DL`——这是本节的核心。** 它们与 release 里的
基因层 `DR`/`DL` **不是同一个数**：curated 表是每个基因一条代表性 gRNA 的拟合，release 是
基因层聚合拟合。在 4465 个共有基因上实测：DR 相关 0.92，但 **DL 只有 0.55**，3824 个基因
差异 > 0.01（只有 53 个完全相同）。直接叫 `DR`/`DL` 会在 join 聚类表时撞名（被加
`_annotation` 后缀），更糟的是会暗示两者可以互换。

实际对比（注释 `gmm_direct9/final_clusters.tsv`）：

```
Systematic ID  Name     DR    DL  gRNA_DR  gRNA_DL  cluster
 SPAC1002.03c  gls2  0.024 0.000    0.067    9.479        9
 SPAC1002.04c taf11  0.733 2.629    0.895    3.399        4
```

`gls2` 基因层 DL 为 0、gRNA 层为 9.479——压成一列会丢掉这个差异。

**只取 DR/DL 两列**（用户明确选择）：不带 `A`/`t50`/`auc`/`R2` 等拟合参数。

**源表意外地已是每基因一行。** 尽管 `ID` 列带 `_43`/`_213` 后缀、名字里有 "gRNA"，
实测 5060 行对应 5060 个唯一系统名，无重复——所以不需要任何跨 gRNA 聚合，直接 join。
仍然保留「重复 ID 就抛错」的检查：若将来换成真正的多 gRNA 表，静默 fan-out 会放大参考表行数。

读取器同时接受 `um`/`lam` 与 `DR`/`DL` 两种列名，以便上游改名后无需改代码。
10 个 curated 表里的 ID 已不在 PomBase 2026-06-01（已退役），被正确排除（5060 → 5050）。

合计 **24 列**注释。

## 7. CLI 接口、错误处理、测试

```bash
python workflow/scripts/annotate/annotate_pombe_genes.py \
    --input my_gene_list.tsv \
    --gene-column gene_systematic_id \
    --annotation-reference results/annotation/2026-06-01/2026-08-11/gene_annotation_reference.parquet \
    --output my_gene_list.annotated.tsv
```

输入支持 tsv/csv/xlsx（走已有的 `read_file`），输出按扩展名定。

> **与原设计的偏差：`--annotation-reference` 是必填，没有默认值。** 原设计想默认「默认
> pombase version + 最新 sgd 目录」，但那要在脚本里扫目录挑「最新」，会让同一条命令在不同
> 时间产出不同注释（新 fetch 一次 SGD 就静默换了数据源）。显式传路径让参考表版本可追溯。

可选开关：`--columns` 只取部分注释列；`--keep-unmatched` / `--drop-unmatched` 控制未匹配行
是保留（注释列为空）还是丢弃，**默认保留**——丢行是静默数据丢失。

**错误处理（真实会踩的坑）**

1. **未匹配 ID 必须报摘要。** 不做 ID 解析（用户表格已是最新系统名），但仍统计并 log：
   匹配 N 行、未匹配 K 行，并列出未匹配的 ID（>20 个只列前 20）。写错列名或混进非编码基因时
   这是唯一的信号，绝不静默。
2. **输入列有重复 ID。** 用 `merge` 而非 `set_index().join()`，重复行各自拿到注释、行数不变，
   并 log 重复数量。
3. **参考表缺列。** 若参考表是旧版本、缺了新增列，直接报错退出，不产出半张表。
4. **列名撞车不覆盖。** 输入表已有同名列（如自己的 `gene_name`）时，保留用户的，注释列加
   `_annotation` 后缀。
5. **实现期发现：错误信息不能被 traceback 埋掉。** `--gene-column` 写错时，原本
   `summarise_match()` 先用裸 pandas 索引取列，抛出的是无信息的 `KeyError: 'WrongName'`，
   而 `@logger.catch` 又先打印了 60 行 traceback。改为两个入口共用一个前置守卫
   `_require_gene_column()`（报出可用列名），并移除 `run()` 上的 `@logger.catch`。
   现在输出是一行可操作信息 + 退出码 1。
6. **实现期发现：计数列的 dtype。** 参考表 join 后，未覆盖基因引入 NaN 会把 int 提升成
   float，导出成 `1.0`。所有 `*_count` 列在组装后统一转为可空 `Int64`。

**测试** — `tests/test_annotation.py`，实测 **66 个测试**，覆盖：

1. ortholog 字符串解析：`YBR265W(N)+YDR302W(C)`、`YEL066W|YPR193C`、`NONE` 三种形态，
   以及同序拼接后 name / essentiality / qualifier 三列长度一致。
2. `conflicting` 判定：构造同一基因既有 `null+inviable` 又有 `null+viable` 的输入，
   断言输出 `conflicting` 且 evidence 计数正确；分级表型（`viability: decreased`）被排除。
3. SGD 残缺文件读取：14 / 15 / 13 字段三种行形态都不丢行。
4. join 行为：重复 ID 不改变行数、未匹配行保留且注释为空、`--drop-unmatched` 正确丢弃、
   行序保留、同名列加后缀不覆盖、错列名报出可用列名。
5. 各注释块：pombe 侧回落与两套必需性、complex/PFAM 去重、GO-slim 三 namespace 分列与排序、
   gRNA 层 um/lam→gRNA_DR/gRNA_DL 改名与重复 ID 报错。
6. 组装：pombe 基因集定义行集合（外部多余基因不引入）、缺失块留空、计数列保持 Int64。

与仓库现有测试一致：用小的手造 DataFrame，不碰 Snakemake、不读真实大文件。

## 8. 验证结果

- **测试：** 66 个注释测试通过；整体 237 passed。3 个 `test_io_duplicate_columns` 失败与
  2 个模块无法 collect（缺 `dill`/`requests`）均为**改动前既有**问题，已用 `git stash`
  对照确认。
- **参考表：** 经 `snakemake --use-conda` 构建为 12,685 基因 × 24 列。
- **端到端：** 注释 `results/clustering/HD_DIT_HAP/gmm_direct9/final_clusters.tsv`，
  4513/4513 行全部匹配。跨物种必需性一致性合理：缺失文库判 E 的里 662 个在酿酒酵母也
  inviable、判 V 的里 1533 个也 viable。
- **数值正确性：** `gRNA_DR`/`gRNA_DL` 与源表在全部 5050 个基因上逐行相同；
  ortholog 对齐不变量在 5113 行上 0 处错位。
- **已知残留：** `snakemake -n` 会因 `biopython.yml` 环境定义变更报一次
  provenance 触发的重建（`--rerun-triggers mtime` 下无事可做）。这与本次改动无关，
  同环境的既有规则一样受影响，跑过一次后即稳定。
