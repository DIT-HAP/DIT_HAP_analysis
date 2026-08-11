# Pombe 基因注释工具（给任意表格加注释列）

**日期：** 2026-08-11
**范围：** 中等。新增一个通用 CLI：给定任意含 pombe 系统名列的表格，追加约 21 列注释
（budding yeast 同源基因及其必需性、pombe 自身必需性、功能注释），帮助判断基因功能与必需性。
注释数据源由一条 Snakemake 规则预先解析成一张参考表缓存，CLI 只做 join。

## 1. 决策与取舍（用户已定）

| 决策点 | 选定方案 | 理由 |
|---|---|---|
| 使用形态 | 通用 CLI + 预构建参考表 | 既能进分析流程，也能随手注释手头任何表格 |
| SGD 数据获取 | 手动 fetch 脚本，不进 DAG | 与 `resources/external/` 现状一致（pombase/biogrid 全是手工填的，git-ignored） |
| 注释范围 | 同源基因块 + pombe 必需性 + 功能注释 | 定量背景（表达量/网络/paralog）暂不做（YAGNI），已有 features 矩阵可查 |
| 多同源处理 | 每个同源基因单独给值，同序拼接 | 「其中一个必需」与「两个都必需」含义不同，汇总会丢信息 |
| ID 解析 | **不过 `update_sysIDs()`**，直接按系统名 join | 用户表格用的都是最新 ID |
| japonicus 同源 | 不纳入 | 主要用于保守性打分，对判断功能帮助不大 |

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
   宽表，index = pombe systematic id，21 列注释
          │
          ▼  通用 CLI（手动跑，不进 DAG）
annotate_pombe_genes.py --input any.tsv --gene-column X --output annotated.tsv
```

新增文件：

- `workflow/scripts/annotate/fetch_sgd_data.sh` — `curl` 两个 SGD 文件到
  `resources/external/sgd/{date}/`。手动跑一次，把「怎么拉的」记录下来。
- `workflow/src/annotation/core.py` — 纯函数库：每个注释块一个 `build_*` 函数，
  加 `build_annotation_reference()` 组装、`annotate_table()` 做 join。
- `workflow/scripts/annotate/build_annotation_reference.py` — 构建参考表的 CLI。
- `workflow/scripts/annotate/annotate_pombe_genes.py` — 注释 CLI。
- `workflow/rules/annotate.smk` — 只放构建参考表那条规则。
- `tests/test_annotation.py` — 见 §6。

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

**`conflicting` 的来历（已核实）。** SGD phenotype 表里 `null` 突变体：1245 个基因标
inviable、5223 个标 viable，**其中 202 个基因两种标注都有**（不同文献/背景结论不一致）。
不强行二分，输出 `conflicting` 并由 evidence 列给出各自条目数，让使用者看到证据强度，
而不是被一个可能错的二分标签误导。

判定规则：取 `mutant_type == "null"` 且 `phenotype` 属于 `{inviable, viable}` 的行，
按基因分组计数 → 只有 inviable → `inviable`；只有 viable → `viable`；两者都有 →
`conflicting`；无记录 → 空。

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

合计 **21 列**注释。

## 6. CLI 接口、错误处理、测试

```bash
python workflow/scripts/annotate/annotate_pombe_genes.py \
    --input my_gene_list.tsv \
    --gene-column gene_systematic_id \
    --output my_gene_list.annotated.tsv
```

`--annotation-reference` 有默认值（默认 pombase version + 最新 sgd 目录），常规用法只需三个参数。
输入支持 tsv/csv/xlsx（走已有的 `read_file`），输出按扩展名定。

可选开关：`--columns` 只取部分注释列；`--keep-unmatched` / `--drop-unmatched` 控制未匹配行
是保留（注释列为空）还是丢弃，**默认保留**——丢行是静默数据丢失。

**错误处理（三个真实会踩的坑）**

1. **未匹配 ID 必须报摘要。** 不做 ID 解析（用户表格已是最新系统名），但仍统计并 log：
   匹配 N 行、未匹配 K 行，并列出未匹配的 ID（>20 个只列前 20）。写错列名或混进非编码基因时
   这是唯一的信号，绝不静默。
2. **输入列有重复 ID。** 用 `merge` 而非 `set_index().join()`，重复行各自拿到注释、行数不变，
   并 log 重复数量。
3. **参考表缺列。** 若参考表是旧版本、缺了新增列，直接报错退出，不产出半张表。

**测试** — `tests/test_annotation.py`，覆盖三处真正会错的逻辑：

1. ortholog 字符串解析：`YBR265W(N)+YDR302W(C)`、`YEL066W|YPR193C`、`NONE` 三种形态，
   以及同序拼接后 name / essentiality / qualifier 三列长度一致。
2. `conflicting` 判定：构造同一基因既有 `null+inviable` 又有 `null+viable` 的输入，
   断言输出 `conflicting` 且 evidence 计数正确。
3. join 行为：重复 ID 不改变行数、未匹配行保留且注释为空、`--drop-unmatched` 正确丢弃。

与仓库现有测试一致：用小的手造 DataFrame，不碰 Snakemake、不读真实大文件。
