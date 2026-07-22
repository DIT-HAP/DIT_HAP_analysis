# 无法迁移到 Parquet 的文件清单

本文档列出了所有**无法**从 pickle 迁移到 parquet 的文件，并详细说明原因。

---

## 不能迁移的文件（3类）

### 1. `genesets.pkl`

**位置**：`results/enrichment/raw/{dataset}/{variant}/_work/genesets.pkl`

**生成脚本**：`workflow/scripts/enrichment/prepare_genesets.py`

**存储内容**：
```python
@dataclass
class ClusterGeneSets:
    cluster_genes: dict[int, list[str]]  # {cluster_id: [gene_ids]}
    bg_genes: list[str]                  # 背景基因列表
    nonwt_bg_genes: list[str]            # 非WT背景基因列表
```

**当前代码**：
```python
# Line 113-114 in prepare_genesets.py
pd.to_pickle(genesets, config.work_dir / "genesets.pkl")
```

**为什么不能用 Parquet**：
- ❌ Parquet 只支持表格数据（行列结构）
- ❌ 不支持嵌套的字典结构 `{int: list[str]}`
- ❌ 不支持 Python dataclass 对象

**可行的替代方案**：

#### 方案 A：重构为多个 DataFrame
```python
# 写入
# 1. cluster_genes: 长格式表
cluster_df = pd.DataFrame([
    {"cluster": cluster_id, "gene_id": gene_id}
    for cluster_id, genes in genesets.cluster_genes.items()
    for gene_id in genes
])
write_parquet(cluster_df, config.work_dir / "cluster_genes.parquet")

# 2. bg_genes: 单列表
write_parquet(
    pd.DataFrame({"gene_id": genesets.bg_genes}),
    config.work_dir / "bg_genes.parquet"
)

# 3. nonwt_bg_genes: 单列表
write_parquet(
    pd.DataFrame({"gene_id": genesets.nonwt_bg_genes}),
    config.work_dir / "nonwt_bg_genes.parquet"
)

# 读取
cluster_df = read_parquet(config.work_dir / "cluster_genes.parquet")
cluster_genes = cluster_df.groupby("cluster")["gene_id"].apply(list).to_dict()
cluster_genes = {int(k): v for k, v in cluster_genes.items()}  # 确保 key 是 int

bg_genes = read_parquet(config.work_dir / "bg_genes.parquet")["gene_id"].tolist()
nonwt_bg_genes = read_parquet(config.work_dir / "nonwt_bg_genes.parquet")["gene_id"].tolist()

genesets = ClusterGeneSets(
    cluster_genes=cluster_genes,
    bg_genes=bg_genes,
    nonwt_bg_genes=nonwt_bg_genes
)
```

**优点**：跨语言兼容，列式存储高效  
**缺点**：代码复杂度增加，需要修改 3 个脚本

#### 方案 B：使用 JSON
```python
import json

# 写入
with open(config.work_dir / "genesets.json", "w") as f:
    json.dump({
        "cluster_genes": genesets.cluster_genes,
        "bg_genes": genesets.bg_genes,
        "nonwt_bg_genes": genesets.nonwt_bg_genes,
    }, f, indent=2)

# 读取
with open(config.work_dir / "genesets.json", "r") as f:
    data = json.load(f)
    # JSON 会将 dict key 转为字符串，需要转回 int
    data["cluster_genes"] = {int(k): v for k, v in data["cluster_genes"].items()}
    genesets = ClusterGeneSets(**data)
```

**优点**：人类可读，版本控制友好，跨语言  
**缺点**：比 pickle/parquet 慢，文件更大

---

### 2. `id2name.pkl`

**位置**：`results/enrichment/raw/{dataset}/{variant}/_work/id2name.pkl`

**生成脚本**：`workflow/scripts/enrichment/prepare_genesets.py`

**存储内容**：
```python
id2name: dict[str, str]  # {systematic_id: gene_name}
```

**当前代码**：
```python
# Line 114 in prepare_genesets.py
pd.to_pickle(gene_meta.id2name, config.work_dir / "id2name.pkl")
```

**为什么不能用 Parquet**：
- ❌ Parquet 需要表格结构（DataFrame/Series）
- ❌ 纯字典对象不符合 parquet 的数据模型

**可行的替代方案**：

#### 方案 A：转为 DataFrame
```python
# 写入
id2name_df = pd.DataFrame(
    list(gene_meta.id2name.items()),
    columns=["systematic_id", "gene_name"]
)
write_parquet(id2name_df, config.work_dir / "id2name.parquet")

# 读取
id2name_df = read_parquet(config.work_dir / "id2name.parquet")
id2name = dict(zip(id2name_df["systematic_id"], id2name_df["gene_name"]))
```

**优点**：标准 parquet 格式，跨语言兼容  
**缺点**：转换开销小（但可忽略）

#### 方案 B：使用 JSON
```python
# 写入
with open(config.work_dir / "id2name.json", "w") as f:
    json.dump(gene_meta.id2name, f, indent=2)

# 读取
with open(config.work_dir / "id2name.json", "r") as f:
    id2name = json.load(f)
```

**优点**：人类可读  
**缺点**：比 parquet 慢

---

### 3. `{GO,FYPO,MONDO}_frames.pkl`

**位置**：`results/enrichment/raw/{dataset}/{variant}/_work/{ontology}_frames.pkl`

**生成脚本**：`workflow/scripts/enrichment/enrich_one_ontology.py`

**存储内容**：
```python
frames: dict[str, pd.DataFrame] = {
    "full": full_enrichment_df,
    "slim": slim_enrichment_df,
    "nonwt_full": nonwt_full_enrichment_df,
    "nonwt_slim": nonwt_slim_enrichment_df,
}
```

**当前代码**：
```python
# Line 127-130 in enrich_one_ontology.py
pd.to_pickle(
    {"full": full, "slim": slim, "nonwt_full": nonwt_full, "nonwt_slim": nonwt_slim},
    config.work_dir / f"{config.ontology}_frames.pkl",
)
```

**为什么不能用 Parquet**：
- ❌ Parquet 一个文件只能存一个 DataFrame
- ❌ 不支持字典包装的多个 DataFrame

**可行的替代方案**：

#### 方案 A：拆分为独立的 parquet 文件
```python
# 写入
base_dir = config.work_dir / f"{config.ontology}_frames"
base_dir.mkdir(exist_ok=True)

write_parquet(full, base_dir / "full.parquet")
write_parquet(slim, base_dir / "slim.parquet")
write_parquet(nonwt_full, base_dir / "nonwt_full.parquet")
write_parquet(nonwt_slim, base_dir / "nonwt_slim.parquet")

# 读取（在 finalize_enrichment.py 中）
frames = {}
for onto in ONTOLOGIES:
    base_dir = config.work_dir / f"{onto}_frames"
    frames[onto] = {
        "full": read_parquet(base_dir / "full.parquet"),
        "slim": read_parquet(base_dir / "slim.parquet"),
        "nonwt_full": read_parquet(base_dir / "nonwt_full.parquet"),
        "nonwt_slim": read_parquet(base_dir / "nonwt_slim.parquet"),
    }
```

**优点**：
- 标准 parquet 格式，跨语言兼容
- 可以单独读取某个变体（不需要加载全部4个）
- 更符合 parquet 的设计理念

**缺点**：
- 目录结构从 1 个文件变为 1 个目录 + 4 个文件
- 需要修改 2 个脚本（写入和读取）

#### 方案 B：合并为单个 DataFrame + 标签列
```python
# 写入
combined = pd.concat([
    full.assign(variant="full"),
    slim.assign(variant="slim"),
    nonwt_full.assign(variant="nonwt_full"),
    nonwt_slim.assign(variant="nonwt_slim"),
], ignore_index=True)
write_parquet(combined, config.work_dir / f"{config.ontology}_frames.parquet")

# 读取
df = read_parquet(config.work_dir / f"{config.ontology}_frames.parquet")
frames = {
    "full": df[df["variant"] == "full"].drop(columns="variant"),
    "slim": df[df["variant"] == "slim"].drop(columns="variant"),
    "nonwt_full": df[df["variant"] == "nonwt_full"].drop(columns="variant"),
    "nonwt_slim": df[df["variant"] == "nonwt_slim"].drop(columns="variant"),
}
```

**优点**：单个文件，简化目录结构  
**缺点**：必须读取全部数据才能获取单个变体

---

## 受影响的脚本

如果要迁移这 3 类文件，需要修改以下脚本：

### 写入端
1. `workflow/scripts/enrichment/prepare_genesets.py`
   - `genesets.pkl` (Line 113)
   - `id2name.pkl` (Line 114)

2. `workflow/scripts/enrichment/enrich_one_ontology.py`
   - `{ontology}_frames.pkl` (Line 127-130)

### 读取端
1. `workflow/scripts/enrichment/enrich_one_ontology.py`
   - `genesets.pkl` (Line 113)
   - `id2name.pkl` (Line 114)

2. `workflow/scripts/enrichment/finalize_enrichment.py`
   - `{ontology}_frames.pkl` (Line 89)

---

## 推荐做法

### 立即采取的行动：❌ 不迁移

**理由**：
1. **投入产出比低**：Enrichment 阶段的这些文件不是性能瓶颈
2. **复杂度增加**：迁移需要修改 3 个脚本，引入额外的复杂性
3. **pickle 已足够**：这些文件只在 Python 环境中使用，不需要跨语言
4. **风险可控**：这些是中间文件，不是最终输出，pickle 的安全风险较低

### 未来考虑的场景：✅ 可选迁移

如果出现以下情况，可以考虑迁移：

1. **需要跨语言访问**：例如，想用 R 或 Julia 读取这些数据
2. **性能成为瓶颈**：enrichment 阶段变慢，需要优化 I/O
3. **版本控制需求**：需要人类可读的格式进行 diff（用 JSON）
4. **标准化要求**：整个项目要求统一的数据格式

---

## 技术债务追踪

如果未来决定迁移这些文件，可以参考本文档中的替代方案。推荐顺序：

1. **优先级 1**: `id2name.pkl` → DataFrame + parquet（最简单，收益明确）
2. **优先级 2**: `{ontology}_frames.pkl` → 拆分为 4 个 parquet（结构清晰）
3. **优先级 3**: `genesets.pkl` → DataFrame + parquet（最复杂，需要重构逻辑）

---

## 总结

| 文件类型 | 不能迁移原因 | 推荐替代方案 | 迁移优先级 |
|---------|-------------|-------------|-----------|
| `genesets.pkl` | 嵌套字典 + dataclass | 多个 DataFrame | 低 (P3) |
| `id2name.pkl` | 纯字典对象 | DataFrame 或 JSON | 中 (P1) |
| `{ontology}_frames.pkl` | 字典包装多个 DataFrame | 拆分为独立文件 | 中 (P2) |

**当前决策**：保持使用 pickle，文档化原因，留待未来优化。

**文档日期**：2026-07-22  
**更新人**：Claude Opus 4.8
