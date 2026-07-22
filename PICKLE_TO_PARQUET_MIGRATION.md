# Pickle 到 Parquet 迁移分析报告

## 概述

本文档分析了项目中所有 pickle 文件的使用情况，并确定哪些可以安全地迁移到 parquet 格式。

## 分类汇总

### ✅ 可以迁移的文件 (DataFrames/Series)

这些文件存储的是 pandas DataFrame 或 Series，可以直接迁移到 parquet：

#### 1. Features 阶段的中间文件
- `results/features/{version}/_levels/dna_features.pkl` → DataFrame
- `results/features/{version}/_levels/rna_features.pkl` → DataFrame
- `results/features/{version}/_levels/protein_features.pkl` → DataFrame
- `results/features/{version}/_levels/evolutionary_features.pkl` → DataFrame
- `results/features/{version}/_levels/network_features.pkl` → DataFrame
- `results/features/{version}/_levels/phenotype_features.pkl` → DataFrame

**迁移理由**: 这些都是特征表，纯 DataFrame 格式，无重复列，无特殊对象。

#### 2. Clustering 数据文件
- `results/clustering/{dataset}/_work/annotated_data.pkl` → DataFrame
  - 完整的拟合统计表 + RevisedDeletion_essentiality，index = systematic ID
- `results/clustering/{dataset}/_work/scaled_data.pkl` → DataFrame
  - 缩放后的 (DR, DL) 矩阵，用于聚类
- `results/clustering/{dataset}/_work/k_sweep_metrics.pkl` → DataFrame
  - KMeans k-sweep 指标（inertia + silhouette/CH/DB per k）
- `results/clustering/{dataset}/{variant}/_labels.pkl` → Series
  - 聚类标签，pd.Series (0-based cluster labels)

**迁移理由**: 都是结构化的表格数据，无复杂对象。

#### 3. ML 准备数据
- `results/ml/models/{dataset}/{version}/_work/modeling_data.pkl` → DataFrame
  - 合并的、DR过滤后的建模表

**迁移理由**: 纯 DataFrame。

---

### ❌ 不能迁移的文件 (复杂Python对象)

这些文件存储的是复杂的 Python 对象，parquet 不支持：

#### 1. Enrichment 中间对象
- `results/enrichment/raw/{dataset}/{variant}/_work/genesets.pkl`
  - 类型: `ClusterGeneSets` dataclass
  - 内容: `{cluster_genes: dict, bg_genes: list, nonwt_bg_genes: list}`
  - **不能迁移原因**: 自定义 dataclass，包含嵌套字典和列表

- `results/enrichment/raw/{dataset}/{variant}/_work/id2name.pkl`
  - 类型: `dict[str, str]` (systematic ID -> gene name mapping)
  - **不能迁移原因**: 纯字典对象，不是 DataFrame

- `results/enrichment/raw/{dataset}/{variant}/_work/{GO,FYPO,MONDO}_frames.pkl`
  - 类型: `dict` 包含4个DataFrame: `{"full": df, "slim": df, "nonwt_full": df, "nonwt_slim": df}`
  - **不能迁移原因**: 字典包装的多个 DataFrame，不是单一 DataFrame

**替代方案**:
- `genesets.pkl`: 可以序列化为 JSON
- `id2name.pkl`: 可以序列化为 JSON 或转换为 DataFrame 后存 parquet
- `{ontology}_frames.pkl`: 可以拆分成4个独立的 parquet 文件

---

## 迁移策略

### 阶段1: 直接替换 (DataFrame/Series → Parquet)

对于所有纯 DataFrame/Series 的文件，执行以下替换：

**代码修改模式**:
```python
# 旧代码
df.to_pickle(path)
df = pd.read_pickle(path)

# 新代码
df.to_parquet(path, engine='pyarrow', compression='snappy', index=True)
df = pd.read_parquet(path, engine='pyarrow')
```

**文件扩展名**: `.pkl` → `.parquet`

**涉及的文件**:
- `workflow/scripts/features/*.py` (6个脚本)
- `workflow/scripts/clustering/prepare_clustering_data.py`
- `workflow/scripts/clustering/cluster_one_method.py`
- `workflow/scripts/clustering/finalize_*.py` (3个脚本)
- `workflow/scripts/ml/prepare_ml_data.py`
- `workflow/scripts/ml/train_automl.py`
- `workflow/rules/features.smk`
- `workflow/rules/clustering.smk`
- `workflow/rules/ml.smk`
- `tests/test_*.py` (相关测试)

### 阶段2: 结构化替换 (复杂对象 → 替代格式)

#### 2.1 ClusterGeneSets (genesets.pkl)

**当前结构**:
```python
@dataclass
class ClusterGeneSets:
    cluster_genes: dict  # {int -> [str]}
    bg_genes: list
    nonwt_bg_genes: list
```

**替代方案**: JSON
```python
# 写入
with open(path.with_suffix('.json'), 'w') as f:
    json.dump({
        'cluster_genes': genesets.cluster_genes,
        'bg_genes': genesets.bg_genes,
        'nonwt_bg_genes': genesets.nonwt_bg_genes
    }, f)

# 读取
with open(path, 'r') as f:
    data = json.load(f)
    # 需要转换 cluster_genes 的 key 为 int
    data['cluster_genes'] = {int(k): v for k, v in data['cluster_genes'].items()}
genesets = ClusterGeneSets(**data)
```

#### 2.2 id2name 字典 (id2name.pkl)

**替代方案1**: JSON
```python
with open(path, 'w') as f:
    json.dump(id2name, f)
```

**替代方案2**: Parquet (转为DataFrame)
```python
# 写入
pd.DataFrame(list(id2name.items()), columns=['systematic_id', 'name']).to_parquet(path)

# 读取
df = pd.read_parquet(path)
id2name = dict(zip(df['systematic_id'], df['name']))
```

#### 2.3 Ontology frames 字典 ({ontology}_frames.pkl)

**当前结构**: `{"full": df, "slim": df, "nonwt_full": df, "nonwt_slim": df}`

**替代方案**: 拆分成独立文件
```python
# 写入
base = work_dir / f"{ontology}_frames"
base.mkdir(exist_ok=True)
frames['full'].to_parquet(base / 'full.parquet')
frames['slim'].to_parquet(base / 'slim.parquet')
frames['nonwt_full'].to_parquet(base / 'nonwt_full.parquet')
frames['nonwt_slim'].to_parquet(base / 'nonwt_slim.parquet')

# 读取
base = work_dir / f"{ontology}_frames"
frames = {
    'full': pd.read_parquet(base / 'full.parquet'),
    'slim': pd.read_parquet(base / 'slim.parquet'),
    'nonwt_full': pd.read_parquet(base / 'nonwt_full.parquet'),
    'nonwt_slim': pd.read_parquet(base / 'nonwt_slim.parquet')
}
```

---

## 收益分析

### Parquet 的优势
1. **跨语言兼容**: R, Python, Java, Spark 等都能读取
2. **更高压缩率**: 列式存储 + Snappy/ZSTD 压缩
3. **更快的读写**: 特别是选择性列读取
4. **类型保留**: 保留 int/float/datetime 等类型信息
5. **Schema 验证**: 内置 schema，防止数据损坏

### Pickle 的劣势
1. **Python 专属**: 只能在 Python 中使用
2. **安全风险**: 反序列化可执行任意代码
3. **版本兼容性**: Python 版本不同可能导致无法读取
4. **不可读**: 二进制格式，无法检查

---

## 不能迁移的理由总结

### ClusterGeneSets (genesets.pkl)
- **为什么不能用 Parquet**: 自定义 dataclass，包含嵌套字典 `{int: [str]}` 和列表
- **Parquet 限制**: 只支持表格数据（行列结构），不支持嵌套字典
- **推荐替代**: JSON (保持结构) 或重构为 DataFrame

### id2name (id2name.pkl)
- **为什么不能用 Parquet**: 纯字典对象
- **Parquet 限制**: 需要是 DataFrame 或类似表格结构
- **推荐替代**: 转为 DataFrame 后存 parquet，或用 JSON

### Ontology frames ({ontology}_frames.pkl)
- **为什么不能用 Parquet**: 字典包装的多个 DataFrame
- **Parquet 限制**: 一个 parquet 文件只能存一个 DataFrame
- **推荐替代**: 拆分成4个独立 parquet 文件（更符合 parquet 的设计理念）

---

## 实施优先级

### P0 (立即执行): 纯 DataFrame/Series 迁移
- Features 阶段 (6个文件)
- Clustering 阶段 (4个文件)
- ML 阶段 (1个文件)

**影响**: 最大收益，无风险，代码改动最小

### P1 (考虑执行): 结构化对象重构
- Enrichment 阶段的3类文件

**影响**: 需要重新设计数据结构，测试成本较高

---

## 实施计划

1. ✅ **创建 worktree** (已完成)
2. ⏳ **迁移 Features 阶段** (6个脚本 + 1个规则文件)
3. ⏳ **迁移 Clustering 阶段** (4个脚本 + 1个规则文件)
4. ⏳ **迁移 ML 阶段** (2个脚本 + 1个规则文件)
5. ⏳ **更新所有测试** (5个测试文件)
6. ⏳ **运行完整测试套件**
7. ⏳ **文档更新** (CLAUDE.md 中关于 pickle 的说明)

---

## 测试策略

对于每个迁移的文件:
1. 运行相关的 pytest 测试
2. 使用 Snakemake dry-run 验证规则完整性
3. 对比迁移前后的输出文件内容 (通过 pandas.testing.assert_frame_equal)
4. 确保字节级精度保持不变 (DR cap, DL divisor 等quirks)

---

## 风险与缓解

### 风险1: Parquet 不支持重复列名
- **影响**: CLAUDE.md 明确提到某些中间文件可能有 "intentional duplicate columns"
- **缓解**: 在迁移前检查所有 DataFrame 是否有重复列，如有则保留 pickle

### 风险2: 浮点精度差异
- **影响**: Parquet 的浮点存储可能与 pickle 有微小差异
- **缓解**: 使用高精度存储选项，测试时使用 `rtol=1e-10` 比较

### 风险3: Index 处理
- **影响**: Parquet 对 MultiIndex 的支持与 pickle 不同
- **缓解**: 确保所有 `to_parquet` 调用明确指定 `index=True`

---

## 后续工作

如果用户希望迁移 Enrichment 阶段的复杂对象，需要:
1. 重构 `ClusterGeneSets` 为多个 DataFrame
2. 将 `id2name` 转换为 DataFrame 格式
3. 拆分 ontology frames 字典为独立文件
4. 更新所有相关的读写逻辑
