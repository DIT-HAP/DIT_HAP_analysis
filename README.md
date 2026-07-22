# Pickle to Parquet Migration - README

本目录包含了将项目中的 pickle 中间文件迁移到 parquet 格式的完整工作。

---

## 📋 文档导航

### 快速开始
- **[MIGRATION_SUMMARY.md](MIGRATION_SUMMARY.md)** ⭐ 
  - 最简洁的总结，推荐优先阅读
  - 包含迁移统计、验证结果、使用说明

### 详细文档
1. **[PICKLE_TO_PARQUET_MIGRATION.md](PICKLE_TO_PARQUET_MIGRATION.md)**
   - 迁移分析报告
   - 分类所有 pickle 文件为"可迁移"和"不可迁移"
   - 技术决策和风险分析

2. **[PICKLE_TO_PARQUET_COMPLETION_REPORT.md](PICKLE_TO_PARQUET_COMPLETION_REPORT.md)**
   - 完整的实施报告
   - 代码变更模式、收益分析、测试策略

3. **[CANNOT_MIGRATE_TO_PARQUET.md](CANNOT_MIGRATE_TO_PARQUET.md)**
   - 无法迁移的文件清单（3类）
   - 详细说明原因和替代方案

---

## ✅ 迁移结果

### 成功迁移：11 类文件
- **Features 阶段**: 6 类 DataFrame（dna, rna, protein, evolutionary, network, phenotype）
- **Clustering 阶段**: 4 类文件（annotated_data, scaled_data, k_sweep_metrics, labels）
- **ML 阶段**: 1 类文件（modeling_data）

### 保持 pickle：3 类文件
- **Enrichment 阶段**: genesets.pkl, id2name.pkl, {ontology}_frames.pkl
- **原因**: 复杂 Python 对象（dataclass, dict, 多 DataFrame）
- **详见**: [CANNOT_MIGRATE_TO_PARQUET.md](CANNOT_MIGRATE_TO_PARQUET.md)

---

## 📊 统计数据

| 指标 | 数量 |
|------|------|
| 修改的脚本 | 14 |
| 修改的规则文件 | 3 |
| 修改的核心库 | 2 |
| 修改的测试 | 4 |
| 新增文档 | 5 |
| **总修改文件** | **28** |

---

## 🚀 快速验证

### 检查 Snakemake DAG
```bash
mamba run -n snakemake snakemake -n
# 应该显示: Building DAG of jobs... (无错误)
```

### 运行测试
```bash
mamba activate snakemake
pytest tests/ -v
```

### 重新生成中间文件
```bash
mamba activate snakemake
snakemake --use-conda --cores 8
```

---

## 🔧 使用新的 Parquet API

### 写入
```python
from workflow.src.io import write_parquet

write_parquet(df, output_path)
```

### 读取
```python
from workflow.src.io import read_parquet

df = read_parquet(input_path)
```

---

## 💡 为什么迁移到 Parquet？

### 优势
1. ✅ **跨语言兼容** - R, Python, Java, Spark 都能读取
2. ✅ **更高压缩率** - 列式存储 + Snappy 压缩，文件减少 30-60%
3. ✅ **更快读写** - 列选择性读取，多核并行
4. ✅ **类型安全** - Schema 验证，防止数据损坏
5. ✅ **安全性** - 无代码执行风险（pickle 可执行任意代码）

---

## 📝 关键决策

### 为什么不迁移 Enrichment 阶段？
1. **技术限制**: 存储的是复杂 Python 对象（dataclass, dict），不是 DataFrame
2. **投入产出比**: 迁移需要重构 3 个脚本，收益不明显
3. **使用场景**: 这些文件只在 Python 中使用，不需要跨语言
4. **风险可控**: 中间文件，非最终输出，安全风险较低

---

## 🎯 验证清单

- [x] ✅ 所有脚本语法正确
- [x] ✅ Snakemake DAG 构建成功
- [x] ✅ 所有路径更新为 .parquet
- [x] ✅ 核心库添加 parquet I/O 函数
- [x] ✅ 测试文件已更新
- [ ] ⏳ 运行完整测试套件（需要 snakemake 环境）
- [ ] ⏳ 重新生成中间文件验证兼容性

---

**最后更新**: 2026-07-22  
**执行环境**: Worktree `.claude/worktrees/pickle-to-parquet`
