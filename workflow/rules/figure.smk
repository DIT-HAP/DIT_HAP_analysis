"""
figure.smk - 绘图规则集中管理
================================

本文件集中管理所有图表生成规则（输出 PDF/PNG 的规则）。

职责分离：
- 计算规则（生成数据）保留在各自的领域 .smk 文件中（clustering.smk, enrichment.smk 等）
- 绘图规则（消费数据，输出图表）统一放在此文件中

命名约定：
- 规则名：plot_<功能>（如 plot_coherence, plot_variant_clusters）
- 环境声明：conda: "../envs/cnsplots.yml"
- 产物路径：直接描述内容，不加类型后缀（coherence.pdf，不是 coherence_figure.pdf）

使用模式：
    rule plot_example:
        conda: "../envs/cnsplots.yml"
        input: "results/example/{dataset}/data.parquet"
        output: "results/example/{dataset}/example.pdf"
        script: "../scripts/example/plot_example.py"

Author: Yusheng Yang
Date: 2026-09-02
"""

# 绘图规则将在此处添加
# 例如：
# rule plot_coherence:
#     conda: "../envs/cnsplots.yml"
#     input: "results/coherence/{dataset}/coherence.parquet"
#     output: "results/coherence/{dataset}/coherence.pdf"
#     script: "../scripts/coherence/plot_coherence.py"
