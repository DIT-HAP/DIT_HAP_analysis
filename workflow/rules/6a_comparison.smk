# =============================================================================
# comparison.smk — Pairwise fitness comparison with other large-scale studies
# =============================================================================
#
# Split into 3 rules so each analysis step is independently re-runnable:
#   prepare_fitness_table    -> fitness_table.parquet intermediate (the merge)
#   compute_comparison_stats -> fitness_correlation_stats.tsv (Pearson + Spearman, BH-FDR)
#   plot_comparison_figures  -> pairwise_fitness_comparison(_pN).pdf + correlation
#                               heatmaps (pearson & spearman)
# plot_comparison_figures reads BOTH the prepared parquet (for the actual data)
# and compute_comparison_stats's stats TSV (for the col_x/col_y pairs that
# survived the per-pair overlap filter), so the PDF panels always match the
# TSV rows even though the two rules now run independently.
#
# Data sources are the shared upstream tables, not per-stage recomputations:
# the DIT-HAP (HD_DIT_HAP_DR) and gRNA (gRNA_DR, sign-flipped by the annotation
# stage itself) metrics come from the 1c_annotate gene annotation reference,
# and the other studies' fitness/depletion columns from the 1b_features merged
# protein-features table. Per-dataset via the annotation reference is not
# applicable (it is HD_DIT_HAP-frozen), so this stage effectively analyses the
# HD_DIT_HAP release regardless of {dataset}.
# Pairwise scatter matrix (multi-page, density-coloured) with Pearson r stats,
# plus clustered 10x10 correlation heatmaps per coefficient.

# Parquet intermediate shared by the stats + figures rules.
_CWORK = "results/comparison/{dataset}/_work"


rule prepare_fitness_table:
    input:
        annotation_reference=_ANNOT_REF,
        protein_features=(
            f"results/1b_features/{DATASETS['reference']['pombase_version']}/"
            "pombe_coding_gene_protein_features.tsv"
        ),
    output:
        fitness_table=f"{_CWORK}/fitness_table.parquet",
    params:
        clip_upper=config.get("comparison", {}).get("clip_upper", 200),
    log:
        "logs/comparison/prepare_fitness_table_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [comparison] Preparing fitness table for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/comparison/prepare_fitness_table.py \
            --annotation-reference {input.annotation_reference} \
            --protein-features {input.protein_features} \
            --clip-upper {params.clip_upper} \
            --output-fitness-table {output.fitness_table} &> {log}
        """


rule compute_comparison_stats:
    input:
        fitness_table=f"{_CWORK}/fitness_table.parquet",
    output:
        stats="results/comparison/{dataset}/fitness_correlation_stats.tsv",
    log:
        "logs/comparison/compute_comparison_stats_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [comparison] Computing correlation stats for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/comparison/compute_comparison_stats.py \
            --fitness-table {input.fitness_table} \
            --output-stats {output.stats} &> {log}
        """


rule plot_comparison_figures:
    input:
        fitness_table=f"{_CWORK}/fitness_table.parquet",
        stats="results/comparison/{dataset}/fitness_correlation_stats.tsv",
    output:
        figures="results/comparison/{dataset}/pairwise_fitness_comparison.pdf",
        heatmap_pearson="results/comparison/{dataset}/correlation_pearson_heatmap.pdf",
        heatmap_spearman="results/comparison/{dataset}/correlation_spearman_heatmap.pdf",
    log:
        "logs/comparison/plot_comparison_figures_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [comparison] Plotting pairwise comparison figures for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/comparison/plot_comparison_figures.py \
            --fitness-table {input.fitness_table} \
            --stats {input.stats} \
            --output-dir results/comparison/{wildcards.dataset} &> {log}
        """
