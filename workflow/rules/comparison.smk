# =============================================================================
# comparison.smk — Pairwise fitness comparison with other large-scale studies
# =============================================================================
#
# Split into 3 rules so each analysis step is independently re-runnable:
#   prepare_fitness_table   -> fitness_table.parquet intermediate (the merge)
#   compute_comparison_stats -> fitness_correlation_stats.tsv (Pearson r/p/n)
#   plot_comparison_figures  -> pairwise_fitness_comparison.pdf
# plot_comparison_figures reads BOTH the prepared parquet (for the actual data)
# and compute_comparison_stats's stats TSV (for the col_x/col_y pairs that
# survived the per-pair overlap filter), so the PDF panels always match the
# TSV rows even though the two rules now run independently.
#
# Batch B. Sources the DIT-HAP gene fitness table from the clustering
# finalize-variant system: final_clusters_path(dataset, selected_variant(dataset))
# (clustering.smk). Buildable variants (direct/auto_merge/grid) produce it under
# results/clustering/final/...; the manual_merge variant needs its curated tsv
# first. Only Systematic ID / A / DR / DL are read here (never the cluster id),
# so any variant's table gives identical results.
# Per-dataset: merges DIT-HAP data with gRNA, Barseq, integration density,
# colony size, growth rate. Pairwise scatter with KDE overlay + Pearson r stats.
# Ported from compare_with_other_large_scale_studies.ipynb.

# Parquet intermediate shared by the stats + figures rules.
_CWORK = "results/comparison/{dataset}/_work"


rule prepare_fitness_table:
    input:
        final_clusters=lambda wc: final_clusters_path(wc.dataset, selected_variant(wc.dataset)),
        protein_features=lambda wc: (
            f"results/features/{DATASETS['reference']['pombase_version']}/"
            "pombe_coding_gene_protein_features.tsv"
        ),
        gRNA_data=config.get("comparison", {}).get(
            "gRNA_data_file",
            "resources/curated/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv"
        ),
    output:
        fitness_table=f"{_CWORK}/fitness_table.parquet",
    params:
        clip_upper=config.get("comparison", {}).get("clip_upper", 200),
    log:
        "logs/comparison/prepare_fitness_table_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [comparison] Preparing fitness table for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/comparison/prepare_fitness_table.py \
            --final-clusters {input.final_clusters} \
            --protein-features {input.protein_features} \
            --grna-data {input.gRNA_data} \
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
        "../envs/statistics_and_figure_plotting.yml"
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
    log:
        "logs/comparison/plot_comparison_figures_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [comparison] Plotting pairwise comparison figures for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/comparison/plot_comparison_figures.py \
            --fitness-table {input.fitness_table} \
            --stats {input.stats} \
            --output-figures {output.figures} &> {log}
        """
