# =============================================================================
# comparison.smk — Pairwise fitness comparison with other large-scale studies
# =============================================================================
#
# Batch B (requires resources/curated/final_clusters.tsv — the un-buildable,
# human-curated cluster table from the manual finalize_gene_clusters step).
# Per-dataset: merges DIT-HAP data with gRNA, Barseq, integration density,
# colony size, growth rate. Pairwise scatter with KDE overlay + Pearson r stats.
#
# Single rule (no prepare/compute split — the merge + correlate is tiny and
# self-contained, same shape as coverage.smk / verification.smk).

rule compare_large_scale_studies:
    input:
        final_clusters="resources/curated/final_clusters.tsv",
        protein_features=lambda wc: (
            f"results/features/{DATASETS['reference']['pombase_version']}/"
            "pombe_coding_gene_protein_features.tsv"
        ),
        gRNA_data=config.get("comparison", {}).get(
            "gRNA_data_file",
            "resources/curated/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv"
        ),
    output:
        stats="results/comparison/{dataset}/fitness_correlation_stats.tsv",
        figures="results/comparison/{dataset}/pairwise_fitness_comparison.pdf",
    params:
        clip_upper=config.get("comparison", {}).get("clip_upper", 200),
    log:
        "logs/comparison/compare_large_scale_studies_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [comparison] Running pairwise fitness comparison for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/comparison/compare_large_scale_studies.py \
            --final-clusters {input.final_clusters} \
            --protein-features {input.protein_features} \
            --grna-data {input.gRNA_data} \
            --clip-upper {params.clip_upper} \
            --output-stats {output.stats} \
            --output-figures {output.figures} &> {log}
        """
