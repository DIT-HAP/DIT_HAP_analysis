# =============================================================================
# verification.smk — Deletion library phenotype verification
# =============================================================================
#
# Per-dataset: merges gene-level DIT-HAP results with Hayles-2013 deletion
# library categories and the curated essentiality verification table.
# Produces the category donut + DR scatter, plus (notebook §4-5) boxplot/violin
# DR comparisons, four critical-gene groups (each a boxplot + verification
# donut + review TSV), and DIT-HAP-vs-gRNA depletion curves. Single-rule shape
# like coverage.smk — data is tiny and self-contained.
#
# gRNA per-timepoint LFC (for the depletion-curve overlay) is HD-only and lives
# in the legacy pipeline repo, so it's sourced from a per-dataset map — same
# pattern as noncoding_rna.smk's _NONCODING_FITTING. Datasets absent from the
# map render DIT-HAP-only curves (the --grna-timepoints flag is omitted).
_GRNA_TIMEPOINT_DATA = {
    "HD_DIT_HAP": "/data/c/yangyusheng_optimized/DIT_HAP_pipeline/resources/HD_gRNA_data.csv",
}


rule compare_deletion_library:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        gene_timepoints=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/gene_level_fitting_statistics.tsv"
        ),
        deletion_library="resources/curated/deletion_library_categories.xlsx",
        essentiality_verification="resources/curated/essentiality_verification.csv",
    output:
        stats="results/verification/{dataset}/verification_stats.tsv",
        figures="results/verification/{dataset}/deletion_library_comparison.pdf",
        boxplots="results/verification/{dataset}/verification_boxplots.pdf",
        depletion_curves="results/verification/{dataset}/verification_depletion_curves.pdf",
        critical_genes_dir=directory("results/verification/{dataset}/critical_genes"),
    params:
        # Optional gRNA overlay: build the flag only when the dataset is in the
        # map, else pass nothing (DIT-HAP-only curves).
        grna_flag=lambda wc: (
            f"--grna-timepoints {_GRNA_TIMEPOINT_DATA[wc.dataset]}"
            if wc.dataset in _GRNA_TIMEPOINT_DATA else ""
        ),
    log:
        "logs/verification/compare_deletion_library_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [verification] Comparing deletion library for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/verification/compare_deletion_library.py \
            --fitting-results {input.fitting_results} \
            --deletion-library {input.deletion_library} \
            --essentiality-verification {input.essentiality_verification} \
            --output-stats {output.stats} \
            --output-figures {output.figures} \
            --gene-timepoints {input.gene_timepoints} \
            {params.grna_flag} \
            --output-boxplots {output.boxplots} \
            --output-depletion-curves {output.depletion_curves} \
            --output-critical-genes-dir {output.critical_genes_dir} &> {log}
        """
