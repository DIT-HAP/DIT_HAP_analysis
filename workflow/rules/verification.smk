# =============================================================================
# verification.smk — Deletion library phenotype verification
# =============================================================================
#
# Per-dataset: merges gene-level DIT-HAP results with Hayles-2013 deletion
# library categories and the curated essentiality verification table.
# Produces a donut chart of phenotype categories + DR scatter (same single-rule
# shape as coverage.smk — data is tiny and self-contained, no prepare/compute
# split needed).

rule compare_deletion_library:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        deletion_library="resources/curated/deletion_library_categories.xlsx",
        essentiality_verification="resources/curated/essentiality_verification.csv",
    output:
        stats="results/verification/{dataset}/verification_stats.tsv",
        figures="results/verification/{dataset}/deletion_library_comparison.pdf",
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
            --output-figures {output.figures} &> {log}
        """
