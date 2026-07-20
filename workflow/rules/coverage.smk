# =============================================================================
# coverage.smk — Gene insertion coverage statistics
# =============================================================================
#
# Per-dataset: computes insertion coverage fractions (in-gene vs intergenic)
# and gene coverage (covered vs not covered) per essentiality class.
# Uses the exact IN_GENE_FILTER string from the source notebook (quirk).
#
# Single rule (no prepare/compute split — data is tiny and self-contained,
# same shape as spikein.smk). Gene-level fitting_results.tsv already carries
# DeletionLibrary_essentiality as a native column, so unlike clustering.smk
# there's no essentiality_verification.csv merge needed here.

rule compute_coverage_stats:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/insertion_level/fitting_results.tsv"
        ),
        annotations=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/insertion_level/annotations.tsv.gz"
        ),
        gene_level=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
    output:
        stats="results/coverage/{dataset}/coverage_stats.tsv",
        figures="results/coverage/{dataset}/coverage_figures.pdf",
    log:
        "logs/coverage/compute_coverage_stats_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coverage] Computing insertion + gene coverage for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coverage/compute_coverage_stats.py \
            --fitting-results {input.fitting_results} \
            --annotations {input.annotations} \
            --gene-level {input.gene_level} \
            --output-stats {output.stats} \
            --output-figures {output.figures} &> {log}
        """
