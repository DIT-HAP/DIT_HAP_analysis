# =============================================================================
# coverage.smk — Gene insertion coverage statistics
# =============================================================================
#
# Per-dataset: computes insertion coverage fractions (in-gene vs intergenic)
# and gene coverage (covered vs not covered) per essentiality class. Uses the
# exact IN_GENE_FILTER string from the source notebook (quirk).
#
# Split into 3 rules so each analysis step is independently re-runnable:
#   prepare_coverage_data   -> annotations / gene_result parquet intermediates
#                              (the single fan-out point)
#   compute_coverage_stats  -> coverage_stats.tsv
#   plot_coverage_figures   -> coverage_figures.pdf
# The two downstream rules depend only on prepare_coverage_data's output, so
# editing the stats table never forces the figures to rebuild (and vice
# versa). Gene-level fitting_results.tsv already carries
# DeletionLibrary_essentiality as a native column, so unlike clustering.smk
# there's no essentiality_verification.csv merge needed here.

_COVWORK = "results/coverage/{dataset}/_work"


rule prepare_coverage_data:
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
        annotations=f"{_COVWORK}/annotations.parquet",
        gene_result=f"{_COVWORK}/gene_result.parquet",
    log:
        "logs/coverage/prepare_coverage_data_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coverage] Preparing annotations + gene-result tables for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coverage/prepare_coverage_data.py \
            --fitting-results {input.fitting_results} \
            --annotations {input.annotations} \
            --gene-level {input.gene_level} \
            --output-annotations {output.annotations} \
            --output-gene-result {output.gene_result} &> {log}
        """


rule compute_coverage_stats:
    input:
        annotations=f"{_COVWORK}/annotations.parquet",
        gene_result=f"{_COVWORK}/gene_result.parquet",
    output:
        stats="results/coverage/{dataset}/coverage_stats.tsv",
    log:
        "logs/coverage/compute_coverage_stats_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coverage] Computing insertion + gene coverage stats for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coverage/compute_coverage_stats.py \
            --annotations {input.annotations} \
            --gene-result {input.gene_result} \
            --output-stats {output.stats} &> {log}
        """


rule plot_coverage_figures:
    input:
        annotations=f"{_COVWORK}/annotations.parquet",
        gene_result=f"{_COVWORK}/gene_result.parquet",
    output:
        figures="results/coverage/{dataset}/coverage_figures.pdf",
    log:
        "logs/coverage/plot_coverage_figures_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coverage] Plotting coverage figures for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coverage/plot_coverage_figures.py \
            --annotations {input.annotations} \
            --gene-result {input.gene_result} \
            --output-figures {output.figures} &> {log}
        """
