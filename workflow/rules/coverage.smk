# =============================================================================
# coverage.smk — Gene insertion coverage statistics
# =============================================================================
#
# Per-dataset: computes insertion coverage fractions (in-gene vs intergenic)
# and gene coverage (covered vs not covered) per essentiality class AND per
# characterisation_status. Uses the exact IN_GENE_FILTER string from the source
# notebook (quirk).
#
# Split into 3 rules:
#   prepare_coverage_data   -> annotations / gene_result parquet intermediates
#                              (the single fan-out point). gene_result is built
#                              from the FULL protein-coding gene universe (PomBase
#                              metadata) left-joined to the fitting results, so
#                              uncovered genes survive as DR=NaN rows; essentiality
#                              for uncovered genes is backfilled from the curated
#                              deletion_library_categories.xlsx.
#   compute_coverage_stats  -> coverage_stats.tsv + detailed_genes.xlsx
#   plot_coverage_figures   -> coverage_figures.pdf
# plot_coverage_figures reads coverage_stats.tsv (donuts/bars) + gene_result
# parquet (per-gene DR/DL histograms), so the figures always match the numbers
# in the stats table. Editing the stats rule therefore rebuilds the figures too —
# a deliberate coupling for figure/table agreement.

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
        gene_metadata=lambda wc: (
            f"resources/external/pombase/{DATASETS['reference']['pombase_version']}/"
            "Gene_metadata/gene_ids_and_details.parquet"
        ),
        deletion_library_xlsx="resources/curated/deletion_library_categories.xlsx",
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
            --gene-metadata {input.gene_metadata} \
            --deletion-library-xlsx {input.deletion_library_xlsx} \
            --output-annotations {output.annotations} \
            --output-gene-result {output.gene_result} &> {log}
        """


rule compute_coverage_stats:
    input:
        annotations=f"{_COVWORK}/annotations.parquet",
        gene_result=f"{_COVWORK}/gene_result.parquet",
        gene_metadata=lambda wc: (
            f"resources/external/pombase/{DATASETS['reference']['pombase_version']}/"
            "Gene_metadata/gene_ids_and_details.parquet"
        ),
    output:
        stats="results/coverage/{dataset}/coverage_stats.tsv",
        detailed_genes_xlsx="results/coverage/{dataset}/detailed_genes.xlsx",
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
            --gene-metadata {input.gene_metadata} \
            --output-stats {output.stats} \
            --output-detailed-genes-xlsx {output.detailed_genes_xlsx} &> {log}
        """



# -----------------------------------------------------------------------------
# Stage 2b: Plot figures
# -----------------------------------------------------------------------------
# Moved to figure.smk: rule plot_coverage_figures
# Reads coverage_stats.tsv + gene_result.parquet -> coverage_figures.pdf


rule plot_coverage_figures:
    input:
        stats="results/coverage/{dataset}/coverage_stats.tsv",
        gene_result=f"{_COVWORK}/gene_result.parquet",
    output:
        figures="results/coverage/{dataset}/coverage_figures.pdf",
    log:
        "logs/coverage/plot_coverage_figures_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [coverage] Plotting coverage figures for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coverage/plot_coverage_figures.py \
            --stats {input.stats} \
            --gene-result {input.gene_result} \
            --output-figures {output.figures} &> {log}
        """
