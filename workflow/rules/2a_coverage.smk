# =============================================================================
# 2a_coverage.smk — Gene insertion coverage statistics
# =============================================================================
#
# Per-dataset: computes insertion coverage fractions (in-gene vs intergenic)
# and gene coverage (covered vs not covered) broken down by each of the three
# annotation dimensions every table now carries — characterisation_status,
# FYPOviability and deletion_essentiality (the names are the annotation
# reference's own; see coverage/core.py's DIMENSION_LABELS). Uses the exact
# IN_GENE_FILTER string from the source notebook (quirk).
#
# Split into 3 rules:
#   prepare_coverage_data   -> annotations / gene_result parquet intermediates
#                              (the single fan-out point). gene_result is built
#                              from the FULL protein-coding gene universe + its
#                              annotation columns, both taken from the gene
#                              annotation reference (1c_annotate.smk), left-joined to
#                              the fitting results so uncovered genes survive as
#                              DR=NaN rows. Only the reference's annotation columns
#                              are read — its baked-in HD_DIT_HAP DR/DL and the
#                              SGD-derived blocks are deliberately not (a dataset's
#                              DR/DL must come from that dataset's own release).
#   compute_coverage_stats  -> coverage_stats.tsv + detailed_genes.xlsx
#   plot_coverage_figures   -> ten single-figure PDFs (cnsplots style): for EACH
#                              dimension a coverage composition plus a DR and a DL
#                              histogram, and insertion placement per chromosome.
# plot_coverage_figures reads coverage_stats.tsv (the composition figures) +
# gene_result parquet (per-gene DR/DL histograms), so the figures always match
# the numbers in the stats table. Editing the stats rule therefore rebuilds the
# figures too — a deliberate coupling for figure/table agreement.

_COVWORK = "results/2a_coverage/{dataset}/_work"


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
        annotation_reference=_ANNOT_REF,
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
            --annotation-reference {input.annotation_reference} \
            --output-annotations {output.annotations} \
            --output-gene-result {output.gene_result} &> {log}
        """


rule compute_coverage_stats:
    input:
        annotations=f"{_COVWORK}/annotations.parquet",
        gene_result=f"{_COVWORK}/gene_result.parquet",
    output:
        stats="results/2a_coverage/{dataset}/coverage_stats.tsv",
        detailed_genes_xlsx="results/2a_coverage/{dataset}/detailed_genes.xlsx",
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
            --output-stats {output.stats} \
            --output-detailed-genes-xlsx {output.detailed_genes_xlsx} &> {log}
        """



# -----------------------------------------------------------------------------
# Stage 2b: Plot figures
# -----------------------------------------------------------------------------
# Reads coverage_stats.tsv + gene_result.parquet -> ten cnsplots figures: one
# coverage composition + a DR and a DL distribution figure for EACH annotation
# dimension in config-visible DIMENSION_LABELS (characterisation_status,
# FYPOviability, deletion_essentiality), plus insertion placement per chromosome.
# The file names are derived from the column names in the script, so the set
# below is just that derivation spelled out for Snakemake.
#
# ponytail: 10 paths here + 10 derived in the script is the one duplication; they
# cannot drift silently (Snakemake fails on a missing declared output).


rule plot_coverage_figures:
    input:
        stats="results/2a_coverage/{dataset}/coverage_stats.tsv",
        gene_result=f"{_COVWORK}/gene_result.parquet",
    output:
        composition_characterisation="results/2a_coverage/{dataset}/coverage_by_characterisation_status.pdf",
        composition_fypoviability="results/2a_coverage/{dataset}/coverage_by_FYPOviability.pdf",
        composition_essentiality="results/2a_coverage/{dataset}/coverage_by_deletion_essentiality.pdf",
        insertion_placement="results/2a_coverage/{dataset}/coverage_insertion_placement.pdf",
        dr_characterisation="results/2a_coverage/{dataset}/coverage_dr_by_characterisation_status.pdf",
        dl_characterisation="results/2a_coverage/{dataset}/coverage_dl_by_characterisation_status.pdf",
        dr_fypoviability="results/2a_coverage/{dataset}/coverage_dr_by_FYPOviability.pdf",
        dl_fypoviability="results/2a_coverage/{dataset}/coverage_dl_by_FYPOviability.pdf",
        dr_essentiality="results/2a_coverage/{dataset}/coverage_dr_by_deletion_essentiality.pdf",
        dl_essentiality="results/2a_coverage/{dataset}/coverage_dl_by_deletion_essentiality.pdf",
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
            --output-dir results/2a_coverage/{wildcards.dataset} &> {log}
        """
