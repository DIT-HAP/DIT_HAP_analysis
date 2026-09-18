# =============================================================================
# 2b_verification.smk — Deletion library phenotype verification
# =============================================================================
#
# Split into 4 rules so each analysis step is independently re-runnable:
#   prepare_verification_table  -> merged / verification parquet intermediates
#                                  (the single fan-out point)
#   verification_category_summary -> stats TSV + deletion-library comparison figure
#   verification_boxplots         -> category box plot + one figure and one review
#                                    TSV per critical-gene group
#   verification_depletion_curves -> one single-page per-gene depletion-curve
#                                    figure per critical-gene group
# The three figure rules depend only on the prepared parquets, so editing e.g.
# the boxplots never forces the depletion curves to rebuild. Ported from
# compare_with_deletion_library.ipynb; the altair charts stay notebook-only.
#
# Figures follow the house cnsplots pipeline (ADR-0001): workflow/src/figure_render/
# verification.py renders, workflow/src/figures.py owns the style, and every rule
# writes <stem>.pdf plus a <stem>.review.png preview. Only the PDFs are declared
# as outputs — the PNG is a review artifact.
#
# prepare_verification_table stays on the statistics env: it reads the curated
# .xlsx and renders nothing (the cnsplots env has no openpyxl).

# gRNA per-timepoint LFC (depletion-curve overlay) is HD-only. The curated
# fitted-parameters table is the project's single gRNA source (1c_annotate.smk and
# comparison.smk read the same file) — a per-dataset map, same pattern as
# noncoding_rna.smk's _NONCODING_FITTING. Datasets absent from the map render
# DIT-HAP-only curves (the --grna-timepoints flag is omitted).
_GRNA_PARAMETERS = "resources/curated/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv"
_GRNA_TIMEPOINT_DATA = {
    "HD_DIT_HAP": _GRNA_PARAMETERS,
}

# Parquet intermediates shared by the three figure rules.
_VWORK = "results/2b_verification/{dataset}/_work"


rule prepare_verification_table:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        deletion_library="resources/curated/deletion_library_categories.xlsx",
        essentiality_verification="resources/curated/essentiality_verification.csv",
    output:
        merged=f"{_VWORK}/merged.parquet",
        verification=f"{_VWORK}/verification.parquet",
    log:
        "logs/2b_verification/prepare_verification_table_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [verification] Preparing merged tables for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/verification/prepare_verification_table.py \
            --fitting-results {input.fitting_results} \
            --deletion-library {input.deletion_library} \
            --essentiality-verification {input.essentiality_verification} \
            --output-merged {output.merged} \
            --output-verification {output.verification} &> {log}
        """


rule verification_category_summary:
    input:
        merged=f"{_VWORK}/merged.parquet",
        verification=f"{_VWORK}/verification.parquet",
    output:
        stats="results/2b_verification/{dataset}/verification_stats.tsv",
        figure="results/2b_verification/{dataset}/deletion_library_comparison.pdf",
    log:
        "logs/2b_verification/verification_category_summary_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [verification] Category summary for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/verification/verification_category_summary.py \
            --merged {input.merged} \
            --verification {input.verification} \
            --output-stats {output.stats} \
            --output-figure results/2b_verification/{wildcards.dataset}/deletion_library_comparison &> {log}
        """


rule verification_boxplots:
    input:
        merged=f"{_VWORK}/merged.parquet",
        verification=f"{_VWORK}/verification.parquet",
    output:
        category_boxplot="results/2b_verification/{dataset}/verification_category_boxplot.pdf",
        critical_genes_dir=directory("results/2b_verification/{dataset}/critical_genes"),
    log:
        "logs/2b_verification/verification_boxplots_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [verification] Box plots + critical-gene figures/TSVs for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/verification/verification_boxplots.py \
            --merged {input.merged} \
            --verification {input.verification} \
            --output-figure results/2b_verification/{wildcards.dataset}/verification_category_boxplot \
            --output-critical-genes-dir results/2b_verification/{wildcards.dataset}/critical_genes &> {log}
        """


rule verification_depletion_curves:
    input:
        merged=f"{_VWORK}/merged.parquet",
        gene_timepoints=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/gene_level_fitting_statistics.tsv"
        ),
    output:
        depletion_curves=directory("results/2b_verification/{dataset}/depletion_curves"),
    params:
        # Optional gRNA overlay: build the flag only when the dataset is in the map.
        grna_flag=lambda wc: (
            f"--grna-timepoints {_GRNA_TIMEPOINT_DATA[wc.dataset]}"
            if wc.dataset in _GRNA_TIMEPOINT_DATA else ""
        ),
    log:
        "logs/2b_verification/verification_depletion_curves_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [verification] Depletion curves for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/verification/verification_depletion_curves.py \
            --merged {input.merged} \
            --gene-timepoints {input.gene_timepoints} \
            {params.grna_flag} \
            --output-dir results/2b_verification/{wildcards.dataset}/depletion_curves &> {log}
        """
