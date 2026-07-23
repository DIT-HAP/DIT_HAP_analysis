# =============================================================================
# verification.smk — Deletion library phenotype verification
# =============================================================================
#
# Split into 4 rules so each analysis step is independently re-runnable:
#   prepare_verification_table  -> merged / final_merged / simplified_verification
#                                  parquet intermediates (the single fan-out point)
#   verification_category_summary -> stats TSV + category donut/scatter PDF
#   verification_boxplots         -> boxplot/violin PDF + per-critical-group TSVs
#   verification_depletion_curves -> DIT-HAP (+gRNA) depletion-curve PDF
# The three figure rules depend only on the prepared parquets, so editing e.g.
# the boxplots never forces the depletion curves to rebuild. Ported from
# compare_with_deletion_library.ipynb; altair charts stay notebook-only.
#
# gRNA per-timepoint LFC (depletion-curve overlay) is HD-only and lives in the
# legacy pipeline repo, so it's sourced from a per-dataset map — same pattern as
# noncoding_rna.smk's _NONCODING_FITTING. Datasets absent from the map render
# DIT-HAP-only curves (the --grna-timepoints flag is omitted).
_GRNA_TIMEPOINT_DATA = {
    "HD_DIT_HAP": "/data/c/yangyusheng_optimized/DIT_HAP_pipeline/resources/HD_gRNA_data.csv",
}

# Parquet intermediates shared by the three figure rules.
_VWORK = "results/verification/{dataset}/_work"


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
        final_merged=f"{_VWORK}/final_merged.parquet",
        simplified_verification=f"{_VWORK}/simplified_verification.parquet",
    log:
        "logs/verification/prepare_verification_table_{dataset}.log",
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
            --output-final-merged {output.final_merged} \
            --output-simplified-verification {output.simplified_verification} &> {log}
        """


rule verification_category_summary:
    input:
        merged=f"{_VWORK}/merged.parquet",
        simplified_verification=f"{_VWORK}/simplified_verification.parquet",
    output:
        stats="results/verification/{dataset}/verification_stats.tsv",
        figures="results/verification/{dataset}/deletion_library_comparison.pdf",
    log:
        "logs/verification/verification_category_summary_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [verification] Category summary for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/verification/verification_category_summary.py \
            --merged {input.merged} \
            --simplified-verification {input.simplified_verification} \
            --output-stats {output.stats} \
            --output-figures {output.figures} &> {log}
        """


rule verification_boxplots:
    input:
        merged=f"{_VWORK}/merged.parquet",
        final_merged=f"{_VWORK}/final_merged.parquet",
        simplified_verification=f"{_VWORK}/simplified_verification.parquet",
    output:
        boxplots="results/verification/{dataset}/verification_boxplots.pdf",
        critical_genes_dir=directory("results/verification/{dataset}/critical_genes"),
    log:
        "logs/verification/verification_boxplots_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [verification] Boxplots + critical-gene TSVs for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/verification/verification_boxplots.py \
            --merged {input.merged} \
            --final-merged {input.final_merged} \
            --simplified-verification {input.simplified_verification} \
            --output-boxplots {output.boxplots} \
            --output-critical-genes-dir {output.critical_genes_dir} &> {log}
        """


rule verification_depletion_curves:
    input:
        merged=f"{_VWORK}/merged.parquet",
        gene_timepoints=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/gene_level_fitting_statistics.tsv"
        ),
    output:
        depletion_curves="results/verification/{dataset}/verification_depletion_curves.pdf",
    params:
        # Optional gRNA overlay: build the flag only when the dataset is in the map.
        grna_flag=lambda wc: (
            f"--grna-timepoints {_GRNA_TIMEPOINT_DATA[wc.dataset]}"
            if wc.dataset in _GRNA_TIMEPOINT_DATA else ""
        ),
    log:
        "logs/verification/verification_depletion_curves_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [verification] Depletion curves for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/verification/verification_depletion_curves.py \
            --merged {input.merged} \
            --gene-timepoints {input.gene_timepoints} \
            {params.grna_flag} \
            --output-depletion-curves {output.depletion_curves} &> {log}
        """
