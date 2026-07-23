# =============================================================================
# noncoding_rna.smk — Non-coding RNA depletion analysis
# =============================================================================
#
# Split into 3 rules so the process and results are separated for
# readability/maintainability, same shape as verification.smk:
#   prepare_ncrna_table  -> combined / nuclear_trnas parquet intermediates
#                           (the single fan-out point)
#   compute_ncrna_stats  -> per-nuclear-tRNA stats TSV
#   plot_ncrna_figures   -> Feature-type donut + tRNA copy-number PDF
# The two downstream rules depend only on prepare_ncrna_table's output, so
# editing e.g. the figures never forces the stats TSV to rebuild. Ported from
# non_coding_RNA_analysis.ipynb.
#
# DATA-SOURCE NOTE (deviation from the plan skeleton — see Task 4 report):
# non-coding-gene curve fitting is NOT part of the DIT_HAP_snakemake release
# contract. The release only fits the 4,513 coding genes (zero tRNA rows), and
# no non-coding fitting exists anywhere under snakemake_repo (neither release/
# nor results/). The only real non-coding fitting is the legacy DIT_HAP_pipeline
# output below, which currently exists for HD_DIT_HAP only and still ships the
# pre-rename um/lam metric headers (normalized to DR/DL in the script). Until a
# native non-coding fitting artifact is added to the snakemake pipeline, the
# ncrna_fitting input is sourced per-dataset from this legacy map; datasets not
# in the map have no non-coding fitting and cannot be targeted.
_NONCODING_FITTING = {
    "HD_DIT_HAP": (
        "/data/c/yangyusheng_optimized/DIT_HAP_pipeline/results/"
        "HD_DIT_HAP_generationRAW/19_insertion_in_non_coding_genes/"
        "Non_coding_genes_Gene_level_statistics_fitted.tsv"
    ),
}

# Parquet intermediates shared by the two downstream rules.
_NCWORK = "results/noncoding_rna/{dataset}/_work"


rule prepare_ncrna_table:
    input:
        ncrna_fitting=lambda wc: _NONCODING_FITTING[wc.dataset],
        ncrna_bed=lambda wc: (
            f"resources/external/pombase/{DATASETS['reference']['pombase_version']}/"
            f"genome_region/non_coding_rna.bed"
        ),
        gtrnadb_bed="resources/external/pombase/schiPomb_972H-tRNAs.bed",
        marguerat_excel="resources/literature/margueratQuantitativeAnalysisFission2012.xlsx",
    output:
        combined=f"{_NCWORK}/combined.parquet",
        nuclear_trnas=f"{_NCWORK}/nuclear_trnas.parquet",
    log:
        "logs/noncoding_rna/prepare_ncrna_table_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [noncoding_rna] Preparing ncRNA tables for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/noncoding_rna/prepare_ncrna_table.py \
            --ncrna-fitting {input.ncrna_fitting} \
            --ncrna-bed {input.ncrna_bed} \
            --gtrnadb-bed {input.gtrnadb_bed} \
            --marguerat-excel {input.marguerat_excel} \
            --output-combined {output.combined} \
            --output-nuclear-trnas {output.nuclear_trnas} &> {log}
        """


rule compute_ncrna_stats:
    input:
        nuclear_trnas=f"{_NCWORK}/nuclear_trnas.parquet",
    output:
        stats="results/noncoding_rna/{dataset}/ncrna_stats.tsv",
    log:
        "logs/noncoding_rna/compute_ncrna_stats_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [noncoding_rna] Computing ncRNA stats for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/noncoding_rna/compute_ncrna_stats.py \
            --nuclear-trnas {input.nuclear_trnas} \
            --output-stats {output.stats} &> {log}
        """


rule plot_ncrna_figures:
    input:
        combined=f"{_NCWORK}/combined.parquet",
        nuclear_trnas=f"{_NCWORK}/nuclear_trnas.parquet",
    output:
        figures="results/noncoding_rna/{dataset}/ncrna_analysis.pdf",
    log:
        "logs/noncoding_rna/plot_ncrna_figures_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [noncoding_rna] Plotting ncRNA figures for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/noncoding_rna/plot_ncrna_figures.py \
            --combined {input.combined} \
            --nuclear-trnas {input.nuclear_trnas} \
            --output-figures {output.figures} &> {log}
        """
