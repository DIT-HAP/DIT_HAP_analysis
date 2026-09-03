# =============================================================================
# figure.smk — Centralized plotting rules (computation-plotting decoupling)
# =============================================================================
#
# This file centralizes all figure-generation rules (rules that output PDF/PNG).
# Separation of concerns: computation rules (data generation) remain in their
# domain .smk files (clustering.smk, enrichment.smk, coherence.smk, etc.);
# plotting rules (consume data, output figures) are unified here.
#
# Naming conventions (aligned with upstream DIT_HAP_snakemake, ADR-0001):
#   Rule names:  plot_<feature> (e.g. plot_coherence, plot_variant_clusters)
#   Environment: conda: "../envs/cnsplots.yml"
#   Output paths: direct content description, no type suffix (coherence.pdf, not coherence_figure.pdf)
#
# Usage pattern:
#   rule plot_example:
#       conda: "../envs/cnsplots.yml"
#       input: "results/example/{dataset}/data.parquet"
#       output: "results/example/{dataset}/example.pdf"
#       script: "../scripts/example/plot_example.py"

# -----------------------------------------------------------------------------
# PCR Quality Control
# -----------------------------------------------------------------------------
_PCRWORK = "results/pcr_qc/_work"

rule plot_pcr_qc:
    input:
        pbl_pbr=f"{_PCRWORK}/pbl_pbr.parquet",
        tech=f"{_PCRWORK}/tech.parquet",
        bio=f"{_PCRWORK}/bio.parquet",
        spikein=f"{_PCRWORK}/spikein.parquet",
    output:
        "results/pcr_qc/PCR_quality_control.pdf",
    log:
        "logs/pcr_qc/plot_pcr_qc.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [pcr_qc] Building 2x2 library-prep QC figure..."
    shell:
        """
        python workflow/scripts/pcr_qc/plot_pcr_qc.py \
            --pbl-pbr {input.pbl_pbr} \
            --tech {input.tech} \
            --bio {input.bio} \
            --spikein {input.spikein} \
            --output {output} &> {log}
        """


# -----------------------------------------------------------------------------
# Clustering
# -----------------------------------------------------------------------------
rule plot_variant_clusters:
    input:
        final_clusters=lambda wc: final_clusters_path(wc.dataset, wc.variant),
    output:
        scatter="results/clustering/{dataset}/{variant}/cluster_scatter.pdf",
    log:
        "logs/clustering/plot_variant_clusters_{dataset}_{variant}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [clustering] Plotting DR/DL cluster scatter for {wildcards.variant} ({wildcards.dataset})..."
    shell:
        """
        python workflow/scripts/clustering/plot_variant_clusters.py \
            --final-clusters {input.final_clusters} \
            --output {output.scatter} \
            --variant-label {wildcards.variant} &> {log}
        """


rule plot_all_variants_grid:
    input:
        final_clusters=lambda wc: all_variant_final_clusters(wc.dataset),
    output:
        scatter="results/clustering/{dataset}/all_variants_cluster_scatter.pdf",
    params:
        variant_labels=lambda wc: buildable_variants(),
    log:
        "logs/clustering/plot_all_variants_grid_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [clustering] Plotting summary grid of all {wildcards.dataset} variants' final clusters..."
    shell:
        """
        python workflow/scripts/clustering/plot_all_variant_clusters.py \
            --final-clusters {input.final_clusters} \
            --variant-labels {params.variant_labels} \
            --dataset {wildcards.dataset} \
            --output {output.scatter} &> {log}
        """


# -----------------------------------------------------------------------------
# Comparison
# -----------------------------------------------------------------------------
_CWORK = "results/comparison/{dataset}/_work"

rule plot_comparison_figures:
    input:
        fitness_table=f"{_CWORK}/fitness_table.parquet",
        stats="results/comparison/{dataset}/fitness_correlation_stats.tsv",
    output:
        figures="results/comparison/{dataset}/pairwise_fitness_comparison.pdf",
    log:
        "logs/comparison/plot_comparison_figures_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [comparison] Plotting pairwise comparison figures for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/comparison/plot_comparison_figures.py \
            --fitness-table {input.fitness_table} \
            --stats {input.stats} \
            --output-figures {output.figures} &> {log}
        """


# -----------------------------------------------------------------------------
# Coverage
# -----------------------------------------------------------------------------
_COVWORK = "results/coverage/{dataset}/_work"

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


# -----------------------------------------------------------------------------
# Non-coding RNA
# -----------------------------------------------------------------------------
_NCWORK = "results/noncoding_rna/{dataset}/_work"

rule plot_ncrna_figures:
    input:
        combined=f"{_NCWORK}/combined.parquet",
        nuclear_trnas=f"{_NCWORK}/nuclear_trnas.parquet",
    output:
        figures="results/noncoding_rna/{dataset}/ncrna_analysis.pdf",
    log:
        "logs/noncoding_rna/plot_ncrna_figures_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [noncoding_rna] Plotting ncRNA figures for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/noncoding_rna/plot_ncrna_figures.py \
            --combined {input.combined} \
            --nuclear-trnas {input.nuclear_trnas} \
            --output-figures {output.figures} &> {log}
        """

