# =============================================================================
# noncoding_rna.smk — Non-coding RNA depletion analysis
# =============================================================================
#
# Per-dataset: merges non-coding-gene depletion stats with GtRNAdb tRNA
# annotations (matched by chr+start+end, NOT by name) and Marguerat 2012 mRNA
# abundance, then characterizes nuclear tRNA depletion (copy number,
# amino-acid/anticodon, DR). Single rule (no prepare/compute split — data is
# tiny and self-contained, same shape as coverage.smk / verification.smk).
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


rule analyze_noncoding_rna:
    input:
        ncrna_fitting=lambda wc: _NONCODING_FITTING[wc.dataset],
        ncrna_bed=lambda wc: (
            f"resources/external/pombase/{DATASETS['reference']['pombase_version']}/"
            f"genome_region/non_coding_rna.bed"
        ),
        gtrnadb_bed="resources/external/pombase/schiPomb_972H-tRNAs.bed",
        marguerat_excel="resources/literature/margueratQuantitativeAnalysisFission2012.xlsx",
    output:
        stats="results/noncoding_rna/{dataset}/ncrna_stats.tsv",
        figures="results/noncoding_rna/{dataset}/ncrna_analysis.pdf",
    log:
        "logs/noncoding_rna/analyze_noncoding_rna_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [noncoding_rna] Analyzing ncRNA depletion for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/noncoding_rna/analyze_noncoding_rna.py \
            --ncrna-fitting {input.ncrna_fitting} \
            --ncrna-bed {input.ncrna_bed} \
            --gtrnadb-bed {input.gtrnadb_bed} \
            --marguerat-excel {input.marguerat_excel} \
            --output-stats {output.stats} \
            --output-figures {output.figures} &> {log}
        """
