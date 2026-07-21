# =============================================================================
# utr.smk — UTR insertion classification (deterministic part of Batch C)
# =============================================================================
#
# Per-dataset: classifies intergenic insertions near gene boundaries as 5UTR or
# 3UTR (strand-aware, distance_threshold=400bp). The release/ annotations model
# each intergenic interval by its two flanking genes (pipe-separated
# Name/Systematic ID/Strand_Interval), so the script resolves each interval to a
# parental gene before the strand-aware call. Merges with insertion- + gene-level
# fitting stats, computes um_ratio (insertion DR / gene DR) and A_ratio.
#
# Single rule (no prepare/compute split — same shape as coverage.smk/verification.smk).
# The human-review notebook is notebooks/domain_analysis/review_utr_insertions.ipynb
# (Task 9), which reads this rule's utr_insertion_stats.tsv.

rule classify_utr_insertions:
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
        stats="results/utr/{dataset}/utr_insertion_stats.tsv",
    params:
        distance_threshold=400,
    log:
        "logs/utr/classify_utr_insertions_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [utr] Classifying UTR insertions for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/utr/classify_utr_insertions.py \
            --fitting-results {input.fitting_results} \
            --annotations {input.annotations} \
            --gene-level {input.gene_level} \
            --distance-threshold {params.distance_threshold} \
            --output-stats {output.stats} &> {log}
        """
