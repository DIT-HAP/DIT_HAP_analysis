# =============================================================================
# domain_differences.smk — intra-gene DR heterogeneity candidates (Batch C)
# =============================================================================
#
# Per-dataset: for genes with high gene-level DR (> 0.15), positions each in-gene
# insertion along the CDS (insertion_fraction = Distance_to_start_codon /
# (Distance_to_start_codon + Distance_to_stop_codon), clamped to [0,1]) and
# reports per-gene distribution statistics (n_insertions, mean/std of
# insertion_fraction, gene_DR), sorted by std descending to surface genes whose
# insertions have the most heterogeneous positional spread — a proxy for
# functional sub-gene domains.
#
# Deterministic distillation of the visualization notebook
# DIT_HAP_pipeline/workflow/notebooks/genes_with_domain_differences.ipynb (the
# notebook itself only scatter-plots insertions + relies on human-curated
# spreadsheets; see the script docstring for the notebook-vs-script deviation).
#
# Single rule (no prepare/compute split — same shape as utr.smk/coverage.smk).

rule compute_domain_stats:
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
        stats="results/domain_differences/{dataset}/domain_candidate_stats.tsv",
    params:
        dr_threshold=0.15,
    log:
        "logs/domain_differences/compute_domain_stats_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [domain_differences] Computing domain candidate stats for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/domain_differences/compute_domain_stats.py \
            --fitting-results {input.fitting_results} \
            --annotations {input.annotations} \
            --gene-level {input.gene_level} \
            --dr-threshold {params.dr_threshold} \
            --output-stats {output.stats} &> {log}
        """
