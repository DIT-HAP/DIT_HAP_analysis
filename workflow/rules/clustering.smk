# =============================================================================
# clustering.smk — Deterministic gene-level clustering (candidate labeling)
# =============================================================================

# Per-dataset: clusters genes in the 2-D depletion feature space (um, lam) from
# gene-level curve fitting. Produces CANDIDATE labels only — the manual 64->9
# merge lives in notebooks/clustering/finalize_gene_clusters.ipynb, which writes
# resources/curated/final_clusters.tsv (design doc §5).
rule generate_candidate_clusters:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        essentiality_verification_csv="resources/curated/essentiality_verification.csv",
    output:
        clusters="results/clustering/candidates/{dataset}/candidate_clusters.tsv",
        metrics="results/clustering/candidates/{dataset}/clustering_metrics.tsv",
    params:
        n_clusters=config.get("clustering", {}).get("n_clusters", 64),
        random_state=config.get("clustering", {}).get("random_state", 42),
    log:
        "logs/clustering/generate_candidate_clusters_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** Generating candidate clusters for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/clustering/generate_candidate_clusters.py \
            --fitting-results {input.fitting_results} \
            --essentiality-verification-csv {input.essentiality_verification_csv} \
            --output {output.clusters} \
            --metrics-output {output.metrics} \
            --n-clusters {params.n_clusters} \
            --random-state {params.random_state} &> {log}
        """
