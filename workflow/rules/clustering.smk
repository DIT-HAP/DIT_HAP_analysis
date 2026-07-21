# =============================================================================
# clustering.smk — Deterministic gene-level clustering (candidate labeling)
# =============================================================================
#
# Per-dataset: clusters genes in the 2-D depletion feature space (DR, DL) from
# gene-level curve fitting. Produces CANDIDATE labels (64), then finalizes to 9.
#
# Candidate stage, split by clustering method, mirroring ml.smk's fan-out:
#   prepare  -> scaled (DR, DL) matrix + k-sweep (the shared "spine")
#   cluster  -> one job per method (kmeans / hierarchical_agg / hierarchical_div / gmm)
#   select   -> attach the pinned best-method (kmeans) labels + aggregate metrics
# Per-method intermediates are pickles under _work/ so label dtype and exact
# metric precision survive round-trip; only the final two files are TSV.
#
# Finalize stage (64 -> 9) has two paths, selected by finalize_mode (design doc §2):
#   auto (default): the auto_finalize_clusters rule below reuses the prepare spine
#     to cluster to k=9 deterministically (lowest-DR cluster = WT = 9).
#   manual: notebooks/clustering/finalize_gene_clusters.ipynb writes the curated
#     resources/curated/final_clusters.tsv (human-judgment 64->9 merge).
# Both emit the unified `cluster` column consumed by enrichment.smk + ml.smk.

_CLUSTER_METHODS = ["kmeans", "hierarchical_agg", "hierarchical_div", "gmm"]
_CWORK = "results/clustering/candidates/{dataset}/_work"

wildcard_constraints:
    method="|".join(_CLUSTER_METHODS),


def final_clusters_path(dataset: str) -> str:
    """Return the final_clusters.tsv path per config.clustering.finalize_mode (+ per-dataset override)."""
    cl = config.get("clustering", {})
    mode = cl.get("finalize_mode_overrides", {}).get(dataset, cl.get("finalize_mode", "auto"))
    if mode == "manual":
        return "resources/curated/final_clusters.tsv"
    if mode == "auto":
        return f"results/clustering/final/{dataset}/final_clusters.tsv"
    raise ValueError(f"Unknown finalize_mode {mode!r} for dataset {dataset!r} (expected auto|manual)")


# --- Preprocessing spine (load/annotate/scale/k-sweep) ---
rule prepare_clustering_data:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        essentiality_verification_csv="resources/curated/essentiality_verification.csv",
    output:
        annotated=f"{_CWORK}/annotated_data.pkl",
        scaled=f"{_CWORK}/scaled_data.pkl",
        ksweep=f"{_CWORK}/k_sweep_metrics.pkl",
    params:
        random_state=config.get("clustering", {}).get("random_state", 42),
        k_min=config.get("clustering", {}).get("k_min", 2),
        k_max=config.get("clustering", {}).get("k_max", 20),
        dr_cap=config.get("clustering", {}).get("dr_cap", 1.3),
        dl_divisor=config.get("clustering", {}).get("dl_divisor", 10),
    log:
        "logs/clustering/prepare_clustering_data_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [clustering] Preparing scaled matrix + k-sweep for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/clustering/prepare_clustering_data.py \
            --fitting-results {input.fitting_results} \
            --essentiality-verification-csv {input.essentiality_verification_csv} \
            --output-annotated {output.annotated} \
            --output-scaled {output.scaled} \
            --output-ksweep {output.ksweep} \
            --random-state {params.random_state} \
            --k-min {params.k_min} \
            --k-max {params.k_max} \
            --dr-cap {params.dr_cap} \
            --dl-divisor {params.dl_divisor} &> {log}
        """


# --- One clustering method (fanned out by the `method` wildcard) ---
rule cluster_one_method:
    input:
        scaled=f"{_CWORK}/scaled_data.pkl",
    output:
        labels=f"{_CWORK}/{{method}}_labels.pkl",
        metrics=f"{_CWORK}/{{method}}_metrics.pkl",
    params:
        n_clusters=config.get("clustering", {}).get("n_clusters", 64),
        random_state=config.get("clustering", {}).get("random_state", 42),
    log:
        "logs/clustering/cluster_{method}_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [clustering] Running {wildcards.method} for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/clustering/cluster_one_method.py \
            --method {wildcards.method} \
            --scaled-data {input.scaled} \
            --output-labels {output.labels} \
            --output-metrics {output.metrics} \
            --n-clusters {params.n_clusters} \
            --random-state {params.random_state} &> {log}
        """


# --- Select best-method labels + aggregate metrics ---
rule select_candidate_clusters:
    input:
        annotated=f"{_CWORK}/annotated_data.pkl",
        ksweep=f"{_CWORK}/k_sweep_metrics.pkl",
        best_labels=f"{_CWORK}/kmeans_labels.pkl",
        method_metrics=expand(f"{_CWORK}/{{method}}_metrics.pkl", method=_CLUSTER_METHODS, allow_missing=True),
    output:
        clusters="results/clustering/candidates/{dataset}/candidate_clusters.tsv",
        metrics="results/clustering/candidates/{dataset}/clustering_metrics.tsv",
    log:
        "logs/clustering/select_candidate_clusters_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [clustering] Selecting candidate clusters for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/clustering/select_candidate_clusters.py \
            --annotated-data {input.annotated} \
            --ksweep {input.ksweep} \
            --best-labels {input.best_labels} \
            --method-metrics {input.method_metrics} \
            --output {output.clusters} \
            --metrics-output {output.metrics} &> {log}
        """


# --- Automatic finalize (deterministic alternative to the manual notebook) ---
# Reuses the prepare spine pickles; clusters to k=9 and DR-numbers (design doc §2-3).
rule auto_finalize_clusters:
    input:
        annotated=f"{_CWORK}/annotated_data.pkl",
        scaled=f"{_CWORK}/scaled_data.pkl",
    output:
        clusters="results/clustering/final/{dataset}/final_clusters.tsv",
    params:
        n_clusters=config.get("clustering", {}).get("final_n_clusters", 9),
        random_state=config.get("clustering", {}).get("random_state", 42),
        wt_cluster=config.get("enrichment", {}).get("wt_cluster", 9),
    log:
        "logs/clustering/auto_finalize_clusters_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [clustering] Auto-finalizing clusters for {wildcards.dataset} (k={params.n_clusters})..."
    shell:
        """
        python workflow/scripts/clustering/auto_finalize_clusters.py \
            --annotated-data {input.annotated} \
            --scaled-data {input.scaled} \
            --output {output.clusters} \
            --n-clusters {params.n_clusters} \
            --random-state {params.random_state} \
            --wt-cluster {params.wt_cluster} &> {log}
        """
