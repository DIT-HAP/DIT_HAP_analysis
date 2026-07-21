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
# Finalize stage (scaled/64 -> 9) is a set of named VARIANTS, one per configured
# strategy (design doc 2026-07-21-clustering-finalize-variants §2). Each variant
# declares a `type`; the three buildable types get a rule below:
#   direct       -> finalize_direct rule       (cluster fresh to k=9, DR-numbered)
#   auto_merge   -> finalize_auto_merge rule    (ward-merge k=64 centroids to 9)
#   grid         -> finalize_grid rule          (axis-cut grid, DR-numbered)
#   manual_merge -> NO rule; notebooks/clustering/finalize_gene_clusters.ipynb
#                   writes the curated per-variant tsv (human-judgment merge).
# All emit the unified `cluster` column consumed by enrichment.smk + ml.smk;
# auto_merge/manual_merge also keep raw_cluster. enrichment fans out over every
# variant (compare); ml.smk uses only selected_variant.

_CLUSTER_METHODS = ["kmeans", "hierarchical_agg", "hierarchical_div", "gmm"]
_CWORK = "results/clustering/candidates/{dataset}/_work"

_CLUSTERING = config.get("clustering", {})
_VARIANTS = _CLUSTERING.get("variants", {})


def _variants_of_type(t: str) -> list[str]:
    """Names of configured finalize variants of a given type."""
    return [name for name, spec in _VARIANTS.items() if spec.get("type") == t]


# Alternation patterns for per-rule variant constraints (never-match sentinel when
# a type has no configured variants, so an empty join can't match a real target).
def _alt(names: list[str]) -> str:
    return "|".join(names) if names else "a^"


wildcard_constraints:
    method="|".join(_CLUSTER_METHODS),
    variant=_alt(list(_VARIANTS)),


def _candidate_labels_pkl(wc) -> str:
    """The k=64 candidate-label pickle an auto_merge variant reuses (from cluster_one_method)."""
    method = _VARIANTS[wc.variant].get("method", "kmeans")
    return f"results/clustering/candidates/{wc.dataset}/_work/{method}_labels.pkl"


def selected_variant(dataset: str) -> str:
    """The single finalize variant ml.smk reads (config.selected_variant + per-dataset override)."""
    return _CLUSTERING.get("selected_variant_overrides", {}).get(
        dataset, _CLUSTERING.get("selected_variant")
    )


def final_clusters_path(dataset: str, variant: str) -> str:
    """Return the final_clusters.tsv path for one variant (buildable results/ vs curated manual)."""
    spec = _VARIANTS.get(variant)
    if spec is None:
        raise ValueError(f"Unknown clustering variant {variant!r} (configure it under clustering.variants)")
    vtype = spec.get("type")
    if vtype == "manual_merge":
        return f"resources/curated/final_clusters/{dataset}/{variant}.tsv"
    if vtype in ("direct", "auto_merge", "grid"):
        return f"results/clustering/final/{dataset}/{variant}/final_clusters.tsv"
    raise ValueError(f"Unknown finalize type {vtype!r} for variant {variant!r} (expected direct|auto_merge|grid|manual_merge)")


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


# --- Finalize: `direct` variant (cluster fresh to k=9, DR-numbered) ---
rule finalize_direct:
    input:
        annotated=f"{_CWORK}/annotated_data.pkl",
        scaled=f"{_CWORK}/scaled_data.pkl",
    output:
        clusters="results/clustering/final/{dataset}/{variant}/final_clusters.tsv",
    params:
        method=lambda wc: _VARIANTS[wc.variant].get("method", "kmeans"),
        n_clusters=_CLUSTERING.get("final_n_clusters", 9),
        random_state=_CLUSTERING.get("random_state", 42),
        wt_cluster=config.get("enrichment", {}).get("wt_cluster", 9),
    wildcard_constraints:
        variant=_alt(_variants_of_type("direct")),
    log:
        "logs/clustering/finalize_direct_{dataset}_{variant}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [clustering] finalize_direct {wildcards.variant} for {wildcards.dataset} (k={params.n_clusters})..."
    shell:
        """
        python workflow/scripts/clustering/finalize_direct_clusters.py \
            --annotated-data {input.annotated} \
            --scaled-data {input.scaled} \
            --output {output.clusters} \
            --method {params.method} \
            --n-clusters {params.n_clusters} \
            --random-state {params.random_state} \
            --wt-cluster {params.wt_cluster} &> {log}
        """


# --- Finalize: `auto_merge` variant (ward-merge the method's k=64 centroids to 9) ---
rule finalize_auto_merge:
    input:
        annotated=f"{_CWORK}/annotated_data.pkl",
        scaled=f"{_CWORK}/scaled_data.pkl",
        candidate_labels=_candidate_labels_pkl,
    output:
        clusters="results/clustering/final/{dataset}/{variant}/final_clusters.tsv",
    params:
        n_clusters=_CLUSTERING.get("final_n_clusters", 9),
        wt_cluster=config.get("enrichment", {}).get("wt_cluster", 9),
    wildcard_constraints:
        variant=_alt(_variants_of_type("auto_merge")),
    log:
        "logs/clustering/finalize_auto_merge_{dataset}_{variant}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [clustering] finalize_auto_merge {wildcards.variant} for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/clustering/finalize_auto_merge_clusters.py \
            --annotated-data {input.annotated} \
            --scaled-data {input.scaled} \
            --candidate-labels {input.candidate_labels} \
            --output {output.clusters} \
            --n-clusters {params.n_clusters} \
            --wt-cluster {params.wt_cluster} &> {log}
        """


# --- Finalize: `grid` variant (axis-cut grid on scaled DR/DL, DR-numbered) ---
rule finalize_grid:
    input:
        annotated=f"{_CWORK}/annotated_data.pkl",
        scaled=f"{_CWORK}/scaled_data.pkl",
    output:
        clusters="results/clustering/final/{dataset}/{variant}/final_clusters.tsv",
    params:
        dr_cuts=lambda wc: " ".join(str(c) for c in _VARIANTS[wc.variant].get("dr_cuts", [])),
        dl_cuts=lambda wc: " ".join(str(c) for c in _VARIANTS[wc.variant].get("dl_cuts", [])),
        n_clusters=_CLUSTERING.get("final_n_clusters", 9),
        wt_cluster=config.get("enrichment", {}).get("wt_cluster", 9),
    wildcard_constraints:
        variant=_alt(_variants_of_type("grid")),
    log:
        "logs/clustering/finalize_grid_{dataset}_{variant}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [clustering] finalize_grid {wildcards.variant} for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/clustering/finalize_grid_clusters.py \
            --annotated-data {input.annotated} \
            --scaled-data {input.scaled} \
            --output {output.clusters} \
            --dr-cuts {params.dr_cuts} \
            --dl-cuts {params.dl_cuts} \
            --n-clusters {params.n_clusters} \
            --wt-cluster {params.wt_cluster} &> {log}
        """
