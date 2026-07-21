# =============================================================================
# clustering.smk — Deterministic gene-level clustering (finalize variants)
# =============================================================================
#
# Per-dataset: clusters genes in the 2-D depletion feature space (DR, DL) from
# gene-level curve fitting, into final_n_clusters groups. There is NO fixed
# "candidate" stage — each configured VARIANT clusters for itself (design doc
# 2026-07-21-clustering-finalize-variants).
#
# Spine (shared by all variants):
#   prepare_clustering_data -> annotated + scaled (DR, DL) matrix (+ k-sweep diagnostic)
#
# Then per variant (config.clustering.variants), keyed by `type`:
#   direct      : cluster_labels (k=final_n_clusters) -> finalize_direct (renumber by DR)
#   auto_merge  : cluster_labels (k=n_intermediate)   -> finalize_auto_merge (ward-merge centroids -> final_n_clusters)
#   grid        : finalize_grid (axis-cut grid on scaled DR/DL, no clustering step)
#   manual_merge: NO rule; notebooks/clustering/finalize_gene_clusters.ipynb writes
#                 the curated per-variant tsv (human-judgment merge).
# Every buildable variant emits final_clusters.tsv (unified `cluster` column,
# 1..final_n_clusters, WT last) + a per-variant metrics.tsv (silhouette/CH/DB) so
# variants can be compared. auto_merge/manual_merge also keep raw_cluster.
# enrichment fans out over every variant; ml.smk uses only selected_variant.
#
# Per-variant intermediate labels are pickles under the variant dir so label dtype
# survives round-trip; only final_clusters.tsv / metrics.tsv are user-facing.

_CLUSTER_METHODS = ["kmeans", "hierarchical_agg", "hierarchical_div", "gmm"]
_CWORK = "results/clustering/{dataset}/_work"

_CLUSTERING = config.get("clustering", {})
_VARIANTS = _CLUSTERING.get("variants", {})
_FINAL_N = _CLUSTERING.get("final_n_clusters", 9)


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


def buildable_variants() -> list[str]:
    """Variant names Snakemake can build end-to-end (every type except manual_merge)."""
    return [n for n, s in _VARIANTS.items() if s.get("type") in ("direct", "auto_merge", "grid")]


def all_variant_metrics(dataset: str) -> list[str]:
    """Per-variant metrics.tsv for every buildable variant (inputs to compare_variants)."""
    return [f"results/clustering/final/{dataset}/{v}/metrics.tsv" for v in buildable_variants()]


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
        random_state=_CLUSTERING.get("random_state", 42),
        k_min=_CLUSTERING.get("k_min", 2),
        k_max=_CLUSTERING.get("k_max", 20),
        dr_cap=_CLUSTERING.get("dr_cap", 1.3),
        dl_divisor=_CLUSTERING.get("dl_divisor", 10),
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


# --- Per-variant clustering: labels at n_intermediate (merge) or final_n_clusters (direct) ---
# Fanned out by the `variant` wildcard; used by direct + auto_merge finalize rules.
rule cluster_variant_labels:
    input:
        scaled=f"{_CWORK}/scaled_data.pkl",
    output:
        labels="results/clustering/final/{dataset}/{variant}/_labels.pkl",
    params:
        method=lambda wc: _VARIANTS[wc.variant].get("method", "kmeans"),
        final_n_clusters=_FINAL_N,
        # merge variants declare n_intermediate; direct variants omit it (cluster to final).
        n_intermediate_flag=lambda wc: (
            f"--n-intermediate {_VARIANTS[wc.variant]['n_intermediate']}"
            if "n_intermediate" in _VARIANTS[wc.variant] else ""
        ),
        random_state=_CLUSTERING.get("random_state", 42),
    wildcard_constraints:
        variant=_alt(_variants_of_type("direct") + _variants_of_type("auto_merge")),
    log:
        "logs/clustering/cluster_variant_labels_{dataset}_{variant}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [clustering] Clustering labels for {wildcards.variant} ({wildcards.dataset})..."
    shell:
        """
        python workflow/scripts/clustering/cluster_one_method.py \
            --method {params.method} \
            --scaled-data {input.scaled} \
            --output-labels {output.labels} \
            --final-n-clusters {params.final_n_clusters} \
            {params.n_intermediate_flag} \
            --random-state {params.random_state} &> {log}
        """


# --- Finalize: `direct` variant (labels already at final_n_clusters -> renumber by DR) ---
rule finalize_direct:
    input:
        annotated=f"{_CWORK}/annotated_data.pkl",
        scaled=f"{_CWORK}/scaled_data.pkl",
        labels="results/clustering/final/{dataset}/{variant}/_labels.pkl",
    output:
        clusters="results/clustering/final/{dataset}/{variant}/final_clusters.tsv",
        metrics="results/clustering/final/{dataset}/{variant}/metrics.tsv",
    params:
        n_clusters=_FINAL_N,
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
            --labels {input.labels} \
            --output {output.clusters} \
            --metrics-output {output.metrics} \
            --n-clusters {params.n_clusters} \
            --wt-cluster {params.wt_cluster} &> {log}
        """


# --- Finalize: `auto_merge` variant (ward-merge n_intermediate centroids -> final_n_clusters) ---
rule finalize_auto_merge:
    input:
        annotated=f"{_CWORK}/annotated_data.pkl",
        scaled=f"{_CWORK}/scaled_data.pkl",
        labels="results/clustering/final/{dataset}/{variant}/_labels.pkl",
    output:
        clusters="results/clustering/final/{dataset}/{variant}/final_clusters.tsv",
        metrics="results/clustering/final/{dataset}/{variant}/metrics.tsv",
    params:
        n_clusters=_FINAL_N,
        wt_cluster=config.get("enrichment", {}).get("wt_cluster", 9),
    wildcard_constraints:
        variant=_alt(_variants_of_type("auto_merge")),
    log:
        "logs/clustering/finalize_auto_merge_{dataset}_{variant}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [clustering] finalize_auto_merge {wildcards.variant} for {wildcards.dataset} (-> k={params.n_clusters})..."
    shell:
        """
        python workflow/scripts/clustering/finalize_auto_merge_clusters.py \
            --annotated-data {input.annotated} \
            --scaled-data {input.scaled} \
            --labels {input.labels} \
            --output {output.clusters} \
            --metrics-output {output.metrics} \
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
        metrics="results/clustering/final/{dataset}/{variant}/metrics.tsv",
    params:
        dr_cuts=lambda wc: " ".join(str(c) for c in _VARIANTS[wc.variant].get("dr_cuts", [])),
        dl_cuts=lambda wc: " ".join(str(c) for c in _VARIANTS[wc.variant].get("dl_cuts", [])),
        n_clusters=_FINAL_N,
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
            --metrics-output {output.metrics} \
            --dr-cuts {params.dr_cuts} \
            --dl-cuts {params.dl_cuts} \
            --n-clusters {params.n_clusters} \
            --wt-cluster {params.wt_cluster} &> {log}
        """


# --- Compare variants: gather every buildable variant's metrics into one table ---
# NOT part of `rule all` (ml only needs selected_variant). Request it explicitly to
# build ALL buildable variants at once and get a side-by-side metrics table for
# choosing which variant to select:
#   snakemake --use-conda results/clustering/final/{dataset}/variant_metrics_comparison.tsv
rule compare_variants:
    input:
        metrics=lambda wc: all_variant_metrics(wc.dataset),
    output:
        table="results/clustering/final/{dataset}/variant_metrics_comparison.tsv",
    params:
        variants=lambda wc: buildable_variants(),
    log:
        "logs/clustering/compare_variants_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [clustering] Comparing {wildcards.dataset} variants: {params.variants}"
    run:
        import pandas as pd
        rows = []
        for variant, path in zip(params.variants, input.metrics):
            df = pd.read_csv(path, sep="\t")
            df.insert(0, "variant", variant)
            rows.append(df)
        out = pd.concat(rows, ignore_index=True)
        out.to_csv(output.table, sep="\t", index=False)
