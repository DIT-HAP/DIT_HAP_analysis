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
#   manual_merge: finalize_manual_merge runs notebooks/clustering/finalize_gene_clusters.ipynb
#                 headlessly (Snakemake native notebook: directive). The notebook holds a
#                 hand-tuned merge_groups dict (the human judgment), but is deterministic,
#                 so it builds like any other variant -> the same results/ final_clusters.tsv.
# Every buildable variant emits final_clusters.tsv (unified `cluster` column,
# 1..final_n_clusters, WT last) + a per-variant metrics.tsv (silhouette/CH/DB) +
# a cluster_scatter.pdf (DR/DL scatter, colored by cluster; a 2nd page for the
# pre-merge raw_cluster when the variant has one) so variants can be compared
# visually and numerically. auto_merge/manual_merge also keep raw_cluster.
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


_FINALIZE_TYPES = ("direct", "auto_merge", "grid", "manual_merge")


def final_clusters_path(dataset: str, variant: str) -> str:
    """Return the final_clusters.tsv path for one variant (all variant types build under results/)."""
    spec = _VARIANTS.get(variant)
    if spec is None:
        raise ValueError(f"Unknown clustering variant {variant!r} (configure it under clustering.variants)")
    vtype = spec.get("type")
    if vtype in _FINALIZE_TYPES:
        return f"results/clustering/{dataset}/{variant}/final_clusters.tsv"
    raise ValueError(f"Unknown finalize type {vtype!r} for variant {variant!r} (expected {'|'.join(_FINALIZE_TYPES)})")


def buildable_variants() -> list[str]:
    """Variant names Snakemake can build end-to-end.

    Every configured type is buildable: direct/auto_merge/grid via deterministic
    finalize scripts, manual_merge by executing finalize_gene_clusters.ipynb through
    the finalize_manual_merge rule (Snakemake's native notebook: directive). The
    notebook still carries the human-tuned merge_groups dict — that judgment is frozen
    into the DAG, but re-running reproduces it deterministically (random_state pinned).
    """
    return [n for n, s in _VARIANTS.items() if s.get("type") in _FINALIZE_TYPES]


def all_variant_metrics(dataset: str) -> list[str]:
    """Per-variant metrics.tsv for every buildable variant (inputs to compare_variants)."""
    return [f"results/clustering/{dataset}/{v}/metrics.tsv" for v in buildable_variants()]


def all_variant_scatters(dataset: str) -> list[str]:
    """Per-variant cluster_scatter.pdf for every buildable variant (request-all convenience target)."""
    return [f"results/clustering/{dataset}/{v}/cluster_scatter.pdf" for v in buildable_variants()]


def all_variant_final_clusters(dataset: str) -> list[str]:
    """Every buildable variant's final_clusters.tsv (inputs to the summary grid scatter)."""
    return [f"results/clustering/{dataset}/{v}/final_clusters.tsv" for v in buildable_variants()]


# --- Preprocessing spine (load/annotate/scale/k-sweep) ---
rule prepare_clustering_data:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        essentiality_verification_csv="resources/curated/essentiality_verification.csv",
    output:
        annotated=f"{_CWORK}/annotated_data.parquet",
        scaled=f"{_CWORK}/scaled_data.parquet",
        ksweep=f"{_CWORK}/k_sweep_metrics.parquet",
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
        scaled=f"{_CWORK}/scaled_data.parquet",
    output:
        labels="results/clustering/{dataset}/{variant}/_labels.parquet",
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
        annotated=f"{_CWORK}/annotated_data.parquet",
        scaled=f"{_CWORK}/scaled_data.parquet",
        labels="results/clustering/{dataset}/{variant}/_labels.parquet",
    output:
        clusters="results/clustering/{dataset}/{variant}/final_clusters.tsv",
        metrics="results/clustering/{dataset}/{variant}/metrics.tsv",
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
        annotated=f"{_CWORK}/annotated_data.parquet",
        scaled=f"{_CWORK}/scaled_data.parquet",
        labels="results/clustering/{dataset}/{variant}/_labels.parquet",
    output:
        clusters="results/clustering/{dataset}/{variant}/final_clusters.tsv",
        metrics="results/clustering/{dataset}/{variant}/metrics.tsv",
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
        annotated=f"{_CWORK}/annotated_data.parquet",
        scaled=f"{_CWORK}/scaled_data.parquet",
    output:
        clusters="results/clustering/{dataset}/{variant}/final_clusters.tsv",
        metrics="results/clustering/{dataset}/{variant}/metrics.tsv",
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


# --- Finalize: `manual_merge` variant (execute the human-judgment notebook headlessly) ---
# Unlike the other finalize rules (shell -> workflow/scripts/*.py), this one uses
# Snakemake's native notebook: directive to run finalize_gene_clusters.ipynb, which
# reads the injected `snakemake` object (dataset/variant/params/inputs/outputs) and
# writes both final_clusters.tsv and metrics.tsv itself. The notebook's merge_groups
# dict is the frozen human judgment; everything else is deterministic. The executed
# notebook is captured via log: notebook=... for review.
rule finalize_manual_merge:
    input:
        annotated=f"{_CWORK}/annotated_data.parquet",
        scaled=f"{_CWORK}/scaled_data.parquet",
    output:
        clusters="results/clustering/{dataset}/{variant}/final_clusters.tsv",
        metrics="results/clustering/{dataset}/{variant}/metrics.tsv",
    params:
        method=lambda wc: _VARIANTS[wc.variant].get("method", "kmeans"),
        n_intermediate=lambda wc: _VARIANTS[wc.variant].get("n_intermediate", 64),
        final_n_clusters=_FINAL_N,
        random_state=_CLUSTERING.get("random_state", 42),
        wt_cluster=config.get("enrichment", {}).get("wt_cluster", 9),
    wildcard_constraints:
        variant=_alt(_variants_of_type("manual_merge")),
    log:
        notebook="logs/clustering/finalize_manual_merge_{dataset}_{variant}.ipynb",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [clustering] finalize_manual_merge {wildcards.variant} for {wildcards.dataset} (executing notebook)..."
    notebook:
        "../../notebooks/clustering/finalize_gene_clusters.ipynb"


# --- Plot: DR/DL scatter of one variant's final clusters (+ intermediate, if any) ---
# Works for every configured variant (uses the global `variant` wildcard_constraints
# above). Every variant, manual_merge included, now builds its final_clusters.tsv
# under results/ via its finalize rule, so this just needs that tsv as input.
rule plot_variant_clusters:
    input:
        final_clusters=lambda wc: final_clusters_path(wc.dataset, wc.variant),
    output:
        scatter="results/clustering/{dataset}/{variant}/cluster_scatter.pdf",
    log:
        "logs/clustering/plot_variant_clusters_{dataset}_{variant}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [clustering] Plotting DR/DL cluster scatter for {wildcards.variant} ({wildcards.dataset})..."
    shell:
        """
        python workflow/scripts/clustering/plot_variant_clusters.py \
            --final-clusters {input.final_clusters} \
            --output {output.scatter} \
            --variant-label {wildcards.variant} &> {log}
        """


# --- Compare variants: gather every buildable variant's metrics into one table ---
# NOT part of `rule all` (ml only needs selected_variant). Request it explicitly to
# build ALL buildable variants at once and get a side-by-side metrics table for
# choosing which variant to select:
#   snakemake --use-conda results/clustering/{dataset}/variant_metrics_comparison.tsv
rule compare_variants:
    input:
        metrics=lambda wc: all_variant_metrics(wc.dataset),
    output:
        table="results/clustering/{dataset}/variant_metrics_comparison.tsv",
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


# --- Plot all: build every buildable variant's cluster_scatter.pdf in one go ---
# Also NOT part of `rule all`. Request it explicitly to render every variant's
# scatter for a visual side-by-side (pairs with compare_variants' metrics table):
#   snakemake --use-conda results/clustering/{dataset}/all_variants_plotted.done
rule plot_all_variants:
    input:
        scatters=lambda wc: all_variant_scatters(wc.dataset),
    output:
        marker=touch("results/clustering/{dataset}/all_variants_plotted.done"),
    message:
        "*** [clustering] Plotted all {wildcards.dataset} variants: {input.scatters}"


# --- Summary grid: every buildable variant's final clusters on one page ---
# Complements the per-variant cluster_scatter.pdf files with a single side-by-side
# grid (one subplot per variant, final `cluster` column, shared DR/DL axes) for a
# quick visual comparison of how the finalize strategies differ. NOT part of
# `rule all`; request explicitly:
#   snakemake --use-conda results/clustering/{dataset}/all_variants_cluster_scatter.pdf
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
        "../envs/statistics_and_figure_plotting.yml"
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
