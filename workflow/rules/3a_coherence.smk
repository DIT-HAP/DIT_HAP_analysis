# =============================================================================
# 3a_coherence.smk — Gene-group coherence in DR-DL fitness space
# =============================================================================
#
# `3a` is the <chapter><section> prefix convention (see the Includes block in
# Snakefile). It carries into results/3a_coherence/ and logs/3a_coherence/ so
# `ls results/` reads in analysis order; the script and library directories keep
# the unprefixed name (workflow/scripts/coherence/, workflow/src/coherence/).
#
# Per (dataset × source): for every group (GO complex / CC / BP term) whose
# DR<threshold members number min_group_size..max_group_size AND whose total
# annotated size is <= max_term_genes, measure how tightly they cluster in
# normalized (DR, DL/10) space vs a seeded permutation null (median pairwise
# distance z-score). Weiszfeld geometric median + seeded MPD permutation test
# live in workflow/src/coherence/metrics.py.
#
# Rules fanning out by {source}, then re-aggregated:
#   prepare_coherence_annotation   -> long-form group annotation (per source)
#   compute_coherence              -> coherence_metrics.parquet
#                                     ({method}_p = one-sided add-one permutation p;
#                                     q_value = per-source Benjamini-Hochberg FDR
#                                     over the primary method's p;
#                                     plus the shared-subunit fraction and, when
#                                     features_panels is on, the abundance /
#                                     conservation uniformity CVs — all computed
#                                     here so the figure only renders)
#   plot_coherence                 -> coherence.pdf (reads the metrics table alone)
#   combine_coherence_metrics      -> coherence_metrics_combined.parquet (all sources;
#                                     FDR re-derived over the union — a q is only
#                                     defined w.r.t. a family, and this table's family
#                                     is its own row set, so cross-source thresholds
#                                     and the dedup representative ranking share one
#                                     q scale. Per-source tables keep their own q.)
#   deduplicate_coherence_terms    -> coherence_terms_deduplicated.tsv (+ _representatives.tsv):
#                                     collapse redundant terms by member overlap + GO DAG
#                                     (display layer; the combined q_value is carried
#                                     through untouched, all terms kept)
#   compute_coherence_attribution  -> incoherence_attribution.tsv + incoherence_split_points.parquet:
#                                     diagnose WHY a complex is dispersed (major/minor GMM
#                                     split, shared subunits, paralog buffering).
#                                     attribution_sources subset.
#   plot_coherence_attribution     -> incoherence_attribution.pdf (renders the persisted split points)
#   plot_coherence_group_scatter   -> group_scatter.pdf (named groups from config)
#
# Data-format rule: per-stage intermediates are Parquet (exact dtypes survive
# round-trip); only the final human-facing artifacts (dedup + attribution tables)
# are TSV.
#
# Plot rules use workflow/src/figures.py's house style and save_dual(), which
# writes a .pdf plus a .review.png sibling — both are declared as outputs so the
# DAG is honest about what lands on disk and `--delete-all-output` cleans it.
#
# Sources are registered in config.coherence.sources (currently: go_macrocomplex,
# go_cc, go_bp). Each source has a loader in workflow/src/coherence/sources.py.
# scatter_groups drives plot_group_scatter (per-source namelist, entries match
# group_name OR group_id).
#
# DATA-PATH NOTE: fitting_results comes from upstream DIT_HAP_snakemake release/
# dirs via DATASETS['datasets'][dataset]['release_dir']; pombase from
# resources/external/pombase/{DATASETS['reference']['pombase_version']}.

import json

_COH = "results/3a_coherence/{dataset}/{source}"
_COH_CFG = config.get("coherence", {})
_COH_SOURCES = _COH_CFG.get("sources", ["go_macrocomplex", "go_cc", "go_bp"])
# Incoherence attribution runs on a (usually smaller) subset of sources: major/minor
# split + shared-subunit are only meaningful for physical complexes, not GO_BP processes.
_COH_ATTR_SOURCES = _COH_CFG.get("attribution_sources", ["go_macrocomplex"])
_COH_FEATURES_PANELS = _COH_CFG.get("features_panels", True)
_COH_FEATURES = (
    f"results/1b_features/{DATASETS['reference']['pombase_version']}"
    f"/pombe_coding_gene_protein_features.tsv"
)

wildcard_constraints:
    source="|".join(_COH_SOURCES),


rule prepare_coherence_annotation:
    input:
        pombase_dir=lambda wc: f"resources/external/pombase/{DATASETS['reference']['pombase_version']}",
    output:
        long_table=f"{_COH}/group_annotation_long.tsv",
    log:
        "logs/3a_coherence/prepare_{dataset}_{source}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [coherence] Preparing {wildcards.source} annotation for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coherence/prepare_annotation.py \
            --source {wildcards.source} \
            --pombase-dir {input.pombase_dir} \
            --output {output.long_table} &> {log}
        """


rule compute_coherence:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        annotation=f"{_COH}/group_annotation_long.tsv",
        features=lambda wc: _COH_FEATURES if _COH_FEATURES_PANELS else [],
    output:
        f"{_COH}/coherence_metrics.parquet",
    params:
        min_size=_COH_CFG.get("min_group_size", 3),
        max_size=_COH_CFG.get("max_group_size", 300),
        max_term_genes=_COH_CFG.get("max_term_genes", 500),
        # Fallback must match config/analysis.yaml and the script default: the DR
        # sign flipped on 2026-09-17 (negative = depleted), so a stale +0.3 here
        # would silently select the *opposite* gene set if the key ever went missing.
        dr_threshold=_COH_CFG.get("dr_threshold", -0.3),
        n_permutations=_COH_CFG.get("n_permutations", 1000),
        random_state=_COH_CFG.get("random_state", 42),
        features_flag=lambda wc, input: (
            f"--features {input.features}" if _COH_FEATURES_PANELS and input.features else ""
        ),
    log:
        "logs/3a_coherence/compute_{dataset}_{source}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coherence] Computing coherence metrics for {wildcards.dataset} × {wildcards.source}..."
    shell:
        """
        python workflow/scripts/coherence/compute_coherence.py \
            --fitting-results {input.fitting_results} \
            --annotation {input.annotation} \
            {params.features_flag} \
            --source {wildcards.source} \
            --min-size {params.min_size} \
            --max-size {params.max_size} \
            --max-term-genes {params.max_term_genes} \
            --dr-threshold {params.dr_threshold} \
            --n-permutations {params.n_permutations} \
            --random-state {params.random_state} \
            --output {output} &> {log}
        """


rule plot_coherence:
    input:
        metrics=f"{_COH}/coherence_metrics.parquet",
    output:
        figure=f"{_COH}/coherence.pdf",
        preview=f"{_COH}/coherence.review.png",
    params:
        label_q_max=_COH_CFG.get("fdr_panel_q_max", 0.05),
        label_quantile=_COH_CFG.get("fdr_panel_label_quantile", 0.05),
        label_max=_COH_CFG.get("fdr_panel_label_max", 5),
    log:
        "logs/3a_coherence/plot_{dataset}_{source}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [coherence] Plotting coherence figure for {wildcards.dataset} × {wildcards.source}..."
    shell:
        """
        python workflow/scripts/coherence/plot_coherence.py \
            --input {input.metrics} \
            --label-q-max {params.label_q_max} \
            --label-quantile {params.label_quantile} \
            --label-max {params.label_max} \
            --output {output.figure} &> {log}
        """


rule combine_coherence_metrics:
    input:
        # One per-source metrics table per registered source (aggregates the
        # {source} fan-out back into a single cross-source target).
        metrics=lambda wc: expand(
            f"results/3a_coherence/{wc.dataset}/{{source}}/coherence_metrics.parquet",
            source=_COH_SOURCES,
        ),
    output:
        combined="results/3a_coherence/{dataset}/coherence_metrics_combined.parquet",
    log:
        "logs/3a_coherence/combine_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coherence] Combining per-source metrics for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coherence/combine_metrics.py \
            --metrics {input.metrics} \
            --output {output.combined} &> {log}
        """


rule deduplicate_coherence_terms:
    input:
        combined="results/3a_coherence/{dataset}/coherence_metrics_combined.parquet",
        obo=lambda wc: (
            f"resources/external/pombase/{DATASETS['reference']['pombase_version']}"
            f"/ontologies_and_associations/go-basic.obo"
        ),
    output:
        all_terms="results/3a_coherence/{dataset}/coherence_terms_deduplicated.tsv",
        representatives="results/3a_coherence/{dataset}/coherence_terms_representatives.tsv",
    params:
        jaccard_threshold=_COH_CFG.get("dedup_jaccard_threshold", 0.5),
        lineage_flag="--merge-dag-lineage" if _COH_CFG.get("dedup_merge_dag_lineage", True) else "--no-merge-dag-lineage",
        scope=_COH_CFG.get("dedup_scope", "pooled"),
        force=lambda wc: " ".join(_COH_CFG.get("dedup_force_representatives", []) or []),
    log:
        "logs/3a_coherence/dedup_{dataset}.log",
    conda:
        # biopython.yml carries goatools (GO DAG depth + is_a/part_of lineage), pandas
        # and scipy (the redundancy graph runs through scipy.sparse.csgraph).
        "../envs/biopython.yml"
    message:
        "*** [coherence] De-duplicating coherence terms for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coherence/deduplicate_terms.py \
            --combined {input.combined} \
            --obo {input.obo} \
            --jaccard-threshold {params.jaccard_threshold} \
            {params.lineage_flag} \
            --scope {params.scope} \
            --force-representatives {params.force} \
            --output-all {output.all_terms} \
            --output-representatives {output.representatives} &> {log}
        """


rule compute_coherence_attribution:
    input:
        metrics=f"{_COH}/coherence_metrics.parquet",
        annotation=f"{_COH}/group_annotation_long.tsv",
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        paralogs="resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv",
    output:
        table=f"{_COH}/incoherence_attribution.tsv",
        points=f"{_COH}/incoherence_split_points.parquet",
    params:
        z_threshold=_COH_CFG.get("attribution_z_threshold", 0.0),
        shared_frac=_COH_CFG.get("attribution_shared_frac_threshold", 0.5),
        paralog_frac=_COH_CFG.get("attribution_paralog_frac_threshold", 0.5),
    log:
        "logs/3a_coherence/attribution_{dataset}_{source}.log",
    conda:
        # stats env: attribution needs sklearn (GMM) + pandas; goatools NOT needed here.
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coherence] Attributing incoherence for {wildcards.dataset} × {wildcards.source}..."
    shell:
        """
        python workflow/scripts/coherence/compute_incoherence_attribution.py \
            --metrics {input.metrics} \
            --annotation {input.annotation} \
            --fitting-results {input.fitting_results} \
            --paralogs {input.paralogs} \
            --z-threshold {params.z_threshold} \
            --shared-frac-threshold {params.shared_frac} \
            --paralog-frac-threshold {params.paralog_frac} \
            --output-table {output.table} \
            --output-points {output.points} &> {log}
        """


rule plot_coherence_attribution:
    input:
        table=lambda wc: f"results/3a_coherence/{wc.dataset}/{wc.source}/incoherence_attribution.tsv",
        points=lambda wc: f"results/3a_coherence/{wc.dataset}/{wc.source}/incoherence_split_points.parquet",
    output:
        figure=f"{_COH}/incoherence_attribution.pdf",
        preview=f"{_COH}/incoherence_attribution.review.png",
    params:
        top_n_plot=_COH_CFG.get("attribution_top_n_plot", 16),
    log:
        "logs/3a_coherence/attribution_plot_{dataset}_{source}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [coherence] Plotting incoherence attribution for {wildcards.dataset} × {wildcards.source}..."
    shell:
        """
        python workflow/scripts/coherence/plot_incoherence_attribution.py \
            --table {input.table} \
            --points {input.points} \
            --top-n-plot {params.top_n_plot} \
            --output {output.figure} &> {log}
        """


rule plot_coherence_group_scatter:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        annotation=f"{_COH}/group_annotation_long.tsv",
        metrics=f"{_COH}/coherence_metrics.parquet",
    output:
        figure=f"{_COH}/group_scatter.pdf",
        preview=f"{_COH}/group_scatter.review.png",
    params:
        # JSON-encode the namelist so it survives shell interpolation as a single
        # parseable literal (Snakemake would otherwise space-join a raw list, and
        # names contain spaces). JSON is double-quoted, so the '...' shell wrapper
        # below stays intact; plot_group_scatter.py parses it with ast.literal_eval.
        groups=lambda wc: json.dumps(_COH_CFG.get("scatter_groups", {}).get(wc.source, [])),
    log:
        "logs/3a_coherence/scatter_{dataset}_{source}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [coherence] Plotting scatter for {wildcards.dataset} × {wildcards.source}..."
    shell:
        """
        python workflow/scripts/coherence/plot_group_scatter.py \
            --fitting-results {input.fitting_results} \
            --annotation {input.annotation} \
            --metrics {input.metrics} \
            --source {wildcards.source} \
            --groups '{params.groups}' \
            --output-figure {output.figure} &> {log}
        """
