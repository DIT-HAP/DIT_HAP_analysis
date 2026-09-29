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
#   plot_coherence_dedup           -> coherence_dedup.pdf (the representative set alone)
#   plot_coherence_by_source       -> coherence_by_source.pdf (every source + the
#                                     representative set on one figure, coloured by source)
#   compute_coherence_attribution  -> incoherence_attribution.tsv + incoherence_split_points.parquet:
#                                     diagnose WHY a complex is dispersed (major/minor GMM
#                                     split, shared subunits, paralog buffering).
#                                     attribution_sources subset.
#   plot_coherence_attribution     -> incoherence_attribution.pdf (renders the persisted split points)
#   compute_coherence_attribution_dedup -> dedup_incoherence_attribution.tsv
#                                     (+ dedup_incoherence_split_points.parquet): the same
#                                     diagnosis over the pooled representative set — every
#                                     source's annotation concatenated, keyed (source, group_id)
#                                     because group_id is NOT unique across sources.
#   plot_coherence_attribution_dedup -> dedup_incoherence_attribution.pdf
#   plot_coherence_group_scatter   -> group_scatter.pdf (named groups from config)
#   plot_coherence_interactive_scatter -> interactive_scatter.html (Altair; term picker
#                                     + per-gene tooltips over the genome-wide cloud)
#   plot_coherence_redundancy_network  -> redundancy_network.html (Altair; the merge
#                                     backbone of each redundancy cluster, edge = Jaccard)
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
# go_cc, go_bp, kegg_brite, kegg_pathway). Each source has a loader in
# workflow/src/coherence/sources.py; the kegg_* ones read the kegg_parser derived
# tables under config.coherence.kegg_dir instead of PomBase. scatter_groups drives
# plot_group_scatter (per-source namelist, entries match group_name OR group_id).
#
# DATA-PATH NOTE: fitting_results comes from upstream DIT_HAP_snakemake release/
# dirs via DATASETS['datasets'][dataset]['release_dir']; pombase from
# resources/external/pombase/{DATASETS['reference']['pombase_version']}.

import json

_COH = "results/3a_coherence/{dataset}/{source}"
_COH_CFG = config.get("coherence", {})
_COH_SOURCES = _COH_CFG.get("sources", ["go_macrocomplex", "go_cc", "go_bp"])
# Incoherence attribution runs on a configurable subset of sources (currently all
# of them): the major/minor + shared-subunit signals are sharpest for physical
# complexes, but the label ladder still separates a genuinely two-population term
# from a broadly-annotated one for the process / pathway sources.
_COH_ATTR_SOURCES = _COH_CFG.get("attribution_sources", ["go_macrocomplex"])
_COH_FEATURES_PANELS = _COH_CFG.get("features_panels", True)
_COH_FEATURES = (
    f"results/1b_features/{DATASETS['reference']['pombase_version']}"
    f"/pombe_coding_gene_protein_features.tsv"
)
# The kegg_* adapters read the kegg_parser derived tables instead of PomBase; the
# GO ones ignore --kegg-dir entirely, so the input stays off for them (same
# optional-input idiom as compute_coherence's features).
_COH_KEGG_DIR = _COH_CFG.get("kegg_dir", "resources/external/kegg/data/derived")
_COH_KEGG_SOURCES = {src for src in _COH_SOURCES if src.startswith("kegg_")}

wildcard_constraints:
    source="|".join(_COH_SOURCES),


rule prepare_coherence_annotation:
    input:
        pombase_dir=lambda wc: f"resources/external/pombase/{DATASETS['reference']['pombase_version']}",
        kegg_dir=lambda wc: _COH_KEGG_DIR if wc.source in _COH_KEGG_SOURCES else [],
    output:
        long_table=f"{_COH}/group_annotation_long.tsv",
    params:
        kegg_flag=lambda wc, input: f"--kegg-dir {input.kegg_dir}" if input.kegg_dir else "",
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
            {params.kegg_flag} \
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
        label_q_max=_COH_CFG.get("coherent_q_max", 0.05),
        label_quantile=_COH_CFG.get("fdr_panel_label_quantile", 0.05),
        label_max=_COH_CFG.get("fdr_panel_label_max", 5),
        label_coherent_z=_COH_CFG.get("coherent_z_threshold", -2.0),
        label_incoherent_z=_COH_CFG.get("incoherent_z_threshold", 1.0),
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
            --label-coherent-z {params.label_coherent_z} \
            --label-incoherent-z {params.label_incoherent_z} \
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
        # Fallbacks must match the script's own defaults, or a config that dropped
        # the key would silently run a different algorithm than a bare CLI call.
        linkage=_COH_CFG.get("dedup_linkage", "complete"),
        lineage_flag="--merge-dag-lineage" if _COH_CFG.get("dedup_merge_dag_lineage", False) else "--no-merge-dag-lineage",
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
            --linkage {params.linkage} \
            {params.lineage_flag} \
            --scope {params.scope} \
            --force-representatives {params.force} \
            --output-all {output.all_terms} \
            --output-representatives {output.representatives} &> {log}
        """


# --- The de-duplicated set as a first-class result --------------------------
# Everything above is per source; the representative set is the cross-source end
# product, so it gets its own coherence figure and its own attribution run. The
# TSV's pooled q_value is carried through untouched from the combined table — see
# combine_metrics.py's FDR note for why re-correcting over the representatives
# would be anti-conservative rather than merely redundant.

rule plot_coherence_dedup:
    input:
        representatives="results/3a_coherence/{dataset}/coherence_terms_representatives.tsv",
    output:
        figure="results/3a_coherence/{dataset}/coherence_dedup.pdf",
        preview="results/3a_coherence/{dataset}/coherence_dedup.review.png",
    params:
        label_q_max=_COH_CFG.get("coherent_q_max", 0.05),
        label_quantile=_COH_CFG.get("fdr_panel_label_quantile", 0.05),
        label_max=_COH_CFG.get("fdr_panel_label_max", 5),
        label_coherent_z=_COH_CFG.get("coherent_z_threshold", -2.0),
        label_incoherent_z=_COH_CFG.get("incoherent_z_threshold", 1.0),
    log:
        "logs/3a_coherence/plot_dedup_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [coherence] Plotting the de-duplicated representative set for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coherence/plot_coherence.py \
            --input {input.representatives} \
            --label-q-max {params.label_q_max} \
            --label-quantile {params.label_quantile} \
            --label-max {params.label_max} \
            --label-coherent-z {params.label_coherent_z} \
            --label-incoherent-z {params.label_incoherent_z} \
            --output {output.figure} &> {log}
        """


rule plot_coherence_by_source:
    input:
        combined="results/3a_coherence/{dataset}/coherence_metrics_combined.parquet",
        representatives="results/3a_coherence/{dataset}/coherence_terms_representatives.tsv",
    output:
        figure="results/3a_coherence/{dataset}/coherence_by_source.pdf",
        preview="results/3a_coherence/{dataset}/coherence_by_source.review.png",
    params:
        label_q_max=_COH_CFG.get("coherent_q_max", 0.05),
        label_quantile=_COH_CFG.get("fdr_panel_label_quantile", 0.05),
        label_max=_COH_CFG.get("fdr_panel_label_max", 5),
        label_coherent_z=_COH_CFG.get("coherent_z_threshold", -2.0),
        label_incoherent_z=_COH_CFG.get("incoherent_z_threshold", 1.0),
    log:
        "logs/3a_coherence/plot_by_source_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [coherence] Plotting the cross-source comparison for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coherence/plot_coherence.py \
            --input {input.combined} \
            --color-by source \
            --dedup-series {input.representatives} \
            --label-q-max {params.label_q_max} \
            --label-quantile {params.label_quantile} \
            --label-max {params.label_max} \
            --label-coherent-z {params.label_coherent_z} \
            --label-incoherent-z {params.label_incoherent_z} \
            --output {output.figure} &> {log}
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
        z_threshold=_COH_CFG.get("incoherent_z_threshold", 1.0),
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
        # The genome-wide cloud behind every panel, and what fixes the panels' shared
        # axis ranges — the same upstream table the attribution itself was computed from.
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
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
            --fitting-results {input.fitting_results} \
            --top-n-plot {params.top_n_plot} \
            --output {output.figure} &> {log}
        """


# The de-duplicated set's attribution. Same diagnosis as the per-source rule above,
# run once over the pooled representatives instead of five times over the full
# per-source tables: the metrics TSV carries the pooled q, and every source's
# annotation is concatenated so a group_id shared by two sources (173 of them)
# still resolves to its own source's members. `--top-n-plot` is still a global
# z-ranking, so the panels can all come from one source — each panel's title now
# names its source, which is what makes that visible rather than misleading.
rule compute_coherence_attribution_dedup:
    input:
        metrics="results/3a_coherence/{dataset}/coherence_terms_representatives.tsv",
        annotations=lambda wc: expand(
            f"results/3a_coherence/{wc.dataset}/{{source}}/group_annotation_long.tsv",
            source=_COH_SOURCES,
        ),
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        paralogs="resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv",
    output:
        table="results/3a_coherence/{dataset}/dedup_incoherence_attribution.tsv",
        points="results/3a_coherence/{dataset}/dedup_incoherence_split_points.parquet",
    params:
        z_threshold=_COH_CFG.get("incoherent_z_threshold", 1.0),
        shared_frac=_COH_CFG.get("attribution_shared_frac_threshold", 0.5),
        paralog_frac=_COH_CFG.get("attribution_paralog_frac_threshold", 0.5),
    log:
        "logs/3a_coherence/attribution_dedup_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coherence] Attributing incoherence for the de-duplicated set of {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coherence/compute_incoherence_attribution.py \
            --metrics {input.metrics} \
            --annotation {input.annotations} \
            --fitting-results {input.fitting_results} \
            --paralogs {input.paralogs} \
            --z-threshold {params.z_threshold} \
            --shared-frac-threshold {params.shared_frac} \
            --paralog-frac-threshold {params.paralog_frac} \
            --output-table {output.table} \
            --output-points {output.points} &> {log}
        """


rule plot_coherence_attribution_dedup:
    input:
        table="results/3a_coherence/{dataset}/dedup_incoherence_attribution.tsv",
        points="results/3a_coherence/{dataset}/dedup_incoherence_split_points.parquet",
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
    output:
        figure="results/3a_coherence/{dataset}/dedup_incoherence_attribution.pdf",
        preview="results/3a_coherence/{dataset}/dedup_incoherence_attribution.review.png",
    params:
        top_n_plot=_COH_CFG.get("attribution_top_n_plot", 16),
    log:
        "logs/3a_coherence/attribution_plot_dedup_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [coherence] Plotting the de-duplicated incoherence attribution for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coherence/plot_incoherence_attribution.py \
            --table {input.table} \
            --points {input.points} \
            --fitting-results {input.fitting_results} \
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


# --- Interactive HTML (Altair) ----------------------------------------------
# These two emit .html rather than the .pdf + .review.png pair, and neither goes
# through figures.py's save_dual(). They are exploration tools, not figure
# panels: one page per source for the scatter, one page per dataset for the
# network, each with a dropdown and hover tooltips.
#
# The scatter page embeds every (term, member) row in the source, which is 1.6k
# rows for go_macrocomplex but 58k for go_bp — a browser-loadable but heavy page.
# It is therefore NOT in `rule all`: build it for the source you are looking at
# (`snakemake --use-conda results/3a_coherence/<dataset>/go_cc/interactive_scatter.html`).
# The network page is a few tens of kilobytes, so it is a default target.

rule plot_coherence_interactive_scatter:
    input:
        metrics=f"{_COH}/coherence_metrics.parquet",
        annotation=f"{_COH}/group_annotation_long.tsv",
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
    output:
        page=f"{_COH}/interactive_scatter.html",
    log:
        "logs/3a_coherence/interactive_scatter_{dataset}_{source}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [coherence] Writing the interactive scatter for {wildcards.dataset} × {wildcards.source}..."
    shell:
        """
        python workflow/scripts/coherence/plot_interactive_scatter.py \
            --metrics {input.metrics} \
            --annotation {input.annotation} \
            --fitting-results {input.fitting_results} \
            --source {wildcards.source} \
            --output {output.page} &> {log}
        """


rule plot_coherence_redundancy_network:
    input:
        combined=f"results/3a_coherence/{{dataset}}/coherence_metrics_combined.parquet",
        deduplicated=f"results/3a_coherence/{{dataset}}/coherence_terms_deduplicated.tsv",
    output:
        page="results/3a_coherence/{dataset}/redundancy_network.html",
        overview="results/3a_coherence/{dataset}/redundancy_overview.html",
    log:
        "logs/3a_coherence/redundancy_network_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [coherence] Writing the redundancy network for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coherence/plot_redundancy_network.py \
            --combined {input.combined} \
            --deduplicated {input.deduplicated} \
            --output {output.page} \
            --output-overview {output.overview} &> {log}
        """
