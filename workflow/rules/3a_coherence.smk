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
# NAMING: the per-source form is the DEFAULT, so a rule that fans out by {source}
# carries no scope suffix (`compute_coherence`, `plot_coherence`, ...). The two
# pooled views are marked on the rule name: `_for_combined` (all sources in one
# table, one q family) and `_for_dedup` (the de-duplicated representative set).
#
# Per-source rules (fan out over {source}):
#   prepare_coherence_annotation   -> group_annotation_long.tsv
#   compute_coherence              -> coherence_metrics.parquet
#                                     ({method}_p = one-sided add-one permutation p;
#                                     q_value = per-source Benjamini-Hochberg FDR
#                                     over the primary method's p;
#                                     plus, when features_panels is on, the abundance /
#                                     conservation uniformity CVs and the
#                                     paralog_fraction of the scored members —
#                                     all computed here so the figure only renders)
#   export_coherence_cohorts       -> {source}/coherence_cohorts.xlsx (thresholded
#                                     straight off the table above — this source's
#                                     own q, no de-duplication columns)
#   plot_coherence                 -> coherence.pdf (reads the metrics table alone)
#   plot_coherence_group_scatter   -> group_scatter.pdf (named groups from config)
#   plot_coherence_interactive_scatter -> interactive_scatter.html (Altair; term picker
#                                     + per-gene tooltips over the genome-wide cloud)
#
# Per-view rules (fan out over {view}: the five sources + combined + dedup):
#   compute_view_fractions         -> {view}/view_fractions.tsv + view_gene_breadth.tsv:
#                                     the per-view paralog / moonlighting fractions and
#                                     the gene-breadth counts behind the latter.
#   plot_fraction_distributions    -> {view}/fraction_distributions.pdf (paralog fraction,
#                                     groups per gene with the moonlighting cut marked,
#                                     moonlighting fraction)
#
# Cross-source rules (the {source} fan-out re-aggregated):
#   combine_coherence_metrics      -> combined/coherence_metrics.parquet (all sources;
#                                     FDR re-derived over the union — a q is only
#                                     defined w.r.t. a family, and this table's family
#                                     is its own row set, so cross-source thresholds
#                                     and the dedup representative ranking share one
#                                     q scale. Per-source tables keep their own q.)
#   deduplicate_coherence_terms    -> dedup/coherence_terms_deduplicated.tsv (+ _representatives.tsv
#                                     + coherence_group_members_long.tsv, the (group, gene)
#                                     long table with the moonlighting columns):
#                                     collapse redundant terms by member overlap + GO DAG
#                                     (display layer; the combined q_value is carried
#                                     through untouched, all terms kept)
#   plot_coherence_for_combined    -> combined/coherence.pdf (every source + the
#                                     representative set on one figure, coloured by source)
#   export_coherence_cohorts_for_combined -> combined/coherence_cohorts.xlsx (the
#                                     full term list on the pooled q, + in_dedup_set,
#                                     moonlighting_fraction)
#   plot_coherence_for_dedup       -> dedup/coherence.pdf (the representative set alone)
#   export_coherence_cohorts_for_dedup -> dedup/coherence_cohorts.xlsx
#   plot_coherence_redundancy_network_for_dedup -> dedup/redundancy_network.html (Altair; the merge
#                                     backbone of each redundancy cluster, edge = Jaccard)
#   plot_coherence_interactive_scatter_for_pooled -> combined|dedup/interactive_scatter.html (the
#                                     per-source Altair explorer re-run over the pooled rows)
#
# VIEW FOLDERS: results/3a_coherence/{dataset}/ holds one folder per view, and the
# folder name — never the file name — says which view a file belongs to:
#   {source}/  (go_macrocomplex, go_bp, ...)  one annotated source's rows
#   combined/                                 all sources pooled, one q family
#   dedup/                                    the de-duplicated representative set
# Every view exposes the same file names where the artifact exists for it, so
# `combined/coherence.pdf` and `go_cc/coherence.pdf` are the same figure over
# different row sets. Nothing stays at the dataset root: every artifact belongs to
# a view, and the paths inside a rule are what keep the views apart.
#
# Data-format rule: per-stage intermediates are Parquet (exact dtypes survive
# round-trip); only the final human-facing artifacts (dedup + view-fraction tables,
# and the cohort workbooks) are TSV / Excel.
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
# The cohort definition, read once. The figures that label with it and the three
# cohort workbooks all consume it, so it lives here rather than as a
# `.get(key, default)` repeated in every params block. Changing the config key
# still rebuilds all of them: Snakemake's default rerun-triggers include `params`.
_COH_COHERENT_Q_MAX = _COH_CFG.get("coherent_q_max", 0.05)
_COH_COHERENT_Z = _COH_CFG.get("coherent_z_threshold", -2.0)
_COH_INCOHERENT_Z = _COH_CFG.get("incoherent_z_threshold", 1.0)
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


# The cohort workbook for one source, straight off the metrics table just computed
# above — no combine step in between. Its q is that source's own BH family, so the
# workbook is NOT a slice of combined/'s (58 of the 322 coherent groups flip
# cohort between the two q families, 2026-09-29), and being one source's raw term
# list it carries no de-duplication column either.
rule export_coherence_cohorts:
    input:
        metrics=f"{_COH}/coherence_metrics.parquet",
    output:
        xlsx=f"{_COH}/coherence_cohorts.xlsx",
    params:
        q_max=_COH_COHERENT_Q_MAX,
        coherent_z=_COH_COHERENT_Z,
        incoherent_z=_COH_INCOHERENT_Z,
    log:
        "logs/3a_coherence/export_cohorts_{dataset}_{source}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coherence] Exporting cohort tables for {wildcards.dataset} × {wildcards.source}..."
    shell:
        """
        python workflow/scripts/coherence/export_cohorts.py \
            --metrics {input.metrics} \
            --source {wildcards.source} \
            --q-max {params.q_max} \
            --coherent-z {params.coherent_z} \
            --incoherent-z {params.incoherent_z} \
            --output {output.xlsx} &> {log}
        """


rule plot_coherence:
    input:
        metrics=f"{_COH}/coherence_metrics.parquet",
    output:
        figure=f"{_COH}/coherence.pdf",
        preview=f"{_COH}/coherence.review.png",
    params:
        label_q_max=_COH_COHERENT_Q_MAX,
        label_quantile=_COH_CFG.get("fdr_panel_label_quantile", 0.05),
        label_max=_COH_CFG.get("fdr_panel_label_max", 5),
        label_coherent_z=_COH_COHERENT_Z,
        label_incoherent_z=_COH_INCOHERENT_Z,
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
        combined="results/3a_coherence/{dataset}/combined/coherence_metrics.parquet",
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
        combined="results/3a_coherence/{dataset}/combined/coherence_metrics.parquet",
        obo=lambda wc: (
            f"resources/external/pombase/{DATASETS['reference']['pombase_version']}"
            f"/ontologies_and_associations/go-basic.obo"
        ),
    output:
        all_terms="results/3a_coherence/{dataset}/dedup/coherence_terms_deduplicated.tsv",
        representatives="results/3a_coherence/{dataset}/dedup/coherence_terms_representatives.tsv",
        group_members="results/3a_coherence/{dataset}/dedup/coherence_group_members_long.tsv",
    params:
        jaccard_threshold=_COH_CFG.get("dedup_jaccard_threshold", 0.5),
        # Fallbacks must match the script's own defaults, or a config that dropped
        # the key would silently run a different algorithm than a bare CLI call.
        linkage=_COH_CFG.get("dedup_linkage", "complete"),
        scope=_COH_CFG.get("dedup_scope", "pooled"),
        force=lambda wc: " ".join(_COH_CFG.get("dedup_force_representatives", []) or []),
    log:
        "logs/3a_coherence/dedup_{dataset}.log",
    conda:
        # biopython.yml carries goatools (GO DAG depth), pandas
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
            --scope {params.scope} \
            --force-representatives {params.force} \
            --output-all {output.all_terms} \
            --output-representatives {output.representatives} \
            --output-group-members {output.group_members} &> {log}
        """


# --- The de-duplicated set as a first-class result --------------------------
# The representative set is the cross-source end product, so it gets its own
# coherence figure and cohort workbook. The
# TSV's pooled q_value is carried through untouched from the combined table — see
# combine_metrics.py's FDR note for why re-correcting over the representatives
# would be anti-conservative rather than merely redundant.

rule plot_coherence_for_dedup:
    input:
        representatives="results/3a_coherence/{dataset}/dedup/coherence_terms_representatives.tsv",
    output:
        figure="results/3a_coherence/{dataset}/dedup/coherence.pdf",
        preview="results/3a_coherence/{dataset}/dedup/coherence.review.png",
    params:
        label_q_max=_COH_COHERENT_Q_MAX,
        label_quantile=_COH_CFG.get("fdr_panel_label_quantile", 0.05),
        label_max=_COH_CFG.get("fdr_panel_label_max", 5),
        label_coherent_z=_COH_COHERENT_Z,
        label_incoherent_z=_COH_INCOHERENT_Z,
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


rule plot_coherence_for_combined:
    input:
        combined="results/3a_coherence/{dataset}/combined/coherence_metrics.parquet",
        representatives="results/3a_coherence/{dataset}/dedup/coherence_terms_representatives.tsv",
    output:
        figure="results/3a_coherence/{dataset}/combined/coherence.pdf",
        preview="results/3a_coherence/{dataset}/combined/coherence.review.png",
    params:
        label_q_max=_COH_COHERENT_Q_MAX,
        label_quantile=_COH_CFG.get("fdr_panel_label_quantile", 0.05),
        label_max=_COH_CFG.get("fdr_panel_label_max", 5),
        label_coherent_z=_COH_COHERENT_Z,
        label_incoherent_z=_COH_INCOHERENT_Z,
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


rule export_coherence_cohorts_for_combined:
    input:
        metrics="results/3a_coherence/{dataset}/combined/coherence_metrics.parquet",
        # The all-terms dedup table, not the representatives-only one: it carries
        # is_representative (the in_dedup_set flag) AND moonlighting_fraction, which
        # the workbook then shows for every term. Same upstream (the dedup rule),
        # so this costs no new dependency.
        dedup_terms="results/3a_coherence/{dataset}/dedup/coherence_terms_deduplicated.tsv",
    output:
        xlsx="results/3a_coherence/{dataset}/combined/coherence_cohorts.xlsx",
    params:
        q_max=_COH_COHERENT_Q_MAX,
        coherent_z=_COH_COHERENT_Z,
        incoherent_z=_COH_INCOHERENT_Z,
    log:
        "logs/3a_coherence/export_cohorts_{dataset}.log",
    conda:
        # stats env: openpyxl for the workbook is in this recipe.
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coherence] Exporting threshold-filtered cohort tables for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coherence/export_cohorts.py \
            --metrics {input.metrics} \
            --dedup-terms {input.dedup_terms} \
            --q-max {params.q_max} \
            --coherent-z {params.coherent_z} \
            --incoherent-z {params.incoherent_z} \
            --output {output.xlsx} &> {log}
        """


# The same cut applied to the de-duplicated representative set (957 groups pooled
# over all sources) instead of the full 2,587. Reads the representatives TSV, whose
# q_value is the pooled one carried through untouched, so this workbook and the
# combined one call the same cohort for the same group. No --dedup-terms flag
# here: every row of that table already IS a representative, so the column would be
# a constant True (it still carries paralog_fraction and moonlighting_fraction, as
# columns of the table itself).
rule export_coherence_cohorts_for_dedup:
    input:
        metrics="results/3a_coherence/{dataset}/dedup/coherence_terms_representatives.tsv",
    output:
        xlsx="results/3a_coherence/{dataset}/dedup/coherence_cohorts.xlsx",
    params:
        q_max=_COH_COHERENT_Q_MAX,
        coherent_z=_COH_COHERENT_Z,
        incoherent_z=_COH_INCOHERENT_Z,
    log:
        "logs/3a_coherence/export_cohorts_dedup_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coherence] Exporting cohort tables for the de-duplicated set of {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coherence/export_cohorts.py \
            --metrics {input.metrics} \
            --q-max {params.q_max} \
            --coherent-z {params.coherent_z} \
            --incoherent-z {params.incoherent_z} \
            --output {output.xlsx} &> {log}
        """


# --- Per-view fraction distributions ----------------------------------------
# `paralog_fraction` (per term, from the feature matrix) and `moonlighting_fraction`
# (per term, from the gene-breadth mode) are descriptive per-view quantities, so
# their distributions are drawn for every view: the five sources, the pooled
# combined table, and the de-duplicated representative set. One rule pair fans out
# over `view` because the computation is identical for all of them — a view is a
# slice of the combined table plus the set of de-duplication clusters its terms
# belong to, and the breadth counts run over the (group, gene) long table.
_VIEWS = [*_COH_SOURCES, "combined", "dedup"]

rule compute_view_fractions:
    input:
        metrics="results/3a_coherence/{dataset}/combined/coherence_metrics.parquet",
        dedup_terms="results/3a_coherence/{dataset}/dedup/coherence_terms_deduplicated.tsv",
        group_members="results/3a_coherence/{dataset}/dedup/coherence_group_members_long.tsv",
    output:
        terms="results/3a_coherence/{dataset}/{view}/view_fractions.tsv",
        genes="results/3a_coherence/{dataset}/{view}/view_gene_breadth.tsv",
    wildcard_constraints:
        view="|".join(_VIEWS),
    log:
        "logs/3a_coherence/view_fractions_{dataset}_{view}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coherence] Computing paralog / moonlighting fractions for {wildcards.dataset} × {wildcards.view}..."
    shell:
        """
        python workflow/scripts/coherence/compute_view_fractions.py \
            --metrics {input.metrics} \
            --dedup-terms {input.dedup_terms} \
            --group-members {input.group_members} \
            --view {wildcards.view} \
            --output-terms {output.terms} \
            --output-genes {output.genes} &> {log}
        """


rule plot_fraction_distributions:
    input:
        terms="results/3a_coherence/{dataset}/{view}/view_fractions.tsv",
        genes="results/3a_coherence/{dataset}/{view}/view_gene_breadth.tsv",
    output:
        figure="results/3a_coherence/{dataset}/{view}/fraction_distributions.pdf",
        preview="results/3a_coherence/{dataset}/{view}/fraction_distributions.review.png",
    wildcard_constraints:
        view="|".join(_VIEWS),
    log:
        "logs/3a_coherence/fraction_distributions_{dataset}_{view}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [coherence] Plotting the fraction distributions for {wildcards.dataset} × {wildcards.view}..."
    shell:
        """
        python workflow/scripts/coherence/plot_fraction_distributions.py \
            --terms {input.terms} \
            --genes {input.genes} \
            --view {wildcards.view} \
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
# The scatter page is a two-panel explorer: the terms' centroid positions on the
# left (click a marker to pick one) and the picked term's genes over the
# genome-wide cloud on the right, both driven by one selection. Its payload was
# cut from 24 MB to 6.9 MB for go_bp by embedding the gene columns once as the
# cloud's dataset and looking them up from the member rows, so it stays out of
# `rule all` only because a 58k-row page is not worth pre-building for four
# sources you are not reading: build it for the source you are looking at
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


# The pooled views get the same page over their own metrics table. One rule for
# both, because the two differ only in the file they read: `--source` is simply
# omitted, which tells the script to take every row of the concatenated annotation
# and to prefix each term with its source (a name is only unique within a source,
# and 173 group_ids are shared by two of them). The scripts/ dir keeps the
# per-source split for rules whose shell command actually differs.
_POOLED_METRICS = {
    "combined": "coherence_metrics.parquet",
    "dedup": "coherence_terms_representatives.tsv",
}

rule plot_coherence_interactive_scatter_for_pooled:
    input:
        metrics=lambda wc: (
            f"results/3a_coherence/{wc.dataset}/{wc.view}/{_POOLED_METRICS[wc.view]}"
        ),
        # Every source's annotation, concatenated by the script: the pooled metrics
        # table carries groups from all of them, and the members have to come along.
        annotations=lambda wc: expand(
            f"results/3a_coherence/{wc.dataset}/{{source}}/group_annotation_long.tsv",
            source=_COH_SOURCES,
        ),
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
    output:
        page="results/3a_coherence/{dataset}/{view}/interactive_scatter.html",
    wildcard_constraints:
        view="combined|dedup",
    log:
        "logs/3a_coherence/interactive_scatter_{dataset}_{view}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [coherence] Writing the interactive scatter for {wildcards.dataset} × {wildcards.view}..."
    shell:
        """
        python workflow/scripts/coherence/plot_interactive_scatter.py \
            --metrics {input.metrics} \
            --annotation {input.annotations} \
            --fitting-results {input.fitting_results} \
            --output {output.page} &> {log}
        """


rule plot_coherence_redundancy_network_for_dedup:
    input:
        combined=f"results/3a_coherence/{{dataset}}/combined/coherence_metrics.parquet",
        deduplicated=f"results/3a_coherence/{{dataset}}/dedup/coherence_terms_deduplicated.tsv",
    output:
        page="results/3a_coherence/{dataset}/dedup/redundancy_network.html",
        overview="results/3a_coherence/{dataset}/dedup/redundancy_overview.html",
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
