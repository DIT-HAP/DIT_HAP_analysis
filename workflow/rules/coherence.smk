# =============================================================================
# coherence.smk — Gene-group coherence in DR-DL fitness space
# =============================================================================
#
# Per (dataset × source): for every group (GO complex / CC / BP term) whose
# DR>dr_threshold members number min_group_size..max_group_size AND whose total
# annotated size is <= max_term_genes, measure how tightly they cluster in
# normalized (DR, DL/10) space vs a seeded permutation null (median pairwise
# distance z-score). Weiszfeld geometric median + seeded MPD permutation test
# live in workflow/src/coherence/metrics.py.
#
# Rules fanning out by {source}, then re-aggregated:
#   prepare_coherence_annotation -> long-form group annotation (per source)
#   compute_coherence             -> coherence_metrics.tsv + coherence_analysis.pdf
#                                    (p_value = one-sided add-one permutation p;
#                                    p_fdr = per-source Benjamini-Hochberg FDR)
#   combine_coherence_metrics     -> coherence_metrics_combined.tsv (all sources)
#   deduplicate_coherence_terms   -> coherence_terms_deduplicated.tsv (+ _representatives.tsv):
#                                    collapse redundant terms by member overlap + GO DAG
#                                    (display layer; full-set p_fdr untouched, all terms kept)
#   attribute_coherence_incoherence -> incoherence_attribution.tsv (+ .pdf): diagnose WHY a
#                                    complex is dispersed (major/minor GMM split, shared
#                                    subunits, paralog buffering). attribution_sources subset.
#   plot_coherence_group_scatter  -> group_scatter.pdf (named groups from config)
#
# Sources are registered in config.coherence.sources (currently: go_macrocomplex,
# go_cc, go_bp). Each source has a loader in workflow/src/coherence/sources.py.
# scatter_groups drives plot_group_scatter (per-source namelist, entries match
# group_name OR group_id).
#
# DATA-PATH NOTE: fitting_results comes from upstream DIT_HAP_snakemake release/
# dirs via DATASETS['datasets'][dataset]['release_dir']; pombase from
# resources/external/pombase/{DATASETS['reference']['pombase_version']}. Features
# (protein/RNA/structural) optionally drive additional diagnostic panels in the
# coherence figure when config.coherence.features_panels is true.

import json

_COH = "results/coherence/{dataset}/{source}"
_COH_CFG = config.get("coherence", {})
_COH_SOURCES = _COH_CFG.get("sources", ["go_macrocomplex", "go_cc", "go_bp"])
# Incoherence attribution runs on a (usually smaller) subset of sources: major/minor
# split + shared-subunit are only meaningful for physical complexes, not GO_BP processes.
_COH_ATTR_SOURCES = _COH_CFG.get("attribution_sources", ["go_macrocomplex"])

wildcard_constraints:
    source="|".join(_COH_SOURCES),


rule prepare_coherence_annotation:
    input:
        pombase_dir=lambda wc: f"resources/external/pombase/{DATASETS['reference']['pombase_version']}",
    output:
        long_table=f"{_COH}/group_annotation_long.tsv",
    log:
        "logs/coherence/prepare_{dataset}_{source}.log",
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
    output:
        f"{_COH}/coherence.parquet",
    params:
        min_size=_COH_CFG.get("min_group_size", 3),
        max_size=_COH_CFG.get("max_group_size", 300),
        max_term_genes=_COH_CFG.get("max_term_genes", 500),
        dr_threshold=_COH_CFG.get("dr_threshold", 0.3),
        n_permutations=_COH_CFG.get("n_permutations", 1000),
        random_state=_COH_CFG.get("random_state", 42),
    log:
        "logs/coherence/compute_{dataset}_{source}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coherence] Computing coherence metrics for {wildcards.dataset} × {wildcards.source}..."
    shell:
        """
        python workflow/scripts/coherence/compute_coherence.py \
            --fitting-results {input.fitting_results} \
            --annotation {input.annotation} \
            --source {wildcards.source} \
            --min-size {params.min_size} \
            --max-size {params.max_size} \
            --max-term-genes {params.max_term_genes} \
            --dr-threshold {params.dr_threshold} \
            --n-permutations {params.n_permutations} \
            --random-state {params.random_state} \
            --output {output} &> {log}
        """


rule combine_coherence_metrics:
    input:
        # One per-source metrics table per registered source (aggregates the
        # {source} fan-out back into a single cross-source target).
        metrics=lambda wc: expand(
            f"results/coherence/{wc.dataset}/{{source}}/coherence_metrics.tsv",
            source=_COH_SOURCES,
        ),
    output:
        combined="results/coherence/{dataset}/coherence_metrics_combined.tsv",
    log:
        "logs/coherence/combine_{dataset}.log",
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
        combined="results/coherence/{dataset}/coherence_metrics_combined.tsv",
        obo=lambda wc: (
            f"resources/external/pombase/{DATASETS['reference']['pombase_version']}"
            f"/ontologies_and_associations/go-basic.obo"
        ),
    output:
        all_terms="results/coherence/{dataset}/coherence_terms_deduplicated.tsv",
        representatives="results/coherence/{dataset}/coherence_terms_representatives.tsv",
    params:
        overlap_threshold=_COH_CFG.get("dedup_overlap_threshold", 0.5),
        lineage_flag="--merge-dag-lineage" if _COH_CFG.get("dedup_merge_dag_lineage", True) else "",
        scope=_COH_CFG.get("dedup_scope", "pooled"),
        force=lambda wc: " ".join(_COH_CFG.get("dedup_force_representatives", []) or []),
    log:
        "logs/coherence/dedup_{dataset}.log",
    conda:
        # biopython.yml carries goatools (GO DAG depth + is_a/part_of lineage) + pandas.
        "../envs/biopython.yml"
    message:
        "*** [coherence] De-duplicating coherence terms for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coherence/deduplicate_terms.py \
            --combined {input.combined} \
            --obo {input.obo} \
            --overlap-threshold {params.overlap_threshold} \
            {params.lineage_flag} \
            --scope {params.scope} \
            --force-representatives {params.force} \
            --output-all {output.all_terms} \
            --output-representatives {output.representatives} &> {log}
        """


rule attribute_coherence_incoherence:
    input:
        metrics=f"{_COH}/coherence_metrics.tsv",
        annotation=f"{_COH}/group_annotation_long.tsv",
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        paralogs="resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv",
    output:
        table=f"{_COH}/incoherence_attribution.tsv",
        figure=f"{_COH}/incoherence_attribution.pdf",
    params:
        z_threshold=_COH_CFG.get("attribution_z_threshold", 0.0),
        top_n_plot=_COH_CFG.get("attribution_top_n_plot", 16),
        shared_frac=_COH_CFG.get("attribution_shared_frac_threshold", 0.5),
        paralog_frac=_COH_CFG.get("attribution_paralog_frac_threshold", 0.5),
    log:
        "logs/coherence/attribution_{dataset}_{source}.log",
    conda:
        # stats env: attribution needs sklearn (GMM) + matplotlib; goatools NOT needed here.
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coherence] Attributing incoherence for {wildcards.dataset} × {wildcards.source}..."
    shell:
        """
        python workflow/scripts/coherence/attribute_incoherence.py \
            --metrics {input.metrics} \
            --annotation {input.annotation} \
            --fitting-results {input.fitting_results} \
            --paralogs {input.paralogs} \
            --z-threshold {params.z_threshold} \
            --top-n-plot {params.top_n_plot} \
            --shared-frac-threshold {params.shared_frac} \
            --paralog-frac-threshold {params.paralog_frac} \
            --output-table {output.table} \
            --output-figure {output.figure} &> {log}
        """


rule plot_coherence_group_scatter:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        annotation=f"{_COH}/group_annotation_long.tsv",
        metrics=f"{_COH}/coherence_metrics.tsv",
    output:
        figure=f"{_COH}/group_scatter.pdf",
    params:
        # JSON-encode the namelist so it survives shell interpolation as a single
        # parseable literal (Snakemake would otherwise space-join a raw list, and
        # names contain spaces). JSON is double-quoted, so the '...' shell wrapper
        # below stays intact; plot_group_scatter.py parses it with ast.literal_eval.
        groups=lambda wc: json.dumps(_COH_CFG.get("scatter_groups", {}).get(wc.source, [])),
    log:
        "logs/coherence/scatter_{dataset}_{source}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
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
