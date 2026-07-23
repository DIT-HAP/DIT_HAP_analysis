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
# Three-rule chain fanning out by {source}:
#   prepare_coherence_annotation -> long-form group annotation (per source)
#   compute_coherence             -> coherence_metrics.tsv + coherence_analysis.pdf
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

_COH = "results/coherence/{dataset}/{source}"
_COH_CFG = config.get("coherence", {})
_COH_SOURCES = _COH_CFG.get("sources", ["go_macrocomplex", "go_cc", "go_bp"])

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
        features=lambda wc: f"results/features/{DATASETS['reference']['pombase_version']}/pombe_coding_gene_protein_features.tsv",
    output:
        metrics=f"{_COH}/coherence_metrics.tsv",
        figure=f"{_COH}/coherence_analysis.pdf",
    params:
        min_size=_COH_CFG.get("min_group_size", 3),
        max_size=_COH_CFG.get("max_group_size", 300),
        max_term_genes=_COH_CFG.get("max_term_genes", 500),
        dr_threshold=_COH_CFG.get("dr_threshold", 0.3),
        n_permutations=_COH_CFG.get("n_permutations", 1000),
        random_state=_COH_CFG.get("random_state", 42),
        features_flag=lambda wc, input: (
            f"--features {input.features}"
            if _COH_CFG.get("features_panels", True) else ""
        ),
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
            {params.features_flag} \
            --output-metrics {output.metrics} \
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
        groups=lambda wc: _COH_CFG.get("scatter_groups", {}).get(wc.source, []),
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
