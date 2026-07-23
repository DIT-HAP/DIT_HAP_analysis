# =============================================================================
# complex.smk — Macromolecular complex coherence analysis
# =============================================================================
#
# Batch B. Sources the DIT-HAP gene fitness table from the clustering
# finalize-variant system via final_clusters_path(dataset, selected_variant(dataset))
# (clustering.smk); only Systematic ID / A / DR / DL are read (never the cluster
# id), so any variant's table gives identical results.
# Two-part: module visualization (named complexes in DR/DL fitness space) +
# coherence permutation test (all complexes 3<=size<=300 whose DR>0.3 members
# form a tighter-than-random cluster). Weiszfeld geometric median + seeded MPD
# permutation test live in the shared workflow/src/coherence/metrics.py.
#
# DATA-PATH NOTE (deviation from the plan skeleton — established pattern):
#   * The plan skeleton used `resources/pombase/{pombase_version}/
#     macromolecular_complex.tsv`. The real PomBase convention in this repo is
#     `resources/external/pombase/{version}/ontologies_and_associations/
#     macromolecular_complex_annotation.tsv` (same layout features.smk /
#     noncoding_rna.smk use). The reference version comes from datasets.yaml
#     (DATASETS["reference"]["pombase_version"] = 2026-06-01), NOT a wildcard —
#     these are per-dataset rules with no {pombase_version} wildcard.
#   * final_clusters is built by the clustering stage for buildable variants; a
#     dry-run chains it off that rule (or reports the curated manual_merge tsv
#     missing if the selected variant is manual).

_CPLX_WORK = "results/complex/{dataset}/_work"


rule analyze_complex_modules:
    input:
        final_clusters=lambda wc: final_clusters_path(wc.dataset, selected_variant(wc.dataset)),
        complex_annotation=lambda wc: (
            f"resources/external/pombase/{DATASETS['reference']['pombase_version']}/"
            "ontologies_and_associations/macromolecular_complex_annotation.tsv"
        ),
    output:
        module_viz=f"{_CPLX_WORK}/module_visualization_done.flag",
        module_figure="results/complex/{dataset}/complex_module_visualization.pdf",
    params:
        modules=config.get("complex", {}).get("modules", {}),
    log:
        "logs/complex/analyze_complex_modules_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [complex] Visualizing named complex modules for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/complex/analyze_complex_modules.py \
            --final-clusters {input.final_clusters} \
            --complex-annotation {input.complex_annotation} \
            --modules "{params.modules}" \
            --output-flag {output.module_viz} \
            --output-figure {output.module_figure} &> {log}
        """


rule compute_complex_coherence:
    input:
        final_clusters=lambda wc: final_clusters_path(wc.dataset, selected_variant(wc.dataset)),
        complex_annotation=lambda wc: (
            f"resources/external/pombase/{DATASETS['reference']['pombase_version']}/"
            "ontologies_and_associations/macromolecular_complex_annotation.tsv"
        ),
    output:
        metrics="results/complex/{dataset}/complex_coherence_metrics.tsv",
        coherence_figure="results/complex/{dataset}/coherence_analysis.pdf",
    params:
        min_size=config.get("complex", {}).get("min_complex_size", 3),
        max_size=config.get("complex", {}).get("max_complex_size", 300),
        dr_threshold=config.get("complex", {}).get("dr_threshold", 0.3),
        n_permutations=config.get("complex", {}).get("n_permutations", 1000),
        random_state=config.get("complex", {}).get("random_state", 42),
    log:
        "logs/complex/compute_complex_coherence_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [complex] Computing coherence metrics for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/complex/compute_complex_coherence.py \
            --final-clusters {input.final_clusters} \
            --complex-annotation {input.complex_annotation} \
            --min-size {params.min_size} \
            --max-size {params.max_size} \
            --dr-threshold {params.dr_threshold} \
            --n-permutations {params.n_permutations} \
            --random-state {params.random_state} \
            --output-metrics {output.metrics} \
            --output-figure {output.coherence_figure} &> {log}
        """
