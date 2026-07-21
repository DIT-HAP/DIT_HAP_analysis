# =============================================================================
# complex.smk — Macromolecular complex coherence analysis
# =============================================================================
#
# Batch B (requires resources/curated/final_clusters.tsv — the un-buildable,
# human-curated cluster table from the manual finalize_gene_clusters step).
# Two-part: module visualization (named complexes in DR/DL fitness space) +
# coherence permutation test (all complexes 3<=size<=300 whose DR>0.3 members
# form a tighter-than-random cluster). Weiszfeld geometric median + seeded MPD
# permutation test live in workflow/src/complex/coherence.py.
#
# DATA-PATH NOTE (deviation from the plan skeleton — established pattern):
#   * The plan skeleton used `resources/pombase/{pombase_version}/
#     macromolecular_complex.tsv`. The real PomBase convention in this repo is
#     `resources/external/pombase/{version}/ontologies_and_associations/
#     macromolecular_complex_annotation.tsv` (same layout features.smk /
#     noncoding_rna.smk use). The reference version comes from datasets.yaml
#     (DATASETS["reference"]["pombase_version"] = 2026-06-01), NOT a wildcard —
#     these are per-dataset rules with no {pombase_version} wildcard.
#   * final_clusters.tsv is absent until the manual finalize step runs; a
#     dry-run reports it as a missing input (expected Batch-B state).

_CPLX_WORK = "results/complex/{dataset}/_work"


rule analyze_complex_modules:
    input:
        final_clusters="resources/curated/final_clusters.tsv",
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
        final_clusters="resources/curated/final_clusters.tsv",
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
