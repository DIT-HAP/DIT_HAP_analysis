# =============================================================================
# Snakefile — DIT-HAP analysis entry point
# =============================================================================

from snakemake.utils import min_version
from pathlib import Path
import yaml

min_version("9.0")

workdir: "/data/c/yangyusheng_optimized/DIT_HAP_analysis"

# This project's analysis parameters (clustering k, enrichment thresholds, ml
# splits, ...). Unlike datasets.yaml (a data registry read directly below), these
# ARE experiment parameters, so they flow through Snakemake's `config` object.
configfile: "config/analysis.yaml"

# ---------------------------------------------------------------------------
# Dataset registry (not a Snakemake configfile — see design doc §8)
# ---------------------------------------------------------------------------
with open("config/datasets.yaml") as f:
    DATASETS = yaml.safe_load(f)

wildcard_constraints:
    dataset="|".join(DATASETS["datasets"].keys()),

# ---------------------------------------------------------------------------
# Includes
# ---------------------------------------------------------------------------
# Rule filenames carry a two-level prefix: <chapter><section>, e.g. 1a_pcr_qc.
# The chapter is the analysis block, the section the position inside it — so
# `ls workflow/rules/` reads in analysis order. Chapters are narrative order,
# NOT dependency depth: same-chapter files have no order between them.
#
#   1  data plausibility + reference layer
#   2  characterisation that does not need clusters
#   3  clustering and what depends on it
#   9  terminal (annotate)
#
# Included in prefix order. clustering.smk must stay ahead of enrichment /
# ml / comparison — they call its selected_variant() / final_clusters_path().
include: "workflow/rules/1a_pcr_qc.smk"
include: "workflow/rules/1b_features.smk"
include: "workflow/rules/2a_coverage.smk"
include: "workflow/rules/annotate.smk"
include: "workflow/rules/clustering.smk"
include: "workflow/rules/enrichment.smk"
include: "workflow/rules/enrichment_network.smk"
include: "workflow/rules/ml.smk"
include: "workflow/rules/verification.smk"
include: "workflow/rules/noncoding_rna.smk"
include: "workflow/rules/comparison.smk"
include: "workflow/rules/coherence.smk"
include: "workflow/rules/utr.smk"
include: "workflow/rules/domain_differences.smk"

# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
# Following the repo convention, per-stage targets are listed but commented —
# uncomment (or pass on the CLI) to run a specific stage. The core chain is:
#   clustering spine -> finalize VARIANT -> enrichment / ml.
# Finalize has named variants (config.clustering.variants), each clustering for
# itself (no fixed candidate stage). Every variant produces
# results/clustering/{dataset}/{variant}/final_clusters.tsv (+ metrics.tsv):
# direct/auto_merge/grid via deterministic scripts, manual_merge via
# finalize_manual_merge (executes notebooks/clustering/finalize_gene_clusters.ipynb
# headlessly). enrichment fans out per variant; ml uses config.clustering.selected_variant.
_REF = DATASETS["reference"]["pombase_version"]
_DATASET = DATASETS["default_dataset"]
_SELECTED_VARIANT = config["clustering"]["selected_variant"]

rule all:
    input:
        # f"results/features/{_REF}/pombe_coding_gene_protein_features.tsv",
        # Selected finalize variant's clusters (buildable variants only):
        # f"results/clustering/{_DATASET}/{_SELECTED_VARIANT}/final_clusters.tsv",
        # Compare ALL buildable variants (builds every variant + a metrics table):
        # f"results/clustering/{_DATASET}/variant_metrics_comparison.tsv",
        # f"results/clustering/{_DATASET}/all_variants_cluster_scatter.pdf",
        # Enrichment (per variant):
        # f"results/enrichment/raw/{_DATASET}/{_SELECTED_VARIANT}/{_REF}/go_enrichment_full_filtered.tsv",
        # Network enrichment (optional, hits STRING/REVIGO — run explicitly):
        # f"results/enrichment/network/{_DATASET}/{_SELECTED_VARIANT}/{_REF}/go_enrichment_full_revigo.tsv",
        # ML AutoML (target x mode; uses selected_variant):
        # f"results/ml/models/{_DATASET}/{_REF}/DR_Explain/metrics.tsv",
        # f"results/ml/models/{_DATASET}/{_REF}/DL_Explain/metrics.tsv",
        # Library-prep QC (no dataset wildcard; spike-in stats feed the figure's
        # panel d, so both targets come from pcr_qc.smk):
        # "results/pcr_qc/PCR_quality_control.pdf",
        # "results/pcr_qc/spike_in_stats.tsv",
        f"results/coverage/{_DATASET}/coverage_stats.tsv",
        f"results/coverage/{_DATASET}/detailed_genes.xlsx",
        # Six single-figure PDFs (cnsplots); see 2a_coverage.smk for the list.
        f"results/coverage/{_DATASET}/coverage_overview.pdf",
        # f"results/coverage/{_DATASET}/coverage_by_deletion_viability.pdf",
        # f"results/coverage/{_DATASET}/coverage_by_characterisation.pdf",
        # f"results/coverage/{_DATASET}/coverage_insertion_placement.pdf",
        # f"results/coverage/{_DATASET}/coverage_dr_by_essentiality.pdf",
        # f"results/coverage/{_DATASET}/coverage_dl_by_essentiality.pdf",
        # f"results/verification/{_DATASET}/verification_stats.tsv",
        # f"results/verification/{_DATASET}/verification_boxplots.pdf",
        # f"results/verification/{_DATASET}/verification_depletion_curves.pdf",
        # f"results/noncoding_rna/{_DATASET}/ncrna_stats.tsv",
        # Batch B (requires resources/curated/final_clusters.tsv):
        # f"results/comparison/{_DATASET}/fitness_correlation_stats.tsv",
        f"results/coherence/{_DATASET}/coherence_terms_representatives.tsv",
        # Incoherence attribution (why complexes are dispersed) for physical-complex sources:
        expand(f"results/coherence/{_DATASET}/{{source}}/incoherence_attribution.tsv",
               source=config["coherence"].get("attribution_sources", ["go_macrocomplex"])),
        # Batch C (requires insertion-level results):
        # f"results/utr/{_DATASET}/utr_insertion_stats.tsv",
        # f"results/domain_differences/{_DATASET}/domain_candidate_stats.tsv",
    message:
        "*** DIT-HAP analysis complete"
