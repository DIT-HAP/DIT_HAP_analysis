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
include: "workflow/rules/features.smk"
include: "workflow/rules/clustering.smk"
include: "workflow/rules/enrichment.smk"
include: "workflow/rules/enrichment_network.smk"
include: "workflow/rules/ml.smk"
include: "workflow/rules/pcr_qc.smk"

# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
# Following the repo convention, per-stage targets are listed but commented —
# uncomment (or pass on the CLI) to run a specific stage. The core chain is:
#   clustering spine -> finalize VARIANT -> enrichment / ml.
# Finalize has named variants (config.clustering.variants), each clustering for
# itself (no fixed candidate stage): buildable ones (direct/auto_merge/grid) produce
# results/clustering/{dataset}/{variant}/final_clusters.tsv (+ metrics.tsv);
# manual_merge is the curated resources/curated/final_clusters/{dataset}/{variant}.tsv
# (run notebooks/clustering/finalize_gene_clusters.ipynb). enrichment fans out per
# variant; ml uses config.clustering.selected_variant.
_REF = DATASETS["reference"]["pombase_version"]
_DATASET = DATASETS["default_dataset"]
_SELECTED_VARIANT = config["clustering"]["selected_variant"]

rule all:
    input:
        # f"results/features/{_REF}/pombe_coding_gene_protein_features.tsv",
        # Selected finalize variant's clusters (buildable variants only):
        # f"results/clustering/{_DATASET}/{_SELECTED_VARIANT}/final_clusters.tsv",
        # Compare ALL buildable variants (builds every variant + a metrics table):
        f"results/clustering/{_DATASET}/variant_metrics_comparison.tsv",
        # Enrichment (per variant; manual_merge needs its curated tsv first):
        # f"results/enrichment/raw/{_DATASET}/{_SELECTED_VARIANT}/{_REF}/go_enrichment_full_filtered.tsv",
        # Network enrichment (optional, hits STRING/REVIGO — run explicitly):
        # f"results/enrichment/network/{_DATASET}/{_SELECTED_VARIANT}/{_REF}/go_enrichment_full_revigo.tsv",
        # ML AutoML (target x mode; uses selected_variant):
        # f"results/ml/models/{_DATASET}/{_REF}/DR_Explain/metrics.tsv",
        # f"results/ml/models/{_DATASET}/{_REF}/DL_Explain/metrics.tsv",
        # PCR / library-prep QC figure (no dataset wildcard):
        # "results/pcr_qc/PCR_quality_control.pdf",
    message:
        "*** DIT-HAP analysis complete"
