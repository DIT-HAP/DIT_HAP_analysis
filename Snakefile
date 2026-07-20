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
include: "workflow/rules/spikein.smk"

# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
# Following the repo convention, per-stage targets are listed but commented —
# uncomment (or pass on the CLI) to run a specific stage. The core chain is:
#   clustering -> [manual: finalize_gene_clusters.ipynb] -> resources/curated/final_clusters.tsv
#   -> enrichment / ml. The manual finalize step is why the chain is not one DAG.
_REF = DATASETS["reference"]["pombase_version"]
_DATASET = DATASETS["default_dataset"]

rule all:
    input:
        # f"results/features/{_REF}/pombe_coding_gene_protein_features.tsv",
        # Clustering candidates (per dataset):
        f"results/clustering/candidates/{_DATASET}/candidate_clusters.tsv",
        # Enrichment (needs resources/curated/final_clusters.tsv from the manual notebook):
        # f"results/enrichment/raw/{_DATASET}/{_REF}/go_enrichment_full_filtered.tsv",
        # Network enrichment (optional, hits STRING/REVIGO — run explicitly):
        # f"results/enrichment/network/{_DATASET}/{_REF}/go_enrichment_full_revigo.tsv",
        # ML AutoML (target x mode):
        # f"results/ml/models/{_DATASET}/{_REF}/DR_Explain/metrics.tsv",
        # f"results/ml/models/{_DATASET}/{_REF}/DL_Explain/metrics.tsv",
        # PCR / library-prep QC figure (no dataset wildcard):
        # "results/pcr_qc/PCR_quality_control.pdf",
    message:
        "*** DIT-HAP analysis complete"
