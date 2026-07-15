# =============================================================================
# Snakefile — DIT-HAP analysis entry point
# =============================================================================

from snakemake.utils import min_version
from pathlib import Path
import yaml

min_version("9.0")

workdir: "/data/c/yangyusheng_optimized/DIT_HAP_analysis"

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
# include: "workflow/rules/features.smk"
# include: "workflow/rules/enrichment.smk"
# include: "workflow/rules/clustering.smk"
# (Phase 1 delivers only features.smk; remaining rules are follow-up work)

# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
rule all:
    input:
        # Uncommented in Task 5 once features.smk exists:
        # f"results/features/{DATASETS['reference']['pombase_version']}/pombe_coding_gene_protein_features.tsv",
    message:
        "*** DIT-HAP analysis complete"
