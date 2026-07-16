# DIT-HAP Analysis

Downstream analysis of DIT-HAP depletion data: gene feature collection, enrichment,
clustering, ML, and thesis figures. Consumes packaged `release/` outputs from
[DIT_HAP_snakemake](../DIT_HAP_snakemake/) via `config/datasets.yaml`.

## Structure

- `config/datasets.yaml` — registry pointing at DIT_HAP_snakemake's per-project `release/` dirs
- `config/analysis.yaml` — this project's own analysis parameters
- `workflow/src/` — shared library: `data_config.py`, `io.py`, `gene_ids.py`, `plotting/`, `enrichment/`, `features/`
- `workflow/rules/` — Snakemake rule files per analysis stage
- `workflow/scripts/` — deterministic per-rule scripts (python-script-conventions)
- `workflow/envs/` — conda environment YAMLs per rule
- `notebooks/` — human-judgment analyses with explicit input/output contracts (see header of each notebook)
- `results/{stage}/` — Snakemake-produced, semantically named, safe to delete and rerun
- `resources/curated/` — human-curated artifacts, version-controlled, NOT reproducible by rerunning Snakemake
- `Snakefile` — entry point

## Requirements

- Python 3.12
- Snakemake 9.0+
- Conda/mamba for environment management

## Usage

```bash
# Activate Snakemake env
mamba activate snakemake

# Specific rule
snakemake --cores 8 --use-conda results/features/2025-10-01/pombe_coding_gene_protein_features.tsv
```

## Core analysis chain (clustering → enrichment / ml)

The chain is not a single DAG: it has one deliberate manual step (cluster merge),
which is a human-judgment decision kept as a notebook (design doc §5).

```
1. Feature collection (dataset-independent, by PomBase version)
   snakemake --use-conda results/features/<version>/pombe_coding_gene_protein_features.tsv

2. Candidate clustering (per dataset, deterministic)
   snakemake --use-conda results/clustering/candidates/<dataset>/candidate_clusters.tsv

3. MANUAL: finalize clusters (human decision — merges 64 candidates → 9)
   Run notebooks/clustering/finalize_gene_clusters.ipynb (set DATASET at top),
   review the feature-space plots, adjust the merge dicts, and write
   resources/curated/final_clusters.tsv  (version-controlled, un-buildable input).

4a. Enrichment (per dataset; needs final_clusters.tsv)
    snakemake --use-conda results/enrichment/raw/<dataset>/<version>/go_enrichment_full_filtered.tsv

4b. ML AutoML (per dataset x target x mode; needs final_clusters.tsv)
    snakemake --use-conda results/ml/models/<dataset>/<version>/um_Explain/metrics.tsv

Optional (hits STRING/REVIGO web APIs; cached under resources/external/enrichment_cache/):
    snakemake --use-conda results/enrichment/network/<dataset>/<version>/go_enrichment_full_revigo.tsv
```

If a rule reports `Missing input files: resources/curated/final_clusters.tsv`,
run the manual finalize notebook (step 3) first — this is intentional.

See `docs/plans/` for design docs and implementation plans.
