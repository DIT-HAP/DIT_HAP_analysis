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

See `docs/plans/` for design docs and implementation plans.
