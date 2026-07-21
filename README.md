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

The finalize step (→ 9 final clusters) is a set of named **variants**, one per
strategy, configured under `config/analysis.yaml` → `clustering.variants`. Each
variant declares a `type`:

| type | how it makes 9 clusters | buildable by Snakemake? |
|------|-------------------------|-------------------------|
| `direct` | cluster fresh to k=9 with `method` | yes |
| `auto_merge` | ward-merge the method's k=64 candidate centroids down to 9 | yes |
| `grid` | axis cuts (`dr_cuts`/`dl_cuts`) on scaled DR/DL form a 9-cell grid | yes |
| `manual_merge` | human 64→9 merge in the notebook | no (curated input) |

Every variant emits the same contract: a `cluster` column with the 1..9 labels
(lowest-mean-DR group = WT = 9). `auto_merge`/`manual_merge` also keep a
`raw_cluster` column (pre-merge labels). Enrichment fans out over **all** variants
so you can compare them; ml/thesis use the single `clustering.selected_variant`
(per-dataset overridable). Numbering is always by mean DR — no variant hand-assigns
final ids (design doc `2026-07-21-clustering-finalize-variants`).

```
1. Feature collection (dataset-independent, by PomBase version)
   snakemake --use-conda results/features/<version>/pombe_coding_gene_protein_features.tsv

2. Candidate clustering (per dataset, deterministic — 64 candidates, 4 methods)
   snakemake --use-conda results/clustering/candidates/<dataset>/candidate_clusters.tsv

3. Finalize clusters (→ 9), per variant
   buildable variants (direct / auto_merge / grid) — no manual step:
     snakemake --use-conda results/clustering/final/<dataset>/<variant>/final_clusters.tsv
   manual_merge variant — human decision via
     notebooks/clustering/finalize_gene_clusters.ipynb (set DATASET/METHOD/VARIANT
     at top), review the feature-space plots, adjust the one merge dict, and write
     resources/curated/final_clusters/<dataset>/<variant>.tsv (version-controlled).

4a. Enrichment (per dataset x variant; needs that variant's final_clusters.tsv)
    snakemake --use-conda results/enrichment/raw/<dataset>/<variant>/<version>/go_enrichment_full_filtered.tsv

4b. ML AutoML (per dataset x target x mode; uses selected_variant's final_clusters.tsv)
    snakemake --use-conda results/ml/models/<dataset>/<version>/DR_Explain/metrics.tsv

Optional (hits STRING/REVIGO web APIs; cached under resources/external/enrichment_cache/):
    snakemake --use-conda results/enrichment/network/<dataset>/<variant>/<version>/go_enrichment_full_revigo.tsv
```

For a `manual_merge` variant only, if a rule reports `Missing input files:
resources/curated/final_clusters/<dataset>/<variant>.tsv`, run the finalize
notebook (step 3) first — this is intentional. The buildable variant types have no
manual step: their final clusters are buildable Snakemake targets.

See `docs/plans/` for design docs and implementation plans.
