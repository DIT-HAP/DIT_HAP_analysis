# Worktree ↔ main compatibility reconciliation

**Date:** 2026-07-21
**Branch:** `worktree-migrate-remaining-notebooks`
**Fork point:** `eb0cef0` (`feat(validation): add FACS depletion curves notebook`)

## Goal

Make the migration branch compatible with the *current* `main` so it merges
cleanly and its Batch B rules run against main's clustering outputs.

**Hard constraint:** only the worktree branch is edited. `main`'s files —
especially `workflow/src/coherence/metrics.py` and everything in the clustering
finalize-variant system — are treated as read-only upstream contracts. We adapt
to them; we never modify them.

## Background: what main changed since the fork

Main advanced ~30 commits, dominated by two threads:

1. **Clustering finalize-variant system.** `final_clusters.tsv` moved from a
   single hand-curated `resources/curated/final_clusters.tsv` to **per-variant**
   paths, resolved by a helper in `clustering.smk`:
   - buildable variants → `results/clustering/final/{dataset}/{variant}/final_clusters.tsv`
   - `manual_merge` variant → `resources/curated/final_clusters/{dataset}/{variant}.tsv`
   - helper: `final_clusters_path(dataset, variant)`; selected variant per
     dataset: `selected_variant(dataset)` (from `config.clustering.selected_variant`
     + `selected_variant_overrides`).
   - The final cluster-id column was renamed repo-wide `revised_cluster` → `cluster`.

2. **Themes A + D verification/coherence.** Added `workflow/src/coherence/`
   (`metrics.py`, `attribution.py`) plus `notebooks/verification_complex/*.py`.
   `metrics.py` ports the SAME `complex_analysis.ipynb` coherence algorithm that
   this branch's Task 6 independently ported into `workflow/src/complex/coherence.py`.

## Incompatibility inventory

Textual overlap (changed on both sides): only `Snakefile`, `config/analysis.yaml`.
Everything else is semantic. Full sweep results:

| # | Issue | Severity | Resolution |
|---|-------|----------|------------|
| 1 | Task 5/6 rule inputs hardcode the OLD `resources/curated/final_clusters.tsv`, which no longer exists on main | **Breaks Batch B DAG** | Rewire inputs through `final_clusters_path(dataset, selected_variant(dataset))` |
| 2 | Task 5/6 docstrings say `revised_cluster` + old path | Cosmetic (scripts read only `Systematic ID`/`DR`/`DL`/`A`, never the cluster id) | Docstring-only cleanup |
| 3 | Duplicate coherence module (`src/complex/coherence.py` vs main's `src/coherence/metrics.py`) | Duplication, no functional break | Delete the worktree module; import main's; adapt caller + tests |

Confirmed NON-issues (swept, no action):
- `workflow/src/plotting/style.py` — main did not touch it; Task scripts' `COLORS`/`AX_*` imports are safe.
- `workflow/src/plotting/gene_level.py::plot_given_genes_on_feature_space` — main changed the file but this function is **byte-identical** to the fork; Task 6's `analyze_complex_modules.py` call site matches the signature.
- `resources/curated/complex_subunit_roles.tsv` (main-added) — Task 6 never references it.
- Snakefile include order: main includes `clustering.smk` at line 31, before the
  line-36+ block where this branch's 8 includes land, so `final_clusters_path` /
  `selected_variant` are defined before Task 5/6 use them.

## Numeric-equivalence evidence for #3

The two coherence implementations were compared offline on a synthetic
background (300 genes, 8-member complex, seed 42), with main's `bg` argument set
to the FULL background (members included), matching Task 6's null construction:

- `observed_mpd` == `pairwise_distance(members, "median")` (exact).
- `z_score` = -0.439671 and `p_value` = 0.353 matched to 6 decimals.
- `rng.choice(2D array)` and `rng.choice(indices)` draw identically for a fixed seed.

So swapping to main's module is numerically exact **provided** the caller keeps
passing the full DR>threshold background (members included) as `bg`.

Known API differences the caller/tests must absorb (main's module, unchanged):
- `compute_distance_zscore(X, bg, method, n_permutations, random_state)` returns
  a `(z, p)` **tuple**, not a dict. `observed_mpd` is obtained separately via
  `pairwise_distance(members, method="median")`.
- `coherence_metrics(X)` returns main's key names (`median_pairwise_distance`,
  `mean_distance_to_centroid`, `mean_knn_distance`, …) — not `mpd`/`median_distance`.
- main's `geometric_median` has no iteration cap and no `n<=1` guard; its null-std
  guard is `if std == 0`. These are main's choices; we do NOT change them. Task 6
  already guards `n_members <= 1` at the caller before calling into the module, so
  the degenerate path never reaches main's ungated code.

## Plan

Executed as commits stacked on top of a merge (not a rebase — one conflict
resolution vs. replaying 16 reviewed commits).

### Step 0 — merge main into the branch
`git merge main`. Resolve the two expected conflicts by union:
- **Snakefile**: keep main's rewritten `rule all` header/targets AND this branch's
  8-line include block (after `pcr_qc.smk`). Confirm the hardcoded
  `workdir: "/data/c/yangyusheng_optimized/DIT_HAP_analysis"` survives (the
  temporary-edit-for-dry-run / revert-before-commit rule still applies).
- **config/analysis.yaml**: keep main's expanded `clustering:` section AND this
  branch's `spikein:` / `comparison:` / `complex:` sections (disjoint — union).

### Step 1 — #1 rewire Batch B inputs (`fix`)
In `workflow/rules/comparison.smk` (1 site) and `workflow/rules/complex.smk`
(2 sites), replace `final_clusters="resources/curated/final_clusters.tsv"` with:
```python
final_clusters=lambda wc: final_clusters_path(wc.dataset, selected_variant(wc.dataset)),
```
Update the rule header comments that describe the old always-absent curated input.

### Step 2 — #3 adopt main's coherence module + #2 docstrings (`refactor`)
- Delete `workflow/src/complex/coherence.py` and `workflow/src/complex/__init__.py`
  IF nothing else imports them (verify: only `compute_complex_coherence.py` +
  the test do). Keep `workflow/src/complex/__init__.py` only if the package is
  still needed for other `complex/` modules.
- Rewrite `workflow/scripts/complex/compute_complex_coherence.py` to
  `from workflow.src.coherence.metrics import coherence_metrics, compute_distance_zscore, pairwise_distance`
  (+ `normalize_dr_dl` if it lets us drop local normalization). Adapt:
  - build `member_points` and pass the full background as `bg`;
  - call `compute_distance_zscore(member_points, bg, method="median_pairwise_distance", ...)`,
    unpack `(z, p)`;
  - compute `observed_mpd = pairwise_distance(member_points, method="median")`;
  - remap any `coherence_metrics` key reads to main's names.
  Preserve the caller-side `n_members <= 1` degenerate guard.
- Docstring cleanup (#2) in `compute_complex_coherence.py`,
  `analyze_complex_modules.py`, and `compare_large_scale_studies.py`:
  `revised_cluster` → `cluster`; example paths → per-variant path.
- Rewrite `tests/test_complex_coherence.py` to import and assert against main's
  API (tuple return, main's metric keys). Drop the `EPSILON` import (main's module
  doesn't export it as a public constant the same way).

### Step 3 — verify
1. `pytest tests/ --ignore=tests/test_features_assembly.py --ignore=tests/test_features_genome.py`
   green (baseline re-established after merge pulls in main's tests).
2. Task 6 end-to-end on real HD_DIT_HAP: confirm coherence z/p unchanged vs.
   pre-swap (regression confirmation of the offline equivalence proof).
3. Temporarily point `workdir` at the worktree; `snakemake -n` the three Batch B
   targets; confirm `final_clusters` resolves to
   `results/clustering/final/HD_DIT_HAP/kmeans_direct9/final_clusters.tsv`; revert `workdir`.
4. `git diff` sanity: no changes under `workflow/src/coherence/`, `clustering.smk`,
   or any other main-owned file.

### Commit sequence
1. merge commit (Snakefile/config union resolution)
2. `fix(batch-b): source final_clusters via per-variant clustering helper`
3. `refactor(complex): adopt shared src/coherence module, sync cluster-column docs`
