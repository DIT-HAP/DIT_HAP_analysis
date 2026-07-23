# Verification Plots Migration — Design

**Date:** 2026-07-22
**Branch:** worktree-migrate-remaining-notebooks
**Source notebook:** `DIT_HAP_pipeline/workflow/notebooks/compare_with_deletion_library.ipynb` §4-5

## Problem

The T3 migration of `compare_with_deletion_library.ipynb` into
`workflow/scripts/verification/compare_deletion_library.py` was a deliberate
"simplified port": it kept only the category donut + DR-by-category scatter,
and dropped the notebook's §4-5 plotting (boxplot+violin comparisons, the
critical-gene verification analysis, and the DIT_HAP-vs-gRNA depletion
curves). This design adds those back as reproducible, deterministic outputs.

Out of scope (stays notebook-only): the altair interactive charts
(`area_fraction_vs_DR`, `area_fraction_on_DR_vs_DL`) — they produce HTML and
are exploratory. The user scoped this migration to matplotlib figures only.

## Scope

Add to the verification stage, for `compare_with_deletion_library.ipynb` §4-5:

1. **Boxplot+violin comparisons** (§4.1): basic per-category DR distribution
   (6 canonical categories) + the "including small colonies (E)" variant.
2. **Critical-gene analysis** (§4.2-4.4): four outlier groups, each with a
   boxplot+violin, a donut of verification-result composition, and a per-group
   gene-detail TSV for human review.
3. **Depletion curves** (§5): per critical-group, DIT_HAP measured points +
   sigmoid (Gompertz) fit + inflection slope, with the HD gRNA curve overlaid
   where gRNA data exists.

## Data availability (verified)

- Basic + critical boxplots need only data the current `run()` already loads
  (gene_level DR/DL + `DeletionLibrary_essentiality`, deletion-library
  `Category`, curated verification table). No new load for these.
- Critical-gene TSVs need the verification table's area columns
  (`median/mean_area_day3..6`) — present in the current
  `essentiality_verification.csv` (411 rows, 20 cols).
- Depletion curves need per-timepoint LFC: DIT_HAP from
  `release/gene_level/gene_level_fitting_statistics.tsv` (has YES0-4 + A/DR/DL,
  every dataset); gRNA from `DIT_HAP_pipeline/resources/HD_gRNA_data.csv`
  (t0-t5, HD-only, in the pipeline repo).
- Critical-group hit counts under canonical-category mapping (HD_DIT_HAP):
  WT2nonWT 231 (224 verified), scE2E 26, sc2E 12, E2V 40 — none empty.

## Category schema drift (must handle)

> **SUPERSEDED (2026-07-22):** the `cat_canon` / `_display_category` folding
> described in this section was later removed by project decision. The
> verification stage now uses the RAW curated `Category` labels verbatim
> everywhere — display, ordering, boxplot/critical grouping, and the outlier
> filters all match the literal labels (`Category == 'WT-like'`,
> `Category in ['spores','germinated','microcolonies']`). Compound
> multi-phenotype labels are therefore NOT folded and do not enter the critical
> groups. Colors are the only place a fallback remains: `_category_color_key()`
> maps a raw label to its phenotype-family color for lookup only, never
> changing the shown text. The rest of this section is kept for history.

The current `deletion_library_categories.xlsx` uses a newer vocabulary than
the notebook: `WT` -> `WT-like`, plus compound labels (`spores, germinated`,
`microcolonies, small colonies`, ...). The notebook's raw filters
(`Category == 'WT'`, `Category in ['spores','germinated','microcolonies']`)
would match zero rows. Fix: derive a `cat_canon` column via the existing
`_display_category` alias map, then apply the notebook's filters against
`cat_canon`.

## gRNA alignment (must handle)

gRNA rows are indexed by gene **Name**, which is NOT unique (5257 rows /
5250 unique Names). The systematic ID embedded in `gRNA_ID`
(`"SPAC1002.02_42"` -> `"SPAC1002.02"`) IS unique (5257) and overlaps 4465 of
the 4513 DIT_HAP genes. So: **align DIT_HAP and gRNA on Systematic ID**
(extract gRNA's via `gRNA_ID.rsplit("_", 1)[0]`), and **title each panel with
gene_name** for readability. This avoids the duplicate-Name `.loc` crash the
notebook is exposed to.

## Architecture

**Reusable plotting primitives** (added to `workflow/src/plotting/`):

- `generic.py::boxplot_with_violinplot(labels, values, ax, colors)` — horizontal
  violin+box composite. Domain-agnostic, sits beside `donut_chart`.
- `gene_level.py::sigmoid_gompertz(x, A, DR, DL)` — the notebook's
  `sigmoid_function` renamed and re-parameterized to DR/DL (release columns;
  notebook used um/lam). `A==0` -> zeros; `alpha=(DR*e)/A`;
  `A*exp(-exp(clip(alpha*(DL-x)+1, -700, 700)))`.
- `gene_level.py::plot_gene_depletion_curve(ax, dit_row, grna_row, dit_gen,
  grna_gen, title)` — one gene per panel: DIT_HAP measured points (YES0-4),
  sigmoid fit line, inflection slope line, and (when `grna_row` is not None)
  the gRNA curve. `title` is the gene_name.

**Analysis layer** (in `compare_deletion_library.py`, reusing already-loaded
data in `run()`):

- `canonicalize_category(merged)` — add `cat_canon` via `_display_category`.
- `build_final_merged(gene_level, deletion_library, verification_full)` —
  reconstruct the notebook's `final_merged`: gene_level (FYPOviability, DR, DL)
  + deletion-library description columns + the full verification table (area
  day3-6). Feeds the critical-gene TSVs. Includes the gpd1 quirk (append
  `SPBC215.05` E/E to the simplified verification table, byte-faithful; the
  reason is documented inline).
- `prepare_verification_data(merged, final_merged, verification, filter, sort)`
  — pure function: select outliers by `filter` -> cross with verification ->
  bucket into `{verified_category: [DR...], "Not verified": [...]}`, and return
  the per-group gene-detail DataFrame. Unit-test focus.

Time-point constants (module-level, notebook-hardcoded):
DIT_HAP `[0.0, 2.352, 5.588, 9.104, 12.48]`, gRNA `[0.0, 4.8, 7.9, 11.4, 14.7,
18.]` — gRNA generations from HD_gRNA_data's `time_points`.

The four critical groups are declared once as `_CRITICAL_GROUPS`:

```python
_CRITICAL_GROUPS = {
    "WT2nonWT": {"filter": "cat_canon == 'WT' and DR > 0.35", "sort": "desc"},
    "scE2E":    {"filter": "cat_canon == 'small colonies' and DR > 0.75 and DeletionLibrary_essentiality == 'E'", "sort": "desc"},
    "sc2E":     {"filter": "cat_canon == 'small colonies' and DR > 0.75 and DeletionLibrary_essentiality != 'E'", "sort": "desc"},
    "E2V":      {"filter": "cat_canon in ['spores','germinated','microcolonies'] and DR < 0.35", "sort": "asc"},
}
```

## Outputs

New outputs, kept separate from the existing `deletion_library_comparison.pdf`
and `verification_stats.tsv` (those are untouched):

1. `verification_boxplots.pdf` — multi-page: basic boxplot (6 canonical
   categories) + the small-colonies-(E) variant + one boxplot per critical
   group + one donut per critical group. All four groups get a donut (the
   notebook only drew three; completeness over notebook fidelity, per user).
2. `verification_depletion_curves.pdf` — one section per critical group;
   within a group the outlier genes are laid out in a 4-column grid; each
   panel is DIT_HAP points + sigmoid fit + inflection slope + (HD) gRNA
   overlay; panel title = gene_name.
3. `critical_genes_{group}.tsv` x4 — per-group gene detail (Systematic ID,
   Name, DR, DL, FYPOviability, verification result, area day3-6) for human
   review.

Genes plotted in the depletion-curve PDF are exactly each group's filtered
outliers (deduped, DR-sorted) — the same genes as that group's boxplot/TSV,
fully reproducible. The notebook's hand-picked lists and `noised_genes`
exclusions are NOT carried over.

## Rule wiring (`verification.smk`, worktree-owned; main untouched)

- New inputs: `gene_timepoints` =
  `{repo}/{release_dir}/gene_level/gene_level_fitting_statistics.tsv`;
  `grna_timepoints` via `lambda wc: _GRNA_TIMEPOINT_DATA.get(wc.dataset)` where
  `_GRNA_TIMEPOINT_DATA = {"HD_DIT_HAP": "<pipeline repo>/resources/HD_gRNA_data.csv"}`
  (same per-dataset-map pattern as noncoding_rna.smk's `_NONCODING_FITTING`).
  When a dataset has no gRNA entry the input is omitted and curves render
  DIT_HAP-only.
- New outputs: `boxplots`, `depletion_curves`, and the four
  `critical_genes_{group}.tsv`.
- Script argparse grows matching `--gene-timepoints`, `--grna-timepoints`
  (optional), `--output-boxplots`, `--output-depletion-curves`,
  `--output-critical-genes-dir`.

## Testing

Append to `tests/test_verification.py` (reuse existing fixtures):

- `test_canonicalize_category` — `WT-like`->`WT`, compound labels fold to the
  canonical bucket.
- `test_prepare_verification_data_buckets` — synthetic outliers bucket into
  `{verified_category: [DR...], "Not verified": [...]}` correctly.
- `test_prepare_verification_data_empty_group` — a zero-hit filter returns an
  empty dict without raising.
- `test_sigmoid_gompertz` — `A==0` -> zeros; known A/DR/DL -> expected values;
  out-of-range exponent is clipped (no overflow).
- `test_boxplot_with_violinplot_smoke` — returns an Axes with the expected
  number of y-ticks for given labels/values.
- `test_grna_sysid_extraction` — `gRNA_ID.rsplit("_", 1)[0]` yields the
  Systematic ID; alignment is on Systematic ID, not the non-unique Name.

Pre-existing ignored collection errors are unchanged (codonbias x2, main's
`test_clustering.py`). Verify the full suite stays green plus the new tests,
and dry-run `results/verification/HD_DIT_HAP/verification_boxplots.pdf` (and the
other new targets) via the temporary-workdir trick, reverting before commit.

## Implementation order

1. Add plotting primitives (`boxplot_with_violinplot`, `sigmoid_gompertz`,
   `plot_gene_depletion_curve`) + their unit tests.
2. Add analysis functions (`canonicalize_category`, `build_final_merged`,
   `prepare_verification_data`) + unit tests.
3. Wire the three new figure/TSV builders into `run()`.
4. Update `verification.smk` inputs/outputs + argparse.
5. Dry-run all new targets; run full test suite; commit.

