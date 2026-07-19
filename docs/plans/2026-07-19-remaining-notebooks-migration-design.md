# Remaining Notebooks Migration Design

**Date**: 2026-07-19
**Scope**: Migrate 11 remaining source notebooks from `DIT_HAP_pipeline/workflow/notebooks/` into the current `DIT_HAP_analysis` project structure.
**Status**: Validated, ready for implementation.

---

## 1. Migration map

| Source notebook | New location | Type |
|---|---|---|
| `spike_in.ipynb` | `workflow/rules/spikein.smk` + `workflow/scripts/spikein/` | Deterministic script |
| `gene_coverage_analysis.ipynb` | `workflow/rules/coverage.smk` + `workflow/scripts/coverage/` | Deterministic script |
| `compare_with_deletion_library.ipynb` | `workflow/rules/verification.smk` + `workflow/scripts/verification/` | Deterministic script |
| `non_coding_RNA_analysis.ipynb` | `workflow/rules/noncoding_rna.smk` + `workflow/scripts/noncoding_rna/` | Deterministic script |
| `compare_with_other_large_scale_studies.ipynb` | `workflow/rules/comparison.smk` + `workflow/scripts/comparison/` | Deterministic script |
| `complex_analysis.ipynb` | `workflow/rules/complex.smk` + `workflow/scripts/complex/` | Deterministic script |
| `upstream_and_downstream_analysis.ipynb` | `workflow/rules/utr.smk` + `notebooks/domain_analysis/review_utr_insertions.ipynb` | Split |
| `genes_with_domain_differences.ipynb` | `workflow/rules/domain_differences.smk` + `notebooks/domain_analysis/review_domain_differences.ipynb` | Split |
| `further_analysis_based_on_enrichment.ipynb` | `notebooks/enrichment/further_analysis.ipynb` | Human notebook |
| `thesis_figures.ipynb` | `notebooks/thesis_figures/figures.ipynb` | Human notebook |
| `visualize_dit_hap_style.ipynb` | `notebooks/reference/visualize_dit_hap_style.ipynb` | Reference only |

---

## 2. Execution batches

### Batch A — no `final_clusters.tsv` dependency

These can be run as soon as upstream pipeline `release/` dirs are available.

**`spikein.smk`**
- Input: `{snakemake_repo}/projects/Spikein/release/insertion_level/raw_reads.filtered.tsv`
- Output: `results/spikein/spike_in_correlation.pdf`, `results/spikein/spike_in_stats.tsv`
- No `{dataset}` wildcard — Spikein is a standalone project.
- Logic: extract 5 known spike-in coordinates, compute per-sample read counts, linear regression against known spike-in ratios, produce correlation plot.

**`coverage.smk`**
- Input: `{release}/insertion_level/fitting_results.tsv`, `{release}/insertion_level/annotations.tsv`
- Output: `results/coverage/{dataset}/coverage_stats.tsv`, `results/coverage/{dataset}/coverage_figures.pdf`
- Logic: count total vs. in-gene insertions (filter: `Type != 'Intergenic region' and Distance_to_stop_codon > 4`), produce donut chart per chromosome and overall.
- Note: original notebook reads `insertion_density_analysis.tsv` from reports — this intermediate is produced by the upstream pipeline's `insertion_density_analysis` rule. Port the relevant computation directly rather than depending on that report file.

**`verification.smk`**
- Input: `{release}/gene_level/fitting_results.tsv`, `resources/curated/deletion_library_categories.xlsx`, `resources/curated/essentiality_verification.csv`
- Output: `results/verification/{dataset}/deletion_library_comparison.pdf`, `results/verification/{dataset}/verification_stats.tsv`
- Logic: merge DIT-HAP gene-level results with deletion library phenotype categories and essentiality verification; plot donut charts + scatter comparison.

**`noncoding_rna.smk`**
- Input: non-coding RNA gene-level statistics from upstream pipeline (`{release}/gene_level/noncoding_rna_fitting_results.tsv` — path to be confirmed with `data_config`), PomBase ncRNA bed, GtRNAdb bed, Marguerat 2012 mRNA abundance Excel.
- Output: `results/noncoding_rna/{dataset}/ncrna_analysis.pdf`, `results/noncoding_rna/{dataset}/ncrna_stats.tsv`
- Logic: merge ncRNA stats with tRNA annotations from GtRNAdb, analyze tRNA vs. other ncRNA depletion, compare with mRNA abundance.

### Batch B — requires `resources/curated/final_clusters.tsv`

**`comparison.smk`**
- Input: `resources/curated/final_clusters.tsv`, `results/features/{pombase_version}/pombe_coding_gene_protein_features.tsv`, gRNA data path (register in `analysis.yaml` as `gRNA_data_file`)
- Output: `results/comparison/{dataset}/pairwise_fitness_comparison.pdf`, `results/comparison/{dataset}/fitness_correlation_stats.tsv`
- Logic: merge DIT-HAP cluster data with gRNA fitted parameters, Barseq (dulab + Koch), integration density, colony size; pairwise scatter plots with KDE overlay; Pearson r for each pair.

**`complex.smk`** (contains coherence analysis)
- Input: `resources/curated/final_clusters.tsv`, PomBase macromolecular complex annotation TSV, gRNA data
- Output:
  - `results/complex/{dataset}/complex_coherence_metrics.tsv`
  - `results/complex/{dataset}/complex_module_visualization.pdf`
  - `results/complex/{dataset}/coherence_analysis.pdf`
- Logic (2 parts):
  1. **Module visualization** (section 4 of source): for named complexes (cytoplasmic translation, kinetochore, mitochondria, vesicle, vacuolar ATPase), plot depletion curves and feature-space positions using `workflow/src/plotting/gene_level.py`.
  2. **Coherence analysis** (section 5): Weiszfeld weighted geometric median algorithm to find the multivariate center of each complex in (DR, DL, A) space; median-polished pairwise distance (MPD) as coherence metric; permutation test for significance; FDR correction (Benjamini-Hochberg). Write `complex_coherence_metrics.tsv`.

### Batch C — split analyses (insertion-level + manual review)

**`utr.smk`** + `notebooks/domain_analysis/review_utr_insertions.ipynb`
- Deterministic part: UTR insertion classification (`assign_UTR_type` function), aggregate UTR insertion statistics per gene, write `results/utr/{dataset}/utr_insertion_stats.tsv`
- Human notebook: loads stats, plots intragenic insertion visualizations for candidate non-essential domain genes using `workflow/src/plotting/gene_level.py::intragenic_insertion_visualization`
- Input for rule: insertion-level fitting results + annotations + gene-level results
- Input for notebook: `results/utr/{dataset}/utr_insertion_stats.tsv` + `resources/curated/non_essential_domain_candidates.xlsx`

**`domain_differences.smk`** + `notebooks/domain_analysis/review_domain_differences.ipynb`
- Deterministic part: select genes with significant intra-gene DR heterogeneity (`DR > 0.15` threshold), compute insertion fraction statistics relative to start/stop codon, write `results/domain_differences/{dataset}/domain_candidate_stats.tsv`
- Human notebook: loads candidates, plots rank-order plots, histograms, and intragenic visualization for curated gene list

### Batch D — pure human notebooks

**`notebooks/enrichment/further_analysis.ipynb`**
- Reads: `results/enrichment/raw/{dataset}/` + `resources/curated/enrichment_term_categorization.xlsx`
- Logic: cluster-level enrichment bar charts, paralog cluster analysis, distribution plots using `workflow/src/plotting/gene_level.py::distribution_bar_for_given_genes`

**`notebooks/thesis_figures/figures.ipynb`**
- Reads: `results/{clustering,enrichment,ml,comparison,complex}/{dataset}/`
- Produces: `reports/thesis/*.pdf` for manuscript figures

**`notebooks/reference/visualize_dit_hap_style.ipynb`**
- Copy as-is, update paths to use `workflow/src/plotting/style.py`

---

## 3. New files to create

```
workflow/rules/
├── spikein.smk
├── coverage.smk
├── verification.smk
├── noncoding_rna.smk
├── comparison.smk
├── complex.smk
├── utr.smk
└── domain_differences.smk

workflow/scripts/
├── spikein/run_spikein_analysis.py
├── coverage/compute_coverage_stats.py
├── verification/compare_deletion_library.py
├── noncoding_rna/analyze_noncoding_rna.py
├── comparison/compare_large_scale_studies.py
└── complex/
    ├── analyze_complex_modules.py
    └── compute_complex_coherence.py

workflow/scripts/
├── utr/classify_utr_insertions.py
└── domain_differences/compute_domain_stats.py

notebooks/
├── domain_analysis/
│   ├── review_utr_insertions.ipynb
│   └── review_domain_differences.ipynb
├── enrichment/
│   └── further_analysis.ipynb
├── thesis_figures/
│   └── figures.ipynb
└── reference/
    └── visualize_dit_hap_style.ipynb
```

Snakefile: add includes for all 8 new `.smk` files.

`config/analysis.yaml`: add `gRNA_data_file` path and `complex_modules` dict (named gene lists for cytoplasmic translation / kinetochore / mitochondria / vesicle / vacuolar ATPase).

---

## 4. Quirks to preserve from source notebooks

- **coverage**: `in_gene_filter = "Type != 'Intergenic region' and Distance_to_stop_codon > 4"` — exact string.
- **verification**: color map for categories (`WT`, `E`, `spores`, etc.) uses specific COLORS indices — keep in `workflow/src/plotting/style.py`.
- **complex coherence**: Weiszfeld algorithm uses component-wise median as initialization; zero-distance guard (`replace any zero distances with small value`); MPD = median of all pairwise distances within complex projected to DR/DL axes.
- **comparison**: `clip(upper=200)` on integration density and ipkm/uipkm columns before plotting.
- **utr**: `distance_threshold=400` bp for UTR classification; `assign_UTR_type` handles strand-aware left/right gene boundary logic.
- **noncoding_rna**: merges GtRNAdb by chr+start+end (not by name), replaces `chrI/II/III` → `I/II/III`.
- **spikein**: 5 hardcoded spike-in coordinates (DY215/217/218/339/348) → move to `config/analysis.yaml`.

---

## 5. Dependency order for implementation

```
Batch A (parallel, no cluster dependency):
  spikein → coverage → verification → noncoding_rna

Batch B (after final_clusters.tsv exists):
  comparison → complex

Batch C (after insertion-level results exist):
  utr → domain_differences (deterministic parts)
  → then human notebooks

Batch D (after results from enrichment/comparison/complex exist):
  further_analysis → thesis_figures → reference
```
