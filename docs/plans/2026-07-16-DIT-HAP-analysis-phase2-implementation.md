# DIT-HAP Analysis Phase 2 Implementation Plan — Core Chain (clustering → enrichment → ml)

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Port the core analysis chain — `gene_level_clustering`, `comprehensive_enrichment_analysis`, `machine_learning_data_preparation`, `machine_learning_analysis` — from `DIT_HAP_pipeline/workflow/notebooks/` into deterministic Snakemake rules + scripts under `DIT_HAP_analysis/`, following design doc §5/§7/§8. Deliver `clustering.smk`, `enrichment.smk`, `ml.smk`, the supporting `workflow/src/` module split, and the manual `finalize_gene_clusters.ipynb` notebook that produces the hub artifact `resources/curated/final_clusters.tsv`.

**Architecture:** Snakemake 9.0+ with conda env-per-rule. Notebook logic ported byte-faithfully (preserving documented quirks: custom clustering scaling, WT-cluster=9 branching, `revised_cluster` column). Deterministic steps enter the DAG; human-judgment steps (cluster merge) stay as notebooks connected by explicit input/output contracts.

**Tech Stack:** Python 3.12, Snakemake 9.13, pandas 3.x, scikit-learn, scipy, goatools, mljar-supervised, loguru, pytest 8.4.

## Binding decisions (from user, 2026-07-16)

1. **No seeding of `final_clusters.tsv`** — strict dependency order. `clustering.smk` produces candidate labels → human runs `finalize_gene_clusters.ipynb` → `resources/curated/final_clusters.tsv` → then enrichment/ml consume it. Do NOT copy the existing `kmeans_cluster_result.tsv` as a shortcut.
2. **Split network enrichment** — goatools GO/FYPO/MONDO is the deterministic main rule (`enrichment.smk`, pure-local, in the DAG). STRING-db + REVIGO go into a separate optional rule (`enrichment_network.smk`) with cached responses.
3. **Both mljar modes** — Explain AND Perform. Perform sets an explicit large `total_time_limit` so the full algorithm list always completes; pin `random_state` and the mljar version.

## Data availability constraint

`HD_DIT_HAP_generationRAW` (the notebooks' hardcoded dataset) has **no packaged `release/` dir yet**. Only `LD_DIT_HAP_generationRAW`, `LD_haploid`, `Spore2YES6_1328` have `release/gene_level/fitting_results.tsv`. **Phase 2 develops/tests against `LD_DIT_HAP_generationRAW`** (closest analog). The `{dataset}` wildcard keeps HD switchable once its release/ is packaged. `default_dataset` in `datasets.yaml` stays `HD_DIT_HAP_generationRAW` (semantic default), but rule invocations in tests target LD.

---

## Task 1: Clustering — deterministic candidate labeling

**Source:** `gene_level_clustering.ipynb` cells 4, 6, 11, 14, 16, 18 (deterministic only).

**Files:**
- Create: `workflow/scripts/clustering/generate_candidate_clusters.py`
- Create: `workflow/rules/clustering.smk`
- Create: `tests/test_clustering.py`

**Inputs:** `{release}/gene_level/fitting_results.tsv` (via `data_config.load_dataset_config(dataset).gene_level.fitting_results`), `resources/curated/essentiality_verification.csv`.

**Outputs:** `results/clustering/candidates/{dataset}/candidate_clusters.tsv` (all fit columns + `RevisedDeletion_essentiality` + `cluster`), `results/clustering/candidates/{dataset}/clustering_metrics.tsv` (k-sweep + per-method silhouette/CH/DB).

**Steps:**
1. Load fitting_results (`index_col=0`), inject `RevisedDeletion_essentiality` at position 3 (verified dict, fallback to `DeletionLibrary_essentiality`) — quirk #2.
2. Custom scaling (quirk #1): `um` capped at 1.3 (`x if x < 1.3 else 1.3`), `lam`/10; `dropna()` defines the clustered set.
3. `evaluate_cluster_numbers(k_range=range(2,21))` → metrics table.
4. 4 algorithms at `n_clusters=64`, `random_state=42`, KMeans `n_init=10`: kmeans, agg-ward-euclidean, div-complete-cityblock (`fcluster-1`), gmm.
5. **Pin `best_method="kmeans"` explicitly** (quirk #3 — drop the fragile `set()[0]` selection); map labels back to `data_df`.

**Quirks to preserve:** custom scaling exactly; 0-based labels across all methods; `random_state=42` + `n_init=10`. Drop all viz cells (8-9,12,20,22,24,27,29,36-38,41) and cell-40 pairwise-distance (100MB exploratory).

**Verify:** unit test on synthetic 2-D data (scaling cap/divide, cluster count); integration run on `LD_DIT_HAP_generationRAW`; `snakemake -n` for the target.

**Commit:** `feat: port deterministic gene-level clustering (Phase 2 Task 1)`


## Task 2: Shared plotting modules + finalize_gene_clusters notebook

**Source:** `subset_visualization.py` (plot helpers) + `gene_level_clustering.ipynb` cell-25 (manual merge).

**Files:**
- Create: `workflow/src/plotting/style.py` (mplstyle load, AX_WIDTH/AX_HEIGHT/COLORS — design doc §7)
- Create: `workflow/src/plotting/gene_level.py` (plot_depletion_curves_for_groups/given_genes, plot_groups_on_feature_space, plot_given_genes_on_feature_space — from subset_visualization.py; hardcodes YES0..YES4 columns)
- Create: `notebooks/clustering/finalize_gene_clusters.ipynb`
- Create: `config/DIT_HAP.mplstyle` (if not already present from Task 1 of Phase 1)

**finalize notebook contract (first markdown cell):**
```
## Inputs
- results/clustering/candidates/{dataset}/candidate_clusters.tsv  (from: workflow/scripts/clustering/generate_candidate_clusters.py)
## Outputs
- resources/curated/final_clusters.tsv  (manually curated — the 64->9 merge, consumed by enrichment.smk + ml.smk)
```

**Steps:**
1. Port `plotting/style.py` and `plotting/gene_level.py`; `generic.py` (scatter/donut, no biology) is optional this phase — port only what the notebook needs.
2. Notebook: read candidates → viz feature-space scatter (cell-24 equivalent) → apply the two hand-built dicts (`reformat_cluster` 64→9, `reorder_reformat_cluster` renumber to 1..9) → write `final_clusters.tsv` with `revised_cluster` column + final viz PDFs.
3. **This is a human step** — the plan delivers the notebook scaffold with the current dicts pre-filled from cell-25 as a starting point; the user re-runs and adjusts against LD data.

**Quirks:** `cluster_minus_one=True` in final plot (revised_cluster is 1-based); 10-color palette assumes ≤10 clusters; WT cluster ends up = 9.

**Verify:** notebook runs top-to-bottom against LD candidates; `final_clusters.tsv` has expected columns (`Systematic ID`, `revised_cluster`, `um`, `lam`, `A`, YES fitted cols). No pytest (human notebook).

**Commit:** `feat: add plotting modules + finalize_gene_clusters notebook (Phase 2 Task 2)`

## Task 3: Enrichment pipeline modules (workflow/src/enrichment/)

**Source:** `enrichment_functions.py` (functions beyond what ontology.py already has).

**Files:**
- Modify: `workflow/src/enrichment/ontology.py` (add `ns2slim_assoc` via `mapslim`/`get_slim_ns2assoc` — the 6-tuple→7-tuple gap; add `GeneMetaData`/`GeneMetaConfig`; add `format_phaf_file`/`format_mondo_gaf_file`)
- Create: `workflow/src/enrichment/pipeline.py` (`ontology_enrichment_pipeline`, `ontology_enrichment`, `format_ontology_enrichment_results`, `assign_term_name`, `create_enrichment_dataframe`)
- Create: `tests/test_enrichment_pipeline.py`

**Steps:**
1. Close the slim-assoc gap: extend `load_ontology_data` (or pipeline) to compute `ns2slim_assoc` so the slim enrichment path works.
2. Port `GeneMetaData`/`GeneMetaConfig` (reads deletion library xlsx → now `resources/curated/deletion_library_categories.xlsx`).
3. Port `format_phaf_file`/`format_mondo_gaf_file` — **redirect their output** to a proper results/intermediate dir (not into `resources/` input tree), and **drop the `date.today()` header stamp** (quirk — breaks reproducibility). Preserve PHAF filter: `Allele type in {deletion,disruption}` AND `Condition contains FYECO:0000005`.
4. Port `ontology_enrichment_pipeline` (goatools `GOEnrichmentStudyNS`, `alpha=0.05`, `methods=["fdr_bh"]`, filter `p_fdr_bh < alpha and enrichment=="e"`, GO `relationships={is_a,part_of}` + `propagate_counts=True`).

**Verify:** unit tests — `ns2slim_assoc` non-empty for GO; `format_phaf_file` output has no date stamp + correct filter; `GeneMetaData` loads. Run against real PomBase 2025-10-01 OBO/GAF.

**Commit:** `feat: port enrichment pipeline modules with slim-assoc fix (Phase 2 Task 3)`

## Task 4: Enrichment — deterministic goatools rule

**Source:** `comprehensive_enrichment_analysis.ipynb` (goatools GO/FYPO/MONDO paths only; exclude STRING/REVIGO/section-9/TEST).

**Files:**
- Create: `workflow/scripts/enrichment/run_ontology_enrichment.py`
- Create: `workflow/rules/enrichment.smk`
- Create: `tests/test_run_ontology_enrichment.py`

**Inputs:** `resources/curated/final_clusters.tsv`, PomBase ontology triples (GO/FYPO/MONDO obo+gaf+slim under `resources/external/pombase/{version}/ontologies_and_associations/`), `resources/curated/deletion_library_categories.xlsx`.

**Outputs (per dataset):**
- `results/enrichment/raw/{dataset}/go_enrichment_full_filtered.tsv` (**design-doc target name**; content = `pop_count<400`, `namespace!='MF'` filtered GO enrichment — the notebook's `_noMF` output renamed to match the doc)
- `results/enrichment/raw/{dataset}/gene_ontology_enrichment_results.xlsx` + FYPO + MONDO workbooks (per-namespace sheets)
- Gene list txt files (`DIT_HAP_all_genes.txt`, per-cluster, `all_coding_genes.txt`, `not_covered_genes.txt`)

**Steps:**
1. Read `final_clusters.tsv`, split by `revised_cluster` (values 1..9). Emit gene lists.
2. Per cluster: GO/FYPO/MONDO enrichment via `ontology_enrichment_pipeline` (bg = all genes). For clusters 1-8 additionally run **non-WT comparison** (bg = genes with `revised_cluster<=8`) — quirk: WT cluster=9, `cluster<=8` branch. Make WT threshold a config param.
3. Produce full + slim results per ontology; write filtered GO tsv + xlsx workbooks.

**Quirks:** WT=9 hardcoded→config; `pop_count<400` + `namespace!='MF'` post-filter; `format_ontology_enrichment_results` try/except fallback (get_goea_nts_prt else create_enrichment_dataframe).

**Verify:** unit test cluster-split + gene-list emission on synthetic clusters; integration on LD `final_clusters.tsv`; confirm `go_enrichment_full_filtered.tsv` non-empty with expected columns; `snakemake -n`.

**Commit:** `feat: port deterministic goatools enrichment rule (Phase 2 Task 4)`

## Task 5: Enrichment — network rule (STRING + REVIGO, optional)

**Source:** `comprehensive_enrichment_analysis.ipynb` STRING/REVIGO paths.

**Files:**
- Create: `workflow/scripts/enrichment/run_network_enrichment.py`
- Create: `workflow/rules/enrichment_network.smk`
- Modify: `workflow/src/enrichment/pipeline.py` (add `stringdb_enrichment`, `revigo_analysis`, `format_string_enrichment_results`, response caching)

**Inputs:** gene lists from Task 4 output, cached API responses under `resources/external/enrichment_cache/{dataset}/`.

**Outputs:** `results/enrichment/network/{dataset}/string_enrichment_results.xlsx`, REVIGO-annotated GO tables (`Representative_*`, `Eliminated_*`, `Dispensability_*` columns).

**Steps:**
1. Port STRING (`stringdb_enrichment`, species 4896, `get_string_ids`+`enrichment`, MAX_RETRIES=5) and REVIGO (`revigo_analysis`, POST revigo.irb.hr, species 284812, cutoffs [0.7,0.5], measure SIMREL).
2. **Cache-first**: check `resources/external/enrichment_cache/` before hitting API; write responses there. This makes re-runs deterministic given cached responses.
3. Rule is **not in `rule all`** by default — invoked explicitly (network dependency). Wrap STRING in try/except (quirk: main STRING call currently unwrapped can hard-fail).

**Quirks:** REVIGO HTML-scrape (`pd.read_html`), GO-id zero-pad reconstruction; per namespace×cutoff×cluster×(full+nonWT) → dozens of calls.

**Verify:** unit test cache hit/miss logic (mock responses, no live network in tests); mark live-API test as `@pytest.mark.network` skipped by default. Document that first run needs network.

**Commit:** `feat: add optional STRING/REVIGO network enrichment rule with caching (Phase 2 Task 5)`

## Task 6: ML — feature/target preparation

**Source:** `machine_learning_data_preparation.ipynb` (DIT_HAP_pipeline canonical version, NOT Feature_organization).

**Files:**
- Create: `workflow/scripts/ml/prepare_features_targets.py`
- Create: `tests/test_prepare_features_targets.py`

**Inputs:** `results/features/{pombase_version}/pombe_coding_gene_protein_features.tsv` (from Phase 1), `resources/curated/final_clusters.tsv`.

**Outputs (per dataset):** `results/ml/features_targets/{dataset}/{split}_transformed_{features,targets,features_and_targets}.csv` for splits `all|um_gt_p35|um_le_p35|nonWT`, plus `all_features_with_target_values.csv`, `missing_value_analysis.csv`.

**Steps:**
1. Left-merge feature matrix (`gene_systematic_id`) with targets from `final_clusters.tsv` (`Systematic ID`→`Systematic_ID`, `revised_cluster`→`DIT_HAP_cluster`, cols `A,um,lam`).
2. Targets: `A,um,lam` (regression) + `DIT_HAP_cluster` (from `revised_cluster` — NOT `cluster`).
3. Splits: `all`, `um>0.35` (`um_gt_p35`), `um<=0.35` (`um_le_p35`), `nonWT` (`DIT_HAP_cluster!=9`). Thresholds→config.
4. Per-feature transform via `selected_features_and_transformations` dict (~83 entries: StandardScaler / PowerTransformer / OneHotEncoder(`dummy_na=True`) / binary-encode `+`→1/`-`→0 / set_index) — apply column-by-column with `.reshape(-1,1)` to match byte-for-byte (quirk #5).
5. **No imputation, no train/test split** — `dropna(how='any')` only (drop the unused SimpleImputer/train_test_split imports). Float format `%.5f` raw, `%.3f` transformed.

**Quirks:** WT=9 (not 1); `revised_cluster` source; scaler fit on full split (leakage — preserve); handle our matrix's `protein_half_life_minutes` + duplicate `DeletionLibrary_essentiality` (drop viz that references `t1/2 (min)`); `dummy_na=True` degenerate `Chromosome_nan` col.

**Verify:** unit test transform dict application (each transform type) on synthetic frame; integration on LD; confirm split row counts + column set; `snakemake -n`.

**Commit:** `feat: port ML feature/target preparation (Phase 2 Task 6)`

## Task 7: ML — AutoML analysis (mljar, Explain + Perform)

**Source:** `machine_learning_analysis.ipynb`.

**Files:**
- Create: `workflow/scripts/ml/train_automl.py`
- Create: `workflow/rules/ml.smk`
- Create: `workflow/envs/machine_learning.yml` (pin `mljar-supervised` version + scikit-learn, xgboost, lightgbm, catboost, shap)
- Create: `tests/test_train_automl.py`

**Inputs:** `results/ml/features_targets/{dataset}/all_transformed_features_and_targets.csv` (from Task 6; replaces the notebook's in-Config merge+`um>0.3` filter).

**Outputs (per dataset × target × mode):** `results/ml/models/{dataset}/{target}_{mode}/` (mljar tree: leaderboard.csv, per-model dirs, predictions), `prediction_and_residuals.pdf`, and a **new** persisted `metrics.tsv` (R2/RMSE/MAE/Pearson) + `features_importance.csv` (aggregated from per-model `learner_fold_*_importance.csv` — the notebook's dead-code bug, quirk #1).

**Steps:**
1. Fan out target∈{um,lam} × mode∈{Explain,Perform}. `AutoML(mode, ml_task="regression", results_path, explain_level=2)` — **explicitly pass `random_state`** and **set `total_time_limit` large** (Perform must finish full algorithm list; quirk — default 3600s can silently skip → non-deterministic).
2. train_test_split (`test_size=0.2`, `random_state=42`); PowerTransform features+target on train only; predict; inverse-transform for original-scale metrics.
3. **Persist fitted scalers** (joblib — needed to inverse-transform future predictions; quirk #3). **Clean output dir before run** (mljar errors/resumes on non-empty; quirk #5). Aggregate feature importance manually → write `features_importance.csv` + top-20 PDF.
4. Replace deprecated `mean_squared_error(squared=False)` → `root_mean_squared_error`.

**Quirks:** `um>0.3` filter (now in Task 6 splits or config); metrics on original scale vs mljar leaderboard on transformed scale (document); `A`/`DIT_HAP_cluster` declared-but-unmodeled.

**Verify:** unit test metrics computation + importance aggregation on synthetic predictions (no mljar training in unit tests — too slow); integration = one Explain run on LD (~60s) checked into a `@pytest.mark.slow` or manual step; `snakemake -n`. Perform mode validated manually (30-60min).

**Commit:** `feat: port mljar AutoML analysis rule (Phase 2 Task 7)`

## Task 8: Integration — Snakefile wiring, config, docs

**Files:**
- Modify: `Snakefile` (activate clustering/enrichment/ml includes; extend `rule all` with commented per-stage targets per repo convention)
- Modify: `config/analysis.yaml` (populate: clustering `n_clusters=64`/`random_state=42`/`k_range`/scaling; enrichment `fdr_threshold=0.05`/`wt_cluster=9`/`pop_count_max=400`; ml `um_split=0.35`/`test_size=0.2`/transform dict/mljar `total_time_limit`+`random_state`)
- Modify: `README.md` (document the clustering→finalize→enrichment/ml chain + the manual notebook step)
- Create: `workflow/src/data_config.py` accessors if needed (e.g. `final_clusters_path()`, `clustering_candidates_path(dataset)`)

**Steps:**
1. Wire includes; add `wildcard_constraints` for `target`/`mode` if needed. Keep `rule all` targets commented (run-on-demand convention).
2. Full config population from the PARAMETERS sections of all 4 analyses.
3. `snakemake -n` for each stage target; full `pytest tests/ -v`.
4. Full wet-run of the deterministic chain on LD: clustering → (manual finalize, or a documented note that it needs human input) → enrichment → ml (Explain). Verify with `--use-conda` in fresh envs.

**Verify:** all unit tests pass; `snakemake -n` clean for clustering/enrichment/ml targets; deterministic wet-run reaches enrichment + ml-Explain outputs (ml-Perform + network enrichment validated separately/manually).

**Commit:** `feat: wire Phase 2 core-chain rules + config + docs (Phase 2 Task 8)`

---

## Dependency order (execution)

```
Task 1 (clustering script) ──┐
Task 2 (plotting + finalize) ─┴─> [HUMAN: run finalize_gene_clusters.ipynb] ──> resources/curated/final_clusters.tsv
                                                                                        │
Task 3 (enrichment modules) ──> Task 4 (goatools rule) <─────────────────────────────┤
                                Task 5 (network rule, optional) <──── Task 4 outputs   │
                                Task 6 (ml prep) <────────────────────────────────────┤
                                Task 7 (ml automl) <──── Task 6 outputs                │
Task 8 (integration) <──── all above
```

Tasks 1-3 have no cross-dependency and could be parallelized. Task 4 needs Task 3 + the human finalize step. Tasks 6-7 need the same human step. **The human finalize step (Task 2 notebook) is the critical serialization point** — enrichment and ml cannot be wet-run end-to-end until `final_clusters.tsv` exists.
