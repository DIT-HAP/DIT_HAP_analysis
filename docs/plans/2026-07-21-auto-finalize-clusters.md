# Auto-Finalize Clusters Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a deterministic, no-human-judgment Snakemake path that produces `final_clusters.tsv` (cluster directly to k=9 with DR-based WT numbering), so the core chain no longer depends only on the manual finalize notebook.

**Architecture:** Reuse the existing `prepare_clustering_data` spine (its `annotated_data.pkl` + `scaled_data.pkl`), then a new `auto_finalize` step runs kmeans at k=9 and deterministically renumbers clusters (lowest mean DR = WT = 9, others 1..8 by ascending mean DR). Downstream enrichment/ml pick the auto vs manual file via a `config.clustering.finalize_mode` switch. The final cluster column is unified to `cluster` for both paths; the manual path keeps its pre-merge labels in a new `raw_cluster` column.

**Tech Stack:** Python 3.12, pandas, scikit-learn (KMeans), Snakemake 9, pytest, loguru.

**Design doc:** `docs/plans/2026-07-21-auto-finalize-clusters-design.md`

**Reference skills:** @superpowers:test-driven-development, @superpowers:verification-before-completion

---

## Notes for the implementer

- All `pytest` runs from the **worktree repo root**; tests add repo root to `sys.path` themselves.
- Follow the existing driver-script shape (see `workflow/scripts/clustering/prepare_clustering_data.py`): frozen `@dataclass` config with `validate()`, `setup_logger()`, `@logger.catch(reraise=True)` on `run()`, argparse `main()` returning int, `sys.path.insert(0, ...parents[3])` before `from workflow.src...`.
- Tasks 1-2 are pure-Python + testable now. Tasks 3-6 are the "final contract rename" (`revised_cluster` -> `cluster`) — they must land together in one commit each so the suite stays green. Task 7-9 are Snakemake wiring + notebook, verified by `snakemake -n`.
- The manual notebook and the 64-candidate pipeline are NOT removed.

---

### Task 1: `auto_finalize()` core function

**Files:**
- Modify: `workflow/src/clustering/candidates.py` (add function + import `pandas`/`numpy` already present)
- Test: `tests/test_clustering.py` (append tests)

**Step 1: Write the failing tests**

Append to `tests/test_clustering.py` (the file already imports `numpy as np`, `pandas as pd`, `pytest`, and from `workflow.src.clustering.candidates`):

```python
from workflow.src.clustering.candidates import auto_finalize


def _toy_annotated_scaled(n_per=20, seed=0):
    """Build a small annotated table + matching scaled matrix with 9 well-separated blobs in (DR, DL)."""
    rng = np.random.default_rng(seed)
    # 9 centers with strictly increasing DR so the ordering is unambiguous.
    centers = [(0.05, 0.1), (0.2, 0.3), (0.35, 0.2), (0.5, 0.5), (0.65, 0.4),
               (0.8, 0.6), (0.95, 0.3), (1.1, 0.7), (1.25, 0.5)]
    drs, dls, idx = [], [], []
    for c, (dr, dl) in enumerate(centers):
        for i in range(n_per):
            drs.append(dr + rng.normal(0, 0.005))
            dls.append(dl + rng.normal(0, 0.005))
            idx.append(f"g{c}_{i}")
    annotated = pd.DataFrame({"DR": drs, "DL": dls, "A": 1.0}, index=idx)
    annotated.index.name = "Systematic ID"
    scaled = annotated[["DR", "DL"]].copy()
    return annotated, scaled


def test_auto_finalize_produces_k_clusters_labelled_1_to_9():
    annotated, scaled = _toy_annotated_scaled()
    out = auto_finalize(annotated, scaled, n_clusters=9, random_state=42, wt_cluster=9)
    assert "cluster" in out.columns
    assert sorted(out["cluster"].unique()) == list(range(1, 10))
    # Full annotated columns are preserved; index name intact.
    assert out.index.name == "Systematic ID"
    assert {"DR", "DL", "A"}.issubset(out.columns)


def test_auto_finalize_assigns_lowest_DR_to_wt_cluster():
    annotated, scaled = _toy_annotated_scaled()
    out = auto_finalize(annotated, scaled, n_clusters=9, random_state=42, wt_cluster=9)
    means = out.groupby("cluster")["DR"].mean()
    # WT (id 9) has the lowest mean DR; ascending DR => monotonic 1..9 mapping.
    assert means.idxmin() == 9
    assert list(means.sort_index().index) == list(range(1, 10))
    assert means.is_monotonic_increasing


def test_auto_finalize_is_deterministic():
    annotated, scaled = _toy_annotated_scaled()
    a = auto_finalize(annotated, scaled, n_clusters=9, random_state=42, wt_cluster=9)
    b = auto_finalize(annotated, scaled, n_clusters=9, random_state=42, wt_cluster=9)
    pd.testing.assert_series_equal(a["cluster"], b["cluster"])


def test_auto_finalize_only_labels_scaled_genes():
    """Genes dropped by scaling (NaN DR/DL) get no cluster (NaN), matching the clustered set."""
    annotated, scaled = _toy_annotated_scaled(n_per=15)
    extra = pd.DataFrame({"DR": [np.nan], "DL": [np.nan], "A": [1.0]}, index=["ghost"])
    extra.index.name = "Systematic ID"
    annotated2 = pd.concat([annotated, extra])
    out = auto_finalize(annotated2, scaled, n_clusters=9, random_state=42, wt_cluster=9)
    assert pd.isna(out.loc["ghost", "cluster"])
    assert out.loc["g0_0", "cluster"] in range(1, 10)
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_clustering.py -k auto_finalize -v`
Expected: FAIL — `ImportError: cannot import name 'auto_finalize'`.

**Step 3: Implement `auto_finalize`**

Add to the end of `workflow/src/clustering/candidates.py`:

```python
# =============================================================================
# AUTOMATIC FINALIZE (deterministic k=9, no human merge)
# =============================================================================
# Number of final clusters for the automatic finalize path (design doc §3).
FINAL_N_CLUSTERS = 9


@logger.catch(reraise=True)
def auto_finalize(
    annotated: pd.DataFrame,
    scaled: pd.DataFrame,
    n_clusters: int = FINAL_N_CLUSTERS,
    random_state: int = 42,
    wt_cluster: int = 9,
) -> pd.DataFrame:
    """Cluster the scaled (DR, DL) matrix to n_clusters via kmeans and deterministically
    renumber to 1..n_clusters: lowest mean DR = WT (assigned wt_cluster), the rest in
    ascending mean-DR order. Returns the annotated table with a final `cluster` column
    (NaN for genes not in the scaled/clustered set). See design doc §3-4.
    """
    raw = pd.Series(
        cluster_one_method(BEST_METHOD, scaled, n_clusters, random_state),
        index=scaled.index,
        name="_raw",
    )
    # Rank raw clusters by mean DR (ascending), tie-broken by mean DL then raw id,
    # so the numbering is fully reproducible across runs.
    stats = (
        annotated.loc[scaled.index, ["DR", "DL"]]
        .assign(_raw=raw)
        .groupby("_raw")
        .agg(mean_dr=("DR", "mean"), mean_dl=("DL", "mean"))
        .reset_index()
        .sort_values(["mean_dr", "mean_dl", "_raw"], kind="stable")
        .reset_index(drop=True)
    )
    # Ascending DR -> final ids 1..n_clusters; the lowest-DR row becomes wt_cluster.
    final_ids = list(range(1, n_clusters + 1))
    # Place WT id at rank 0 (lowest DR); remaining ranks fill the other ids in order.
    remaining = [i for i in final_ids if i != wt_cluster]
    ordered_ids = [wt_cluster] + remaining
    raw_to_final = {row._raw: ordered_ids[rank] for rank, row in enumerate(stats.itertuples(index=False))}

    out = annotated.copy()
    out["cluster"] = raw.map(raw_to_final)
    logger.info(f"Auto-finalized {raw.notna().sum()} genes into {n_clusters} clusters (WT={wt_cluster})")
    return out
```

**IMPORTANT ordering note:** the tests assert `means.is_monotonic_increasing` AND `idxmin() == 9`. That requires the mapping to be: rank 0 (lowest DR) -> id 9, rank 1 -> id 1, rank 2 -> id 2, ..., rank 8 -> id 8. So mean DR by final id is: id1<id2<...<id8, and id9 is the lowest of all. `means.sort_index()` = [id1..id9]; is it monotonic increasing? id1..id8 increasing, but id9 (lowest) sits last and would break monotonicity. **Resolve before implementing:** change the test to assert only `means.idxmin() == 9` and that ids 1..8 are monotonic in DR (drop the all-9 monotonic assertion), OR order so WT is the max id but still lowest DR (impossible to be monotonic over all 9). Keep the design's rule (lowest DR = 9) and adjust the test in Step 1 to:

```python
    means = out.groupby("cluster")["DR"].mean()
    assert means.idxmin() == 9                       # WT = lowest DR
    non_wt = means.drop(index=9).sort_index()
    assert non_wt.is_monotonic_increasing            # ids 1..8 ascend in DR
```

Apply that corrected assertion in Step 1 before running.

**Step 4: Run tests to verify they pass**

Run: `pytest tests/test_clustering.py -k auto_finalize -v`
Expected: PASS (4 tests).

**Step 5: Commit**

```bash
git add workflow/src/clustering/candidates.py tests/test_clustering.py
git commit -m "feat(clustering): add deterministic auto_finalize (k=9, DR-based WT numbering)"
```

---

### Task 2: `auto_finalize_clusters.py` driver script

**Files:**
- Create: `workflow/scripts/clustering/auto_finalize_clusters.py`
- Test: `tests/test_clustering.py` (append a config + end-to-end-on-pickles test)

**Step 1: Write the failing test**

Append to `tests/test_clustering.py`:

```python
from workflow.scripts.clustering.auto_finalize_clusters import AutoFinalizeConfig, run as run_auto_finalize


def test_auto_finalize_driver_writes_final_tsv(tmp_path):
    annotated, scaled = _toy_annotated_scaled()
    ap = tmp_path / "annotated.pkl"
    sp = tmp_path / "scaled.pkl"
    annotated.to_pickle(ap)
    scaled.to_pickle(sp)
    out = tmp_path / "final_clusters.tsv"

    cfg = AutoFinalizeConfig(
        annotated_data=ap, scaled_data=sp, output=out,
        n_clusters=9, random_state=42, wt_cluster=9,
    )
    run_auto_finalize(cfg)

    df = pd.read_csv(out, sep="\t")
    assert "Systematic ID" in df.columns          # index written with its name
    assert "cluster" in df.columns
    assert "raw_cluster" not in df.columns          # auto path has no pre-merge labels
    assert sorted(df["cluster"].dropna().unique()) == list(range(1, 10))
```

**Step 2: Run to verify it fails**

Run: `pytest tests/test_clustering.py -k auto_finalize_driver -v`
Expected: FAIL — `ModuleNotFoundError`.

**Step 3: Implement the driver**

Create `workflow/scripts/clustering/auto_finalize_clusters.py` (mirror `prepare_clustering_data.py`'s structure exactly):

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Automatic Cluster Finalization (deterministic, no human merge)
================================================================

The auto alternative to notebooks/clustering/finalize_gene_clusters.ipynb: reads
the prepare_clustering_data spine pickles, clusters the scaled (DR, DL) matrix to
k=9 with kmeans, and deterministically renumbers clusters (lowest mean DR = WT).
Writes final_clusters.tsv with the unified `cluster` column consumed by
enrichment.smk + ml.smk (design doc §2-4).

Input
-----
- annotated_data.pkl (from prepare_clustering_data)
- scaled_data.pkl (from prepare_clustering_data)

Output
------
- final_clusters.tsv: full annotated table + final `cluster` (1..9, WT=9);
  index = systematic ID. No raw_cluster (auto path has no pre-merge labels).

Usage
-----
    python auto_finalize_clusters.py \\
        --annotated-data results/clustering/candidates/{dataset}/_work/annotated_data.pkl \\
        --scaled-data    results/clustering/candidates/{dataset}/_work/scaled_data.pkl \\
        --output         results/clustering/final/{dataset}/final_clusters.tsv \\
        --n-clusters 9 --random-state 42 --wt-cluster 9

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-21
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.clustering.candidates import FINAL_N_CLUSTERS, auto_finalize


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class AutoFinalizeConfig:
    """Inputs, output, and clustering params for the automatic finalize path."""
    annotated_data: Path
    scaled_data: Path
    output: Path
    n_clusters: int = FINAL_N_CLUSTERS
    random_state: int = 42
    wt_cluster: int = 9

    def validate(self) -> None:
        """Raise ValueError if any required input is missing, then ensure output dir exists."""
        for path in [self.annotated_data, self.scaled_data]:
            if not path.exists():
                raise ValueError(f"Required input not found: {path}")
        self.output.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(sys.stdout, format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}", level=log_level, colorize=False)


# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch(reraise=True)
def run(config: AutoFinalizeConfig) -> None:
    """Cluster to k=9, renumber deterministically, write final_clusters.tsv."""
    config.validate()
    annotated = pd.read_pickle(config.annotated_data)
    scaled = pd.read_pickle(config.scaled_data)
    out = auto_finalize(
        annotated, scaled,
        n_clusters=config.n_clusters,
        random_state=config.random_state,
        wt_cluster=config.wt_cluster,
    )
    out.to_csv(config.output, sep="\t", index=True)
    logger.success(f"Wrote {len(out)} genes ({out['cluster'].notna().sum()} clustered) to {config.output}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Automatic deterministic cluster finalization (k=9)")
    parser.add_argument("--annotated-data", type=Path, required=True, help="Annotated data pickle (from prepare)")
    parser.add_argument("--scaled-data", type=Path, required=True, help="Scaled (DR, DL) matrix pickle (from prepare)")
    parser.add_argument("--output", type=Path, required=True, help="Output final_clusters.tsv")
    parser.add_argument("--n-clusters", type=int, default=FINAL_N_CLUSTERS, help=f"Final cluster count (default {FINAL_N_CLUSTERS})")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed (default 42)")
    parser.add_argument("--wt-cluster", type=int, default=9, help="WT/background cluster id (default 9)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: build config, run auto-finalize, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        config = AutoFinalizeConfig(
            annotated_data=args.annotated_data,
            scaled_data=args.scaled_data,
            output=args.output,
            n_clusters=args.n_clusters,
            random_state=args.random_state,
            wt_cluster=args.wt_cluster,
        )
        run(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1
    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
```

**Step 4: Run to verify it passes**

Run: `pytest tests/test_clustering.py -k auto_finalize -v`
Expected: PASS (all auto_finalize tests).

**Step 5: Commit**

```bash
git add workflow/scripts/clustering/auto_finalize_clusters.py tests/test_clustering.py
git commit -m "feat(clustering): add auto_finalize_clusters driver script"
```

---

### Task 3: Rename enrichment cluster column contract (`revised_cluster` -> `cluster`)

**Files:**
- Modify: `workflow/src/enrichment/cluster_enrichment.py:71`
- Test: `tests/test_cluster_enrichment.py:24-36`

**Step 1: Update the test to the new contract (write failing test)**

In `tests/test_cluster_enrichment.py`, change the helper + calls to use `cluster` instead of `revised_cluster`:
- The docstring/comment "Systematic ID + revised_cluster" -> "Systematic ID + cluster".
- The three dict rows `{"Systematic ID": ..., "revised_cluster": N}` -> `{"Systematic ID": ..., "cluster": N}`.
- `load_cluster_genesets(fc, cluster_column="revised_cluster", wt_cluster=9)` -> `cluster_column="cluster"` (or drop the kwarg if you make `CLUSTER_COLUMN` the default — see Step 3).

**Step 2: Run to verify it fails**

Run: `pytest tests/test_cluster_enrichment.py -v`
Expected: FAIL — the fixture now writes `cluster`, but production `CLUSTER_COLUMN` is still `revised_cluster`, so `groupby("revised_cluster")` raises `KeyError`.

**Step 3: Flip the production constant**

In `workflow/src/enrichment/cluster_enrichment.py`:
```python
CLUSTER_COLUMN = "cluster"   # final contract column (auto or manual); was revised_cluster
```

**Step 4: Run to verify it passes**

Run: `pytest tests/test_cluster_enrichment.py -v`
Expected: PASS.

**Step 5: Commit**

```bash
git add workflow/src/enrichment/cluster_enrichment.py tests/test_cluster_enrichment.py
git commit -m "refactor(enrichment): read final cluster column as `cluster` (was revised_cluster)"
```

---

### Task 4: Rename ml `load_modeling_data` cluster column

**Files:**
- Modify: `workflow/src/ml/data.py:56-58` and the docstring at `:15`
- Test: `tests/test_train_automl.py:75`

**Step 1: Update the test to the new contract**

In `tests/test_train_automl.py`, change the fixture DataFrame key `"revised_cluster"` -> `"cluster"` (line ~75).

**Step 2: Run to verify it fails**

Run: `pytest tests/test_train_automl.py -v`
Expected: FAIL — `load_modeling_data` renames `revised_cluster`->`DIT_HAP_cluster`, but the column is now `cluster`, so the `[["Systematic_ID","A","DR","DL","DIT_HAP_cluster"]]` selection raises `KeyError`.

**Step 3: Update production**

In `workflow/src/ml/data.py`, the rename dict:
```python
    targets = pd.read_csv(final_clusters, sep="\t").rename(
        columns={"Systematic ID": "Systematic_ID", "cluster": "DIT_HAP_cluster"}
    )[["Systematic_ID", "A", "DR", "DL", "DIT_HAP_cluster"]]
```
Also update the docstring line 15: `... A, DR, DL, revised_cluster)` -> `... A, DR, DL, cluster)`.

**Step 4: Run to verify it passes**

Run: `pytest tests/test_train_automl.py -v`
Expected: PASS.

**Step 5: Commit**

```bash
git add workflow/src/ml/data.py tests/test_train_automl.py
git commit -m "refactor(ml): read final cluster column as `cluster` in load_modeling_data"
```

---

### Task 5: Rename ml `prepare_features_targets` cluster column

**Files:**
- Modify: `workflow/scripts/ml/prepare_features_targets.py:188-190` (+ docstring `:21`)
- Test: `tests/test_prepare_features_targets.py:79`

**Step 1: Update the test to the new contract**

In `tests/test_prepare_features_targets.py` (line ~79) change fixture key `"revised_cluster"` -> `"cluster"`. Keep the test name `test_merge_uses_revised_cluster_as_dit_hap_cluster` or rename to `..._uses_cluster_...` for clarity (rename is optional; if renaming, do it here).

**Step 2: Run to verify it fails**

Run: `pytest tests/test_prepare_features_targets.py -v`
Expected: FAIL — `merge_features_targets` selects `revised_cluster` which no longer exists -> `KeyError`.

**Step 3: Update production**

In `workflow/scripts/ml/prepare_features_targets.py`:
```python
    # DIT_HAP_cluster comes from the final `cluster` column (auto or manual finalize).
    target_values = metrics[["Systematic_ID", "A", "DR", "DL", "cluster"]].rename(
        columns={"cluster": "DIT_HAP_cluster"}
    )
```
Update docstring line 21: `... A, DR, DL, revised_cluster)` -> `... A, DR, DL, cluster)`.

**Step 4: Run to verify it passes**

Run: `pytest tests/test_prepare_features_targets.py -v`
Expected: PASS.

**Step 5: Commit**

```bash
git add workflow/scripts/ml/prepare_features_targets.py tests/test_prepare_features_targets.py
git commit -m "refactor(ml): read final cluster column as `cluster` in prepare_features_targets"
```

---

### Task 6: Full suite green after the rename

**Step 1: Run the whole suite**

Run: `pytest -q`
Expected: PASS (no `revised_cluster` KeyErrors anywhere).

**Step 2: Grep for stragglers**

Run: `grep -rn "revised_cluster" workflow/ tests/ --include=*.py`
Expected: no matches in `workflow/` or `tests/` (the only remaining `revised_cluster` will be in the notebook, handled in Task 9, and in `prepare_ml_data.py`/`ml/data.py` docstrings — fix any doc stragglers found).

**Step 3: Commit any doc stragglers**

```bash
git add -A && git commit -m "docs(ml): scrub remaining revised_cluster references"
```
(Skip if grep is already clean.)

---

### Task 7: config switch + `auto_finalize_clusters` Snakemake rule

**Files:**
- Modify: `config/analysis.yaml` (clustering block)
- Modify: `workflow/rules/clustering.smk` (append rule)

**Step 1: Add config keys**

In `config/analysis.yaml`, under `clustering:` add:
```yaml
  finalize_mode: auto            # auto | manual — which final_clusters.tsv downstream reads
  finalize_mode_overrides: {}    # per-dataset override, e.g. {HD_DIT_HAP: manual}
  final_n_clusters: 9            # k for the automatic finalize path
```

**Step 2: Add the rule**

Append to `workflow/rules/clustering.smk`:
```python
# --- Automatic finalize (deterministic alternative to the manual notebook) ---
# Reuses the prepare spine pickles; clusters to k=9 and DR-numbers (design doc §2-3).
rule auto_finalize_clusters:
    input:
        annotated=f"{_CWORK}/annotated_data.pkl",
        scaled=f"{_CWORK}/scaled_data.pkl",
    output:
        clusters="results/clustering/final/{dataset}/final_clusters.tsv",
    params:
        n_clusters=config.get("clustering", {}).get("final_n_clusters", 9),
        random_state=config.get("clustering", {}).get("random_state", 42),
        wt_cluster=config.get("enrichment", {}).get("wt_cluster", 9),
    log:
        "logs/clustering/auto_finalize_clusters_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [clustering] Auto-finalizing clusters for {wildcards.dataset} (k={params.n_clusters})..."
    shell:
        """
        python workflow/scripts/clustering/auto_finalize_clusters.py \
            --annotated-data {input.annotated} \
            --scaled-data {input.scaled} \
            --output {output.clusters} \
            --n-clusters {params.n_clusters} \
            --random-state {params.random_state} \
            --wt-cluster {params.wt_cluster} &> {log}
        """
```

**Step 3: Verify the DAG builds**

Run: `snakemake -n results/clustering/final/HD_DIT_HAP/final_clusters.tsv`
Expected: dry-run lists `prepare_clustering_data` -> `auto_finalize_clusters` (2 jobs), no errors. (If the upstream `fitting_results.tsv` is absent on disk, the dry-run may report it as a missing input for `prepare_clustering_data` — that is fine; the rule graph itself must resolve without syntax/wildcard errors.)

**Step 4: Commit**

```bash
git add config/analysis.yaml workflow/rules/clustering.smk
git commit -m "feat(clustering): add auto_finalize_clusters rule + finalize_mode config"
```

---

### Task 8: downstream `finalize_mode` selection in enrichment.smk + ml.smk

**Files:**
- Modify: `workflow/rules/enrichment.smk` (helper + `prepare_genesets` input)
- Modify: `workflow/rules/ml.smk` (helper + `prepare_ml_data` input)

**Step 1: Add a shared selector helper**

At the top of `workflow/rules/clustering.smk` (so both includes can use it — it is included before enrichment/ml in the Snakefile), add:
```python
def final_clusters_path(dataset: str) -> str:
    """Return the final_clusters.tsv path per config.clustering.finalize_mode (+ per-dataset override)."""
    cl = config.get("clustering", {})
    mode = cl.get("finalize_mode_overrides", {}).get(dataset, cl.get("finalize_mode", "auto"))
    if mode == "manual":
        return "resources/curated/final_clusters.tsv"
    if mode == "auto":
        return f"results/clustering/final/{dataset}/final_clusters.tsv"
    raise ValueError(f"Unknown finalize_mode {mode!r} for dataset {dataset!r} (expected auto|manual)")
```

**Step 2: Wire enrichment**

In `workflow/rules/enrichment.smk`, `rule prepare_genesets`, replace:
```python
        final_clusters="resources/curated/final_clusters.tsv",
```
with:
```python
        final_clusters=lambda wc: final_clusters_path(wc.dataset),
```
Also update the module header comment block (lines ~16-19) that says "run notebooks/clustering/finalize_gene_clusters.ipynb first" to note the auto path is now default (manual only when `finalize_mode: manual`).

**Step 3: Wire ml**

In `workflow/rules/ml.smk`, `rule prepare_ml_data`, replace:
```python
        final_clusters="resources/curated/final_clusters.tsv",
```
with:
```python
        final_clusters=lambda wc: final_clusters_path(wc.dataset),
```

**Step 4: Verify both DAGs resolve in each mode**

Run (auto is default):
```
snakemake -n results/enrichment/raw/HD_DIT_HAP/2026-06-01/go_enrichment_full_filtered.tsv
snakemake -n results/ml/models/HD_DIT_HAP/2026-06-01/DR_Explain/metrics.tsv
```
Expected: both dry-runs route `final_clusters` through `results/clustering/final/HD_DIT_HAP/...` (i.e. `auto_finalize_clusters` appears in the job list).

Then confirm manual mode flips the source:
```
snakemake -n --config clustering='{"finalize_mode":"manual"}' results/ml/models/HD_DIT_HAP/2026-06-01/DR_Explain/metrics.tsv
```
Expected: `final_clusters` resolves to `resources/curated/final_clusters.tsv` and `auto_finalize_clusters` does NOT appear. (If `--config` nested-dict override is awkward, instead temporarily set `finalize_mode: manual` in `config/analysis.yaml`, run the dry-run, then revert.)

**Step 5: Commit**

```bash
git add workflow/rules/clustering.smk workflow/rules/enrichment.smk workflow/rules/ml.smk
git commit -m "feat(clustering): route enrichment/ml final_clusters via finalize_mode switch"
```

---

### Task 9: update the manual notebook to the new column contract

**Files:**
- Modify: `notebooks/clustering/finalize_gene_clusters.ipynb` (cells `2d207fbf`, `bdc825ba`, `496f2cb6`, `e6d2c6ba`, and the markdown headers)

**Step 1: Rename raw label column on load (cell `2d207fbf`)**

After reading candidates, rename the incoming `cluster` (the 64 raw candidate labels) to `raw_cluster`:
```python
data_df = pd.read_csv(CANDIDATES, sep="\t", index_col=0).rename(columns={"cluster": "raw_cluster"})
data_df.shape, sorted(data_df['raw_cluster'].unique())[:10]
```

**Step 2: Map from raw_cluster into final `cluster` (cell `bdc825ba`)**

```python
data_df['cluster'] = data_df['raw_cluster'].map(reformat_cluster)
data_df['cluster'] = data_df['cluster'].map(reorder_reformat_cluster)
data_df['cluster'].value_counts().sort_index()
```

**Step 3: Update review plot (cells `00d22296`, `496f2cb6`)**

- `00d22296`: `visualize_cluster_on_feature_space(data_df, 'raw_cluster')` (review the 64 candidates).
- `496f2cb6`: `visualize_cluster_on_feature_space(data_df, 'cluster', show_box=True, legend=True, cluster_minus_one=True)`.

**Step 4: Output note (cell `e6d2c6ba` / markdown)**

The `to_csv` cell is unchanged (writes the whole frame, now with both `cluster` and `raw_cluster`). Update the markdown header cell `513050a3` Outputs bullet to note the file now uses the unified `cluster` column plus `raw_cluster` (pre-merge labels), and that this is the *manual* alternative to the default auto path.

**Step 5: Verify the notebook is coherent (no execution needed if upstream data absent)**

Run: `jupyter nbconvert --to script --stdout notebooks/clustering/finalize_gene_clusters.ipynb | grep -n "cluster"`
Expected: no remaining `revised_cluster`; `raw_cluster` used for the 64 labels, `cluster` for the final 1..9.

**Step 6: Commit**

```bash
git add notebooks/clustering/finalize_gene_clusters.ipynb
git commit -m "refactor(clustering): manual notebook emits unified `cluster` + `raw_cluster`"
```

---

### Task 10: docs + CLAUDE.md sync

**Files:**
- Modify: `README.md` (core-chain section), `CLAUDE.md` (the finalize step description)

**Step 1: Update README core-chain**

In `README.md`, the "Core analysis chain" section: reflect that step 3 is now optional — the default `auto` path builds `results/clustering/final/{dataset}/final_clusters.tsv` with no manual step, and the `Missing input` note applies only under `finalize_mode: manual`.

**Step 2: Update CLAUDE.md**

In `CLAUDE.md`, the "core chain is intentionally NOT one DAG" block: note there are now two finalize paths (auto default / manual optional) selected by `config.clustering.finalize_mode`, and that the final contract column is `cluster` (manual also keeps `raw_cluster`).

**Step 3: Commit**

```bash
git add README.md CLAUDE.md
git commit -m "docs: document auto vs manual finalize paths + cluster column contract"
```
(Note: `CLAUDE.md` is git-ignored per `.gitignore`; the `git add` will no-op it. That is expected — edit it anyway for the working tree.)

---

## Final verification

- `pytest -q` — full suite green.
- `grep -rn "revised_cluster" workflow/ tests/ --include=*.py` — no matches.
- `snakemake -n results/clustering/final/HD_DIT_HAP/final_clusters.tsv` — auto path DAG resolves.
- `snakemake -n results/ml/models/HD_DIT_HAP/2026-06-01/DR_Explain/metrics.tsv` — routes through auto path by default.
- If upstream `fitting_results.tsv` is available on disk: actually build `results/clustering/final/HD_DIT_HAP/final_clusters.tsv` and confirm 9 clusters, `cluster` in 1..9, WT (id 9) has lowest mean DR.
