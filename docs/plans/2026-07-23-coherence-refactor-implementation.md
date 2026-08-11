# Coherence Refactor Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Split `complex.smk` into a clearer `coherence.smk` stage that fans out by annotation **source** (`go_macrocomplex` / `go_cc` / `go_bp`), decoupled from clustering, so new grouping databases are added by writing one adapter + one config line.

**Architecture:** Three deterministic stages per (dataset × source): `prepare_annotation` (source adapter → unified long-table) → `compute_coherence` (byte-faithful seeded permutation, reads long-table + fitting_results) → `plot_group_scatter` (config-namelist-driven per-group scatter). The byte-faithful coherence math in `workflow/src/coherence/metrics.py` is untouched. Design: [2026-07-23-coherence-refactor-design.md](2026-07-23-coherence-refactor-design.md).

**Tech Stack:** Python 3.12, pandas, numpy, scipy, goatools (GO GAF propagation, reused from `workflow/src/enrichment/ontology.py`), matplotlib, Snakemake 9, pytest.

---

## Conventions for the executor

- **Worktree:** `.worktrees/optimize-complex-smk` (branch `optimize-complex-smk`). Run all commands from there.
- **Test env:** `statistics_and_figure_plotting` (has scipy/sklearn/pandas/numpy/matplotlib; `pytest` already pip-installed). Activate with:
  `source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate statistics_and_figure_plotting`
- **goatools tests:** the `go_cc`/`go_bp` GAF loader needs `goatools`, which is NOT in the statistics env. Those loader tests build a tiny fake OBO+GAF fixture and are `pytest.importorskip("goatools")`-guarded so the suite still passes in the statistics env; verify them in the enrichment/`biopython.yml` env where goatools exists.
- **Baseline:** 24 tests pass in `tests/test_complex_coherence.py` (8) + `tests/test_domain_differences.py` (16). 8 unrelated collection errors (missing gffutils/Bio/requests/mljar) are pre-existing, NOT regressions — run coherence tests explicitly, don't judge on full-suite collection.
- **Do NOT touch** `workflow/src/coherence/metrics.py` (byte-faithful Weiszfeld + MPD permutation, seed 42).
- **Do NOT wire** `workflow/src/coherence/attribution.py` this round.
- Commit after every green step. Keep git history via `git mv` where noted.

---

## Task 1: Source-adapter module — unified long-table contract + macrocomplex loader

**Files:**
- Create: `workflow/src/coherence/sources.py`
- Test: `tests/test_coherence_sources.py`

The unified long-table contract (every loader returns exactly these columns, in order):
`source, group_id, group_name, Systematic ID, Name, n_group_genes`.
- `Name` missing → filled with the row's `Systematic ID` (never blank).
- `n_group_genes` = total annotated member count of that term (before any DR filter), per group.

`go_macrocomplex` reads PomBase `macromolecular_complex_annotation.tsv` (columns include
`complex_term_id`, `GO_term_name`, `systematic_id`, `symbol`). This is the flat-TSV path,
byte-faithful to the old `compute_complex_coherence.py` grouping (keyed on `GO_term_name`).

**Step 1: Write the failing test**

```python
# tests/test_coherence_sources.py
"""Contract tests for coherence source adapters (unified long-table schema)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.src.coherence.sources import (
    LONG_TABLE_COLUMNS,
    load_macrocomplex,
)

LONG_TABLE_COLUMNS_EXPECTED = [
    "source", "group_id", "group_name", "Systematic ID", "Name", "n_group_genes",
]


def _write_macrocomplex(tmp_path: Path) -> Path:
    """A tiny macromolecular_complex_annotation.tsv: 2 complexes, one member with no symbol."""
    df = pd.DataFrame(
        {
            "complex_term_id": ["GO:0001", "GO:0001", "GO:0002", "GO:0002"],
            "GO_term_name": ["alpha complex", "alpha complex", "beta complex", "beta complex"],
            "systematic_id": ["SPAC1", "SPAC2", "SPBC1", "SPBC2"],
            "symbol": ["gene1", None, "gene3", "gene4"],
        }
    )
    d = tmp_path / "ontologies_and_associations"
    d.mkdir(parents=True)
    path = tmp_path / "ontologies_and_associations" / "macromolecular_complex_annotation.tsv"
    df.to_csv(path, sep="\t", index=False)
    return tmp_path


def test_long_table_columns_constant():
    assert list(LONG_TABLE_COLUMNS) == LONG_TABLE_COLUMNS_EXPECTED


def test_macrocomplex_returns_contract_columns(tmp_path):
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_macrocomplex(pombase_dir)
    assert list(out.columns) == LONG_TABLE_COLUMNS_EXPECTED
    assert set(out["source"]) == {"go_macrocomplex"}


def test_macrocomplex_fills_missing_name_with_systematic_id(tmp_path):
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_macrocomplex(pombase_dir)
    row = out[out["Systematic ID"] == "SPAC2"].iloc[0]
    assert row["Name"] == "SPAC2"  # no symbol -> filled with systematic id


def test_macrocomplex_n_group_genes_is_per_term_total(tmp_path):
    pombase_dir = _write_macrocomplex(tmp_path)
    out = load_macrocomplex(pombase_dir)
    alpha = out[out["group_name"] == "alpha complex"]
    assert (alpha["n_group_genes"] == 2).all()
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_coherence_sources.py -v`
Expected: FAIL — `ImportError: cannot import name 'LONG_TABLE_COLUMNS'`.

**Step 3: Write minimal implementation**

```python
# workflow/src/coherence/sources.py
"""Source adapters mapping PomBase grouping databases to a unified coherence long-table.

Every adapter returns the same contract columns (LONG_TABLE_COLUMNS), so the
downstream compute/plot stages are source-agnostic. Add a database = add one
adapter here + one entry in SOURCE_LOADERS + one line in config.coherence.sources.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

LONG_TABLE_COLUMNS = [
    "source", "group_id", "group_name", "Systematic ID", "Name", "n_group_genes",
]

# PomBase macrocomplex annotation column -> canonical contract name.
_MACRO_RENAME = {
    "complex_term_id": "group_id",
    "GO_term_name": "group_name",
    "systematic_id": "Systematic ID",
    "symbol": "Name",
}


def _finalize(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Fill missing Name with Systematic ID, add source + n_group_genes, order columns."""
    df = df.copy()
    df["source"] = source
    df["Name"] = df["Name"].fillna(df["Systematic ID"])
    df = df.drop_duplicates(subset=["group_id", "Systematic ID"])
    counts = df.groupby("group_id")["Systematic ID"].transform("size")
    df["n_group_genes"] = counts
    return df[LONG_TABLE_COLUMNS].reset_index(drop=True)


def load_macrocomplex(pombase_dir: Path) -> pd.DataFrame:
    """Flat PomBase macromolecular_complex_annotation.tsv -> unified long-table."""
    path = Path(pombase_dir) / "ontologies_and_associations" / "macromolecular_complex_annotation.tsv"
    raw = pd.read_csv(path, sep="\t").rename(columns=_MACRO_RENAME)
    for required in ["group_id", "group_name", "Systematic ID"]:
        if required not in raw.columns:
            raise ValueError(f"macrocomplex annotation missing '{required}' (have: {list(raw.columns)})")
    if "Name" not in raw.columns:
        raw["Name"] = pd.NA
    return _finalize(raw[["group_id", "group_name", "Systematic ID", "Name"]], "go_macrocomplex")
```

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_coherence_sources.py -v`
Expected: PASS (3 tests + the constant test).

**Step 5: Commit**

```bash
git add workflow/src/coherence/sources.py tests/test_coherence_sources.py
git commit -m "feat(coherence): add source-adapter module + go_macrocomplex loader"
```

---

## Task 2: GAF-namespace loader (go_cc / go_bp) + SOURCE_LOADERS registry

**Files:**
- Modify: `workflow/src/coherence/sources.py`
- Test: `tests/test_coherence_sources.py`

`load_gaf_namespace(pombase_dir, namespace)` reuses the goatools loading path from
`workflow/src/enrichment/ontology.py` (`OntologyDataConfig(...).load_data()` +
`load_ontology_data(..., relationships={"is_a","part_of"}, propagate_counts=True)`),
then keeps only terms whose `dag[term].namespace` matches the requested namespace
(`cellular_component` for CC, `biological_process` for BP), and expands the propagated
`go2genes` term→geneset dict into the unified long-table. `group_name` = `dag[term].name`.

Because goatools is absent from the statistics env, these tests are guarded with
`pytest.importorskip("goatools")` and build a tiny OBO+GAF fixture on disk.

**Step 1: Write the failing test** (append to `tests/test_coherence_sources.py`)

```python
def _write_go_fixture(tmp_path: Path) -> Path:
    """A minimal go-basic.obo + gene_ontology_annotation.gaf.tsv with 1 CC + 1 BP term."""
    d = tmp_path / "ontologies_and_associations"
    d.mkdir(parents=True, exist_ok=True)
    obo = d / "go-basic.obo"
    obo.write_text(
        "format-version: 1.2\n\n"
        "[Term]\nid: GO:0000100\nname: test cc complex\nnamespace: cellular_component\n\n"
        "[Term]\nid: GO:0000200\nname: test bp process\nnamespace: biological_process\n\n"
    )
    gaf = d / "gene_ontology_annotation.gaf.tsv"
    # GAF 2.1: 17 tab-separated columns; col2=DB_Object_ID(gene), col5=GO_ID, col9=Aspect.
    lines = ["!gaf-version: 2.1"]
    def row(gene, go, aspect):
        cols = ["PomBase", gene, gene, "", go, "PMID:1", "IDA", "", aspect,
                "", "", "gene", "taxon:4896", "20250101", "PomBase", "", ""]
        return "\t".join(cols)
    lines += [row("SPAC1", "GO:0000100", "C"), row("SPAC2", "GO:0000100", "C"),
              row("SPBC1", "GO:0000200", "P"), row("SPBC2", "GO:0000200", "P")]
    gaf.write_text("\n".join(lines) + "\n")
    return tmp_path


def test_gaf_namespace_cc_only_returns_cc_terms(tmp_path):
    pytest.importorskip("goatools")
    from workflow.src.coherence.sources import load_gaf_namespace
    pombase_dir = _write_go_fixture(tmp_path)
    out = load_gaf_namespace(pombase_dir, "CC")
    assert list(out.columns) == LONG_TABLE_COLUMNS_EXPECTED
    assert set(out["group_id"]) == {"GO:0000100"}
    assert set(out["source"]) == {"go_cc"}


def test_gaf_namespace_bp_only_returns_bp_terms(tmp_path):
    pytest.importorskip("goatools")
    from workflow.src.coherence.sources import load_gaf_namespace
    pombase_dir = _write_go_fixture(tmp_path)
    out = load_gaf_namespace(pombase_dir, "BP")
    assert set(out["group_id"]) == {"GO:0000200"}
    assert set(out["source"]) == {"go_bp"}


def test_source_loaders_registry_has_three_sources():
    from workflow.src.coherence.sources import SOURCE_LOADERS
    assert set(SOURCE_LOADERS) == {"go_macrocomplex", "go_cc", "go_bp"}
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_coherence_sources.py -v`
Expected: the registry test FAILs (`SOURCE_LOADERS` missing); the two goatools tests either FAIL (import) or SKIP if goatools absent.

**Step 3: Write minimal implementation** (append to `workflow/src/coherence/sources.py`)

```python
# --- GO GAF namespace loader (go_cc / go_bp) -------------------------------
# goatools namespace string per short code; also the config source name.
_NS_LONG = {"CC": "cellular_component", "BP": "biological_process"}
_NS_SOURCE = {"CC": "go_cc", "BP": "go_bp"}
# Match enrichment.smk's GO propagation exactly.
_GO_LOAD_KWARGS = {"relationships": {"is_a", "part_of"}, "propagate_counts": True,
                   "load_obsolete": False, "prt": None}


def load_gaf_namespace(pombase_dir: Path, namespace: str) -> pd.DataFrame:
    """GO GAF for one namespace (CC/BP), goatools-propagated, -> unified long-table.

    Reuses enrichment/ontology.py's OBO+GAF loading (is_a/part_of propagation,
    propagate_counts=True), then keeps only terms in the requested namespace and
    expands the propagated go2genes dict.
    """
    from workflow.src.enrichment.ontology import OntologyDataConfig, load_ontology_data

    if namespace not in _NS_LONG:
        raise ValueError(f"namespace must be one of {sorted(_NS_LONG)}, got {namespace!r}")
    od = Path(pombase_dir) / "ontologies_and_associations"
    data = OntologyDataConfig(
        ontology_obo=od / "go-basic.obo",
        ontology_association_gaf=od / "gene_ontology_annotation.gaf.tsv",
        slim_terms_table=[],  # slim table not needed for raw term->gene expansion
    ).load_data()
    dag, _objanno, _ns2assoc, _gene2go, go2genes, _slim = load_ontology_data(data, **_GO_LOAD_KWARGS)

    ns_long = _NS_LONG[namespace]
    rows = []
    for term, genes in go2genes.items():
        rec = dag.get(term)
        if rec is None or rec.namespace != ns_long:
            continue
        for gene in genes:
            rows.append({"group_id": term, "group_name": rec.name,
                         "Systematic ID": gene, "Name": pd.NA})
    df = pd.DataFrame(rows, columns=["group_id", "group_name", "Systematic ID", "Name"])
    return _finalize(df, _NS_SOURCE[namespace])


SOURCE_LOADERS = {
    "go_macrocomplex": load_macrocomplex,
    "go_cc": lambda d: load_gaf_namespace(d, "CC"),
    "go_bp": lambda d: load_gaf_namespace(d, "BP"),
}
```

**Step 4: Run test to verify it passes**

Run (statistics env — goatools tests SKIP): `pytest tests/test_coherence_sources.py -v`
Expected: registry test PASS; goatools tests SKIPPED.
Run (biopython/enrichment env — goatools present): same command → all PASS.

**Step 5: Commit**

---

## Task 3: `prepare_annotation.py` script (stage [1])

**Files:**
- Create: `workflow/scripts/coherence/prepare_annotation.py`
- Test: `tests/test_coherence_prepare.py`

Thin CLI: `--source {go_macrocomplex|go_cc|go_bp}` + `--pombase-dir` + `--output`.
Looks up `SOURCE_LOADERS[source]`, runs it, writes the long-table TSV. Follows the
`python-script-conventions` skill (argparse + loguru + `@logger.catch(reraise=True)`,
`sys.path.insert` then `from workflow.src...`). Keep it thin — logic lives in `sources.py`.

The unit-testable core is `prepare(source, pombase_dir) -> pd.DataFrame` (calls the
registry). Test it against the macrocomplex fixture (no goatools needed).

**Step 1: Write the failing test**

```python
# tests/test_coherence_prepare.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import pandas as pd
import pytest


def _write_macrocomplex(tmp_path: Path) -> Path:
    df = pd.DataFrame({
        "complex_term_id": ["GO:0001", "GO:0001", "GO:0002"],
        "GO_term_name": ["alpha complex", "alpha complex", "beta complex"],
        "systematic_id": ["SPAC1", "SPAC2", "SPBC1"],
        "symbol": ["gene1", "gene2", "gene3"],
    })
    d = tmp_path / "ontologies_and_associations"
    d.mkdir(parents=True)
    df.to_csv(d / "macromolecular_complex_annotation.tsv", sep="\t", index=False)
    return tmp_path


def test_prepare_macrocomplex_returns_long_table(tmp_path):
    from prepare_annotation import prepare
    out = prepare("go_macrocomplex", _write_macrocomplex(tmp_path))
    assert list(out.columns) == ["source", "group_id", "group_name",
                                 "Systematic ID", "Name", "n_group_genes"]
    assert len(out) == 3


def test_prepare_unknown_source_raises(tmp_path):
    from prepare_annotation import prepare
    with pytest.raises(ValueError, match="unknown source"):
        prepare("kegg_nope", tmp_path)
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_coherence_prepare.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'prepare_annotation'`.

**Step 3: Write minimal implementation**

Create `workflow/scripts/coherence/prepare_annotation.py` following the conventions header. Core:

```python
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.coherence.sources import SOURCE_LOADERS  # noqa: E402


def prepare(source: str, pombase_dir: Path) -> pd.DataFrame:
    """Dispatch to the source adapter and return the unified long-table."""
    if source not in SOURCE_LOADERS:
        raise ValueError(f"unknown source {source!r} (have: {sorted(SOURCE_LOADERS)})")
    return SOURCE_LOADERS[source](Path(pombase_dir))


@logger.catch(reraise=True)
def run(source: str, pombase_dir: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    table = prepare(source, pombase_dir)
    table.to_csv(output, sep="\t", index=False)
    logger.success(f"[{source}] {len(table):,} rows, "
                   f"{table['group_id'].nunique():,} groups -> {output}")
```

Plus argparse (`--source`, `--pombase-dir`, `--output`, `-v`) + `main()` returning int + `setup_logger`.

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_coherence_prepare.py -v`
Expected: PASS (2 tests).

**Step 5: Commit**

```bash
git add workflow/scripts/coherence/prepare_annotation.py tests/test_coherence_prepare.py
git commit -m "feat(coherence): add prepare_annotation.py (source -> long-table)"
```

---

## Task 4: `git mv` + refactor `compute_complex_coherence.py` -> `compute_coherence.py`

**Files:**
- Rename: `workflow/scripts/complex/compute_complex_coherence.py` -> `workflow/scripts/coherence/compute_coherence.py` (via `git mv`)
- Modify: the renamed script
- Test: `tests/test_coherence_compute.py` (new, filter-logic unit tests)

**Step 0: Rename first (keep blame)**

```bash
git mv workflow/scripts/complex/compute_complex_coherence.py workflow/scripts/coherence/compute_coherence.py
git commit -m "refactor(coherence): git mv compute_complex_coherence.py -> coherence/compute_coherence.py"
```

**What changes in the refactor:**
- Inputs: replace `--final-clusters` with `--fitting-results` (release `gene_level/fitting_results.tsv`)
  and `--complex-annotation` with `--annotation` (the long-table from Task 3). Add `--source`.
  Add `--max-term-genes` (default 500). Rename `--min-size`/`--max-size` semantics unchanged
  (now `min_group_size`/`max_group_size`). Add optional `--features` (Task 6 panels; wire the flag
  now, default None).
- `load_final_clusters(...)` -> `load_fitting_results(path, dr_threshold)`: same body (legacy um/lam
  rename, inf→nan dropna, norm_DR/norm_DL, DR>threshold background) but reads `fitting_results.tsv`.
  Keep index handling: `fitting_results.tsv` has the systematic id as index col 0 (see
  `clustering/candidates.py` `read_file(..., index_col=[0])`) — reset it into a `Systematic ID` column.
- Replace `load_complex_annotation` + `build_complex_groups` (which did an internal GO groupby) with
  `load_long_table(path)` + `build_groups(background, long_table, min_group_size, max_group_size, max_term_genes)`:
  merge background onto the long-table on `Systematic ID`, group by `group_id`, keep groups where
  `min_group_size <= (#DR>thr members) <= max_group_size` AND `n_group_genes <= max_term_genes`.
- Output `coherence_metrics.tsv`: add `source`, `group_id`, `n_group_genes`; rename `complex`->`group_name`.
  Keep all existing metric columns + `observed_mpd`/`z_score`/`p_value`/`n_permutations` untouched.
- Byte-faithful core (`coherence_metrics` local fn, `compute_distance_zscore` call,
  `geometric_median`, normalization) unchanged.

**Step 1: Write the failing test** (filter logic is the new, testable surface)

```python
# tests/test_coherence_compute.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import numpy as np
import pandas as pd


def _background(ids):
    return pd.DataFrame({
        "Systematic ID": ids,
        "norm_DR": np.linspace(0.4, 0.9, len(ids)),
        "norm_DL": np.linspace(0.1, 0.5, len(ids)),
    })


def _long(rows):
    return pd.DataFrame(rows, columns=["source", "group_id", "group_name",
                                       "Systematic ID", "Name", "n_group_genes"])


def test_build_groups_respects_min_and_max_size():
    from compute_coherence import build_groups
    bg = _background([f"g{i}" for i in range(10)])
    long = _long([("go_cc", "GO:1", "big", f"g{i}", f"g{i}", 4) for i in range(4)]
                 + [("go_cc", "GO:2", "tiny", "g5", "g5", 1)])
    groups = build_groups(bg, long, min_group_size=3, max_group_size=300, max_term_genes=500)
    assert "GO:1" in groups          # 4 members, ok
    assert "GO:2" not in groups      # only 1 member < min


def test_build_groups_drops_broad_term_by_max_term_genes():
    from compute_coherence import build_groups
    bg = _background([f"g{i}" for i in range(10)])
    # 3 DR-members but the term annotates 999 genes total -> broad parent, dropped.
    long = _long([("go_bp", "GO:9", "broad", f"g{i}", f"g{i}", 999) for i in range(3)])
    groups = build_groups(bg, long, min_group_size=3, max_group_size=300, max_term_genes=500)
    assert "GO:9" not in groups
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_coherence_compute.py -v`
Expected: FAIL — `cannot import name 'build_groups'` (until refactor renames the function).

**Step 3: Write minimal implementation**

Perform the refactor described above. `build_groups` signature:

```python
def build_groups(background, long_table, min_group_size, max_group_size, max_term_genes):
    merged = long_table.merge(
        background[["Systematic ID", "norm_DR", "norm_DL"]], on="Systematic ID", how="inner"
    )
    groups = {}
    for group_id, grp in merged.groupby("group_id"):
        grp = grp.drop_duplicates(subset="Systematic ID")
        n_total = int(grp["n_group_genes"].iloc[0])
        if n_total > max_term_genes:
            continue
        if min_group_size <= len(grp) <= max_group_size:
            groups[group_id] = grp
    return groups
```

Thread `group_id`/`group_name`/`n_group_genes`/`source` into the output rows in `compute_coherence_table`.

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_coherence_compute.py tests/test_complex_coherence.py -v`
Expected: PASS — new filter tests + the 8 byte-faithful regressions still green.

**Step 5: Commit**

---

## Task 5: Summary figure — centroid-position plot + coherence×X panels (in compute)

**Files:**
- Modify: `workflow/scripts/coherence/compute_coherence.py`
- Test: `tests/test_coherence_panels.py`

Replace the old `plot_coherence` (size hist + volcano) with a multi-panel figure written to
`coherence_analysis.pdf`. All panels operate on the already-computed metrics table + the
long-table + (optionally) the features table. Panels:
1. **Size histogram** + **z-score histogram** (2 hist axes).
2. **Centroid-position plot**: x=`centroid_x` (typical DR), y=`centroid_y` (typical DL),
   color=`z_score` (coolwarm_r), size ∝ `term_size`.
3. **Panel A — coherence × shared-subunit fraction**: per group, fraction of its DR-members
   that appear in >=1 OTHER group of the SAME source (computed from the long-table's
   `n_groups_per_gene` within source). Always drawn.
4. **Panel B — coherence × abundance uniformity** and **Panel D — coherence × conservation
   uniformity**: only if `--features` given. Per group, CV of member abundance /
   `evolutionary_rate`. Missing feature column -> skip that panel + `logger.warning`.

The testable pure helpers: `shared_subunit_fraction(long_table)` and
`member_feature_cv(long_table, features, column)` returning `{group_id: value}`.
Keep plotting itself thin; test the numeric helpers, not the matplotlib output.

**Step 1: Write the failing test**

```python
# tests/test_coherence_panels.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import numpy as np
import pandas as pd


def _long(rows):
    return pd.DataFrame(rows, columns=["source", "group_id", "group_name",
                                       "Systematic ID", "Name", "n_group_genes"])


def test_shared_subunit_fraction_counts_cross_group_members():
    from compute_coherence import shared_subunit_fraction
    # gA shared between GO:1 and GO:2; gB only in GO:1.
    long = _long([
        ("go_cc", "GO:1", "one", "gA", "gA", 2),
        ("go_cc", "GO:1", "one", "gB", "gB", 2),
        ("go_cc", "GO:2", "two", "gA", "gA", 1),
    ])
    frac = shared_subunit_fraction(long)
    assert frac["GO:1"] == 0.5   # 1 of 2 members (gA) is shared
    assert frac["GO:2"] == 1.0   # gA shared


def test_member_feature_cv_computes_per_group_cv():
    from compute_coherence import member_feature_cv
    long = _long([
        ("go_cc", "GO:1", "one", "gA", "gA", 2),
        ("go_cc", "GO:1", "one", "gB", "gB", 2),
    ])
    features = pd.DataFrame({"Systematic ID": ["gA", "gB"], "abundance": [10.0, 30.0]})
    cv = member_feature_cv(long, features, "abundance")
    # CV = std/mean of [10,30] = 10/20 = 0.5 (population std) — assert close, ddof-agnostic-ish
    assert 0.4 < cv["GO:1"] < 0.8


def test_member_feature_cv_missing_column_returns_empty():
    from compute_coherence import member_feature_cv
    long = _long([("go_cc", "GO:1", "one", "gA", "gA", 1)])
    features = pd.DataFrame({"Systematic ID": ["gA"], "other": [1.0]})
    assert member_feature_cv(long, features, "abundance") == {}
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_coherence_panels.py -v`
Expected: FAIL — `cannot import name 'shared_subunit_fraction'`.

**Step 3: Write minimal implementation**

Add the two helpers + rewrite `plot_coherence(table, long_table, features=None)` to draw the panels.
`shared_subunit_fraction`: within each source, count how many of a group's members appear in other
groups; divide by the group's member count. `member_feature_cv`: merge features onto long-table by
`Systematic ID`, group by `group_id`, return `std/mean` per group; return `{}` if column absent.

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_coherence_panels.py -v`
Expected: PASS (3 tests).

**Step 5: Commit**

```bash
git add workflow/scripts/coherence/compute_coherence.py tests/test_coherence_panels.py
git commit -m "feat(coherence): centroid-position plot + coherence-vs-biology panels (A/B/D)"
```

---

## Task 6: `git mv` + rewrite module viz -> `plot_group_scatter.py` (stage [3])

**Files:**
- Rename: `workflow/scripts/complex/analyze_complex_modules.py` -> `workflow/scripts/coherence/plot_group_scatter.py` (via `git mv`)
- Modify: the renamed script (full rewrite of logic, keep header/conventions)
- Test: `tests/test_coherence_scatter.py`

**Step 0: Rename first (keep blame)**

```bash
git mv workflow/scripts/complex/analyze_complex_modules.py workflow/scripts/coherence/plot_group_scatter.py
git commit -m "refactor(coherence): git mv analyze_complex_modules.py -> coherence/plot_group_scatter.py"
```

**What it becomes:** a generic per-group scatter driven by a config namelist. Inputs:
`--fitting-results`, `--annotation` (long-table), `--metrics` (Task 4 TSV), `--groups`
(a dict-literal `{source: [names_or_ids]}` rendered by Snakemake, parsed like the old
`parse_modules_arg` via `ast.literal_eval`), `--source`, `--output-figure`. For the given
source, resolve each namelist entry to a `group_id` (match on `group_name` OR `group_id`),
then one subplot per resolved group: background cloud + members highlighted via
`plot_given_genes_on_feature_space`, annotated with `group_name (group_id)`, n_members, and
the group's `z_score`/`p_value` looked up from the metrics table.

Testable pure helper: `resolve_groups(long_table, source, names) -> list[(group_id, group_name, [systematic_ids])]`.

**Step 1: Write the failing test**

```python
# tests/test_coherence_scatter.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import pandas as pd


def _long():
    return pd.DataFrame({
        "source": ["go_cc"] * 3,
        "group_id": ["GO:1", "GO:1", "GO:2"],
        "group_name": ["alpha", "alpha", "beta"],
        "Systematic ID": ["gA", "gB", "gC"],
        "Name": ["gA", "gB", "gC"],
        "n_group_genes": [2, 2, 1],
    })


def test_resolve_groups_by_name():
    from plot_group_scatter import resolve_groups
    out = resolve_groups(_long(), "go_cc", ["alpha"])
    assert len(out) == 1
    gid, gname, genes = out[0]
    assert gid == "GO:1" and gname == "alpha" and set(genes) == {"gA", "gB"}


def test_resolve_groups_by_id():
    from plot_group_scatter import resolve_groups
    out = resolve_groups(_long(), "go_cc", ["GO:2"])
    assert out[0][0] == "GO:2"


def test_resolve_groups_skips_unknown():
    from plot_group_scatter import resolve_groups
    out = resolve_groups(_long(), "go_cc", ["nonexistent"])
    assert out == []
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_coherence_scatter.py -v`
Expected: FAIL — `cannot import name 'resolve_groups'`.

**Step 3: Write minimal implementation**

Rewrite the script; `resolve_groups` filters the long-table to `source`, then for each requested
name matches rows where `group_name == entry` or `group_id == entry`, returning
`(group_id, group_name, sorted(member systematic ids))`; skip entries with no match.

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_coherence_scatter.py -v`
Expected: PASS (3 tests).

**Step 5: Commit**

---

## Task 7: Rename the regression test file + config `complex:` -> `coherence:`

**Files:**
- Rename: `tests/test_complex_coherence.py` -> `tests/test_coherence.py` (via `git mv`)
- Modify: `config/analysis.yaml`

**Step 1: git mv the test (keep blame), then update config**

```bash
git mv tests/test_complex_coherence.py tests/test_coherence.py
```

Replace the `complex:` block (lines ~99-117) in `config/analysis.yaml` with:

```yaml
# --- Coherence of gene groups in DR-DL fitness space (coherence.smk) ---
# Per (dataset x source): for every group (GO complex / CC / BP term) whose
# DR>dr_threshold members number min..max AND whose total annotated size is
# <= max_term_genes, measure how tightly they cluster in normalized (DR, DL/10)
# space vs a seeded permutation null (median pairwise distance z-score). Fans out
# by `source`; add a database = add a loader in workflow/src/coherence/sources.py
# + a name here. scatter_groups drives plot_group_scatter (per-source namelist,
# entries match group_name OR group_id).
coherence:
  sources: [go_macrocomplex, go_cc, go_bp]
  min_group_size: 3
  max_group_size: 300
  max_term_genes: 500
  dr_threshold: 0.3
  n_permutations: 1000
  random_state: 42
  features_panels: true
  scatter_groups:
    go_macrocomplex: [kinetochore, "mitochondrial large ribosomal subunit"]
    go_cc: []
    go_bp: []
```

**Step 2: Verify tests still collect under the new name**

Run: `pytest tests/test_coherence.py -v`
Expected: PASS — the 8 byte-faithful regressions, now under the renamed file.

**Step 3: Commit**

```bash
git add config/analysis.yaml tests/test_coherence.py
git commit -m "refactor(coherence): rename test file + config complex: -> coherence:"
```

---

## Task 8: `git mv` complex.smk -> coherence.smk; write the three rules; update Snakefile

**Files:**
- Rename: `workflow/rules/complex.smk` -> `workflow/rules/coherence.smk` (via `git mv`)
- Modify: the renamed rules file (full rewrite)
- Modify: `Snakefile` (include line + rule all comment)

**Step 0: Rename first**

```bash
git mv workflow/rules/complex.smk workflow/rules/coherence.smk
git commit -m "refactor(coherence): git mv complex.smk -> coherence.smk"
```

**Step 1: Rewrite `coherence.smk`** with three rules fanning out over `{source}`. Reference the
existing enrichment.smk/clustering.smk patterns (wildcard constraints, `resources/external/pombase/{ref}`,
`config.get`). Sketch:

```python
# coherence.smk — gene-group coherence in DR-DL fitness space (per dataset x source)
_COH = "results/coherence/{dataset}/{source}"
_COH_CFG = config.get("coherence", {})
_COH_SOURCES = _COH_CFG.get("sources", ["go_macrocomplex", "go_cc", "go_bp"])

wildcard_constraints:
    source="|".join(_COH_SOURCES),

rule prepare_coherence_annotation:
    input:
        pombase_dir=f"resources/external/pombase/{DATASETS['reference']['pombase_version']}",
    output:
        long_table=f"{_COH}/group_annotation_long.tsv",
    log: "logs/coherence/prepare_{dataset}_{source}.log"
    conda: "../envs/biopython.yml"   # goatools lives here (go_cc/go_bp)
    shell:
        "python workflow/scripts/coherence/prepare_annotation.py "
        "--source {wildcards.source} --pombase-dir {input.pombase_dir} "
        "--output {output.long_table} &> {log}"

rule compute_coherence:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"),
        annotation=f"{_COH}/group_annotation_long.tsv",
        features=f"results/features/{DATASETS['reference']['pombase_version']}/pombe_coding_gene_protein_features.tsv",
    output:
        metrics=f"{_COH}/coherence_metrics.tsv",
        figure=f"{_COH}/coherence_analysis.pdf",
    params:
        min_size=_COH_CFG.get("min_group_size", 3),
        max_size=_COH_CFG.get("max_group_size", 300),
        max_term_genes=_COH_CFG.get("max_term_genes", 500),
        dr_threshold=_COH_CFG.get("dr_threshold", 0.3),
        n_permutations=_COH_CFG.get("n_permutations", 1000),
        random_state=_COH_CFG.get("random_state", 42),
        features_flag=lambda wc: ("--features results/features/"
            f"{DATASETS['reference']['pombase_version']}/pombe_coding_gene_protein_features.tsv"
            if _COH_CFG.get("features_panels", True) else ""),
    log: "logs/coherence/compute_{dataset}_{source}.log"
    conda: "../envs/statistics_and_figure_plotting.yml"
    shell:
        "python workflow/scripts/coherence/compute_coherence.py "
        "--fitting-results {input.fitting_results} --annotation {input.annotation} "
        "--source {wildcards.source} --min-size {params.min_size} --max-size {params.max_size} "
        "--max-term-genes {params.max_term_genes} --dr-threshold {params.dr_threshold} "
        "--n-permutations {params.n_permutations} --random-state {params.random_state} "
        "{params.features_flag} --output-metrics {output.metrics} --output-figure {output.figure} &> {log}"

rule plot_coherence_group_scatter:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"),
        annotation=f"{_COH}/group_annotation_long.tsv",
        metrics=f"{_COH}/coherence_metrics.tsv",
    output:
        figure=f"{_COH}/group_scatter.pdf",
    params:
        groups=lambda wc: _COH_CFG.get("scatter_groups", {}).get(wc.source, []),
    log: "logs/coherence/scatter_{dataset}_{source}.log"
    conda: "../envs/statistics_and_figure_plotting.yml"
    shell:
        "python workflow/scripts/coherence/plot_group_scatter.py "
        "--fitting-results {input.fitting_results} --annotation {input.annotation} "
        "--metrics {input.metrics} --source {wildcards.source} "
        "--groups \"{params.groups}\" --output-figure {output.figure} &> {log}"
```

Notes:
- `compute_coherence` lists `features` as an input so the DAG rebuilds it if features change;
  the `features_flag` param passes the same path only when `features_panels` is true. If you prefer
  no hard features dependency, drop it from `input:` and keep only the param — decide based on whether
  the features table is expected to always exist (it is a normal Snakemake target, so listing it is fine).
- Confirm the release-dir key: `DATASETS['datasets'][wc.dataset]['release_dir']` matches how
  clustering.smk/utr.smk reference `fitting_results.tsv`.

**Step 2: Update `Snakefile`**
- Line 41: `include: "workflow/rules/complex.smk"` -> `include: "workflow/rules/coherence.smk"`.
- In `rule all`, replace the commented `# f"results/complex/{_DATASET}/complex_coherence_metrics.tsv",`
  with a commented coherence target, e.g.:
  `# expand(f"results/coherence/{_DATASET}/{{source}}/coherence_metrics.tsv", source=config["coherence"]["sources"]),`

**Step 3: Dry-run the DAG**

Run:
```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate snakemake
snakemake -n results/coherence/HD_DIT_HAP/go_macrocomplex/coherence_metrics.tsv 2>&1 | tail -30
```
Expected: a clean DAG (jobs for prepare_coherence_annotation + compute_coherence, chaining off
features/fitting_results). If it reports the pombase resource or fitting_results missing, that is a
data-availability issue (inputs not present in this checkout), not a rule error — the rule graph
resolving without wildcard/syntax errors is the pass criterion. Note the Snakefile hardcodes
`workdir` to the MAIN tree (see [[followup-analysis-expansion]]), so a dry-run reflects the main
tree's files; that's expected.

**Step 4: Commit**

```bash
git add workflow/rules/coherence.smk Snakefile
git commit -m "feat(coherence): three-rule coherence.smk fanning out by source; wire Snakefile"
```

---

## Task 9: Remove the dead complex test scaffolding + full verification

**Files:**
- Delete (if now empty): `workflow/scripts/complex/` directory
- Verify: whole coherence test set + baseline unaffected

**Step 1: Confirm the old complex script dir is empty and remove it**

```bash
ls workflow/scripts/complex/    # should be empty after the two git mv's
git rm -r workflow/scripts/complex/ 2>/dev/null || rmdir workflow/scripts/complex
```

**Step 2: Run the full coherence test set**

Run:
```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate statistics_and_figure_plotting
pytest tests/test_coherence.py tests/test_coherence_sources.py tests/test_coherence_prepare.py \
       tests/test_coherence_compute.py tests/test_coherence_panels.py tests/test_coherence_scatter.py \
       tests/test_domain_differences.py -v
```
Expected: all PASS (goatools loader tests SKIP in this env). The 8 original byte-faithful
regressions in `test_coherence.py` must remain green.

**Step 3: Verify goatools loaders in the biopython env**

Run:
```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate <env-with-goatools>
pytest tests/test_coherence_sources.py -v
```
Expected: the two `go_cc`/`go_bp` tests now PASS (no longer skipped). If no goatools env exists,
note it and rely on the smk conda env at run time.

**Step 4: Confirm no lingering references to the old names**

Run: `grep -rn "complex_coherence\|analyze_complex_modules\|results/complex\|complex.smk" workflow/ Snakefile config/ tests/ | grep -v docs/`
Expected: no matches (all migrated). Fix any stragglers.

**Step 5: Final commit**

```bash
git add -A
git commit -m "chore(coherence): remove empty complex/ script dir; migration complete"
```

---

## Done criteria

- `results/coherence/{dataset}/{source}/` produces `group_annotation_long.tsv`,
  `coherence_metrics.tsv`, `coherence_analysis.pdf`, `group_scatter.pdf` for each of
  `go_macrocomplex` / `go_cc` / `go_bp`.
- Adding a new database is: write a loader in `sources.py`, add it to `SOURCE_LOADERS`,
  add its name to `config.coherence.sources`. No change to compute/plot.
- Byte-faithful coherence math unchanged; original regressions green.
- `attribution.py` still unwired (deferred to next PR).
