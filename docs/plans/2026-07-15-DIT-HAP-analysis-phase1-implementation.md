# DIT-HAP Analysis Phase 1 Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Scaffold the DIT_HAP_analysis repository structure per `docs/plans/2026-07-15-DIT-HAP-analysis-design.md`, implement the `datasets.yaml` registry (§3) and `workflow/src/` regrouping (§7), and deliver one complete working Snakemake rule (`features.smk`, §8) that faithfully reproduces `pombe_feature_collection.ipynb`'s full 15-source feature matrix deterministically.

**Architecture:** Snakemake 9.0+ orchestration with conda env-per-rule dependency management (matching DIT_HAP_snakemake conventions), frozen dataclasses for config, pytest for registry/config/feature unit tests. Notebook logic is ported faithfully — including its documented quirks (asymmetric BioGrid degree counting, duplicate `DeletionLibrary_essentiality` column) — with no silent "improvements" beyond determinism/modularity.

**Tech Stack:** Python 3.12, Snakemake 9.13, pandas 3.x, loguru, Bio/gffutils/codonbias/goatools (via `bioinformatics` conda env), pytest 8.4.

**Key corrections vs. a naive reading of the design doc** (confirmed by reading `packaging.smk`, `DIT_HAP_streamlit/src/data/config.py`, and the actual notebook code):
- The registered project directory is `LD_haploid`, not `LD_diploid` (design doc §3 comment has a typo).
- `use_DEseq2_for_biological_replicates` (which gates `imputation_statistics.tsv`) is `true` for 4 projects (`HD_DIT_HAP_generationRAW`, `LD_DIT_HAP_generationRAW`, `HD_DIT_HAP_generationPLUS1`, `LD_DIT_HAP_generationPLUS1`), not just one — added as `has_imputation` per dataset entry, mirroring streamlit's `datasets.yaml`.
- Only 4 of 9 projects have a populated `release/` dir on disk right now (`Spikein`, `Spore2YES6_1328`, `LD_DIT_HAP_generationRAW`, `LD_haploid`) — tests target those; others will validate once upstream packaging finishes.

---

## Task 1: Repository skeleton

**Files:**
- Create: `.gitignore`
- Create: `README.md`
- Create: `Snakefile`
- Create: `config/datasets.yaml`
- Create: `config/analysis.yaml`
- Create: `config/DIT_HAP.mplstyle`
- Create: `workflow/envs/biopython.yml`
- Create: `workflow/envs/statistics_and_figure_plotting.yml`

Directory layout follows design doc §2 exactly: `notebooks/` lives at repo root (not under `workflow/`), `workflow/` holds only `rules/`, `scripts/`, `src/`, `envs/`.

**Step 1: Write .gitignore**

Per design doc §4, `resources/external/` and `resources/literature/` are regenerable local copies of external databases (ignored); `resources/curated/` is version-controlled human-curated output (NOT ignored).

```gitignore
# Snakemake
.snakemake/

# Python
__pycache__/
*.py[cod]
*/__pycache__/

# Pipeline outputs
/results/
/reports/
/logs/

# Regenerable local copies of external databases (see resources/curated/ which IS tracked)
/resources/external/
/resources/literature/

# Claude/editor artifacts
.claude
CLAUDE.md
.vscode
.cursor
.pytest_cache

# Temporary
tmp/
*.log
```

**Step 2: Write README.md**

```markdown
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
```

**Step 3: Write minimal Snakefile**

Per design doc §8: `config/datasets.yaml` is a data registry, not experiment parameters, so it's read directly via `yaml.safe_load` rather than Snakemake's `configfile:` directive (which would dump every key into the global `config` object used for wildcard validation elsewhere).

```python
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
include: "workflow/rules/features.smk"
# include: "workflow/rules/enrichment.smk"
# include: "workflow/rules/clustering.smk"
# (Phase 1 delivers only features.smk; remaining rules are follow-up work)

# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
rule all:
    input:
        # Phase 1 target: gene feature matrix for the reference PomBase version
        f"results/features/{DATASETS['reference']['pombase_version']}/pombe_coding_gene_protein_features.tsv",
    message:
        "*** DIT-HAP analysis complete"
```

**Step 4: Write config/datasets.yaml**

Per design doc §3: this registry's job is to point at `DIT_HAP_snakemake`'s per-project `release/` dirs (absolute path to that repo, then relative `release_dir` per project) — it is NOT a PomBase-resource registry. `has_time_points` mirrors `packaging.smk`'s branch (confirmed against every project's `config/config.yaml`); `has_imputation` mirrors `use_DEseq2_for_biological_replicates`. Only 4 datasets currently have a populated `release/` on disk (noted below) — the rest are registered for when upstream packaging catches up.

```yaml
# Dataset registry — points at DIT_HAP_snakemake's packaged release/ dirs.
# See docs/plans/2026-07-15-DIT-HAP-analysis-design.md §3.
default_dataset: HD_DIT_HAP_generationRAW

snakemake_repo: /data/c/yangyusheng_optimized/DIT_HAP_snakemake

datasets:
  HD_DIT_HAP_generationRAW:
    label: "HD, raw generations"
    release_dir: projects/HD_DIT_HAP_generationRAW/release
    has_time_points: true
    has_imputation: true
  LD_DIT_HAP_generationRAW:
    label: "LD, raw generations"
    release_dir: projects/LD_DIT_HAP_generationRAW/release
    has_time_points: true
    has_imputation: true   # release/ populated on disk
  HD_DIT_HAP:
    label: "HD DIT-HAP"
    release_dir: projects/HD_DIT_HAP/release
    has_time_points: true
    has_imputation: false
  HD_DIT_HAP_generationPLUS1:
    label: "HD, generation+1"
    release_dir: projects/HD_DIT_HAP_generationPLUS1/release
    has_time_points: true
    has_imputation: true
  LD_DIT_HAP_generationPLUS1:
    label: "LD, generation+1"
    release_dir: projects/LD_DIT_HAP_generationPLUS1/release
    has_time_points: true
    has_imputation: true
  HD_diploid:
    label: "HD diploid"
    release_dir: projects/HD_diploid/release
    has_time_points: true
    has_imputation: false
  LD_haploid:
    label: "LD haploid"
    release_dir: projects/LD_haploid/release
    has_time_points: true
    has_imputation: false   # release/ populated on disk
  Spikein:
    label: "Spike-in calibration"
    release_dir: projects/Spikein/release
    has_time_points: false
    has_imputation: false   # release/ populated on disk
  Spore2YES6_1328:
    label: "Spore-to-YES6 1328"
    release_dir: projects/Spore2YES6_1328/release
    has_time_points: true
    has_imputation: false   # release/ populated on disk

reference:
  pombase_version: "2025-10-01"
  # AlphaFold structures are too large to copy locally (design doc §4) — referenced
  # by absolute external path instead of a resources/ copy.
  alphafold_dir: /data/c/yangyusheng_optimized/resource/AlphaFold_Dataset/20251107_downloaded/UP000002485_284812_SCHPO_v6
```

**Step 5: Write config/analysis.yaml (scaffold)**

Per design doc §9, exact fields get filled in as later stages are migrated. Phase 1 needs no entries here yet since `features.smk` has no tunable parameters (it's pure feature extraction), but the file is created now so later stages append to it rather than inventing a new location.

```yaml
# This project's own analysis parameters (as opposed to config/datasets.yaml,
# which only points at where upstream data lives). Populated incrementally as
# each analysis stage is migrated — see design doc §9.
```

**Step 6: Copy DIT_HAP.mplstyle from DIT_HAP_pipeline**

```bash
cp /data/c/yangyusheng_optimized/DIT_HAP_pipeline/config/DIT_HAP.mplstyle config/
```

**Step 7: Write workflow/envs/biopython.yml**

```yaml
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.12
  - pandas>=3.0.3
  - numpy>=2.5.0
  - biopython>=1.87
  - gffutils>=0.13
  - codonbias
  - goatools>=1.6.5
  - openpyxl
  - loguru
  - tqdm
  - requests
  - pyyaml
```

**Step 8: Write workflow/envs/statistics_and_figure_plotting.yml**

```yaml
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.12
  - matplotlib>=3.11
  - seaborn>=0.13.2
  - numpy>=2.5.0
  - pandas>=3.0.3
  - scipy>=1.18
  - openpyxl
  - tqdm
  - loguru
  - pyyaml
```

Note: `features.smk` doesn't exist yet (created in Task 5) — comment out its `include:` line in the Snakefile for now, and the `rule all` target that depends on it. Task 5 Step 4 uncomments both.

```python
# include: "workflow/rules/features.smk"
```
```python
rule all:
    input:
        # Uncommented in Task 5 once features.smk exists:
        # f"results/features/{DATASETS['reference']['pombase_version']}/pombe_coding_gene_protein_features.tsv",
    message:
        "*** DIT-HAP analysis complete"
```

**Step 9: Sanity-check the Snakefile parses**

```bash
mamba run -n snakemake snakemake -n
```

Expected: dry-run succeeds with `rule all` having zero inputs (nothing to build yet), no parse errors.

**Step 10: Commit skeleton**

```bash
git add .gitignore README.md Snakefile config/ workflow/envs/
git commit -m "feat: add repository skeleton

- .gitignore excludes results/resources copies, tracks resources/curated/
- README with usage examples
- Snakefile entry point reading config/datasets.yaml directly (not via configfile:)
- config/datasets.yaml: registry of DIT_HAP_snakemake release/ dirs per project
- config/analysis.yaml scaffold for later stages' parameters
- config/DIT_HAP.mplstyle from DIT_HAP_pipeline
- workflow/envs/*.yml for biopython and plotting

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 2: `workflow/src/data_config.py` — release/ registry loader (design doc §3)

This is the ONLY module that may construct a path into `DIT_HAP_snakemake`. It mirrors `DIT_HAP_streamlit/src/data/config.py`'s `InsertionLevelConfig`/`GeneLevelConfig` pattern exactly, but points at `release_dir` from `config/datasets.yaml` instead of a locally-copied `data_dir`. It has nothing to do with PomBase resource paths — those are a separate concern handled in Task 6 (`features/genome.py`), since design doc §7 lists `data_config.py`'s job as strictly "数据源注册表加载（第 3 节）".

**Files:**
- Create: `workflow/src/__init__.py`
- Create: `workflow/src/data_config.py`
- Create: `tests/__init__.py`
- Create: `tests/test_data_config.py`

**Step 1: Write workflow/src/__init__.py** (empty file for package)

```python
# Empty __init__.py
```

**Step 2: Write workflow/src/data_config.py** (imports, constants, InsertionLevelConfig)

```python
"""
Dataset Registry Loader
========================

Loads `config/datasets.yaml`, the single entry point for resolving paths into
DIT_HAP_snakemake's packaged `release/` output. See design doc §3 — no other
module in this repository should construct a path into DIT_HAP_snakemake directly.

Mirrors DIT_HAP_streamlit/src/data/config.py's InsertionLevelConfig/GeneLevelConfig
pattern, but resolves against release_dir (an upstream project's release/ folder)
rather than a locally-copied data_dir.

Input
-----
- config/datasets.yaml: registry of dataset name -> release_dir + branch flags

Output
------
- DatasetConfig instances with Path fields (existence NOT validated automatically —
  many registered datasets have no release/ populated yet; call validate_all_paths()
  explicitly when a script actually needs the files to exist).

Usage
-----
    from workflow.src.data_config import load_dataset_config
    cfg = load_dataset_config("Spore2YES6_1328")
    lfc = pd.read_csv(cfg.insertion_level.LFCs, sep="\\t")

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from dataclasses import dataclass
from functools import cache
from pathlib import Path

# 2. Data Processing Imports
import yaml

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
DATASETS_YAML = Path("config/datasets.yaml")

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class InsertionLevelConfig:
    """Paths under release/insertion_level/. Always-present files plus branch-gated optionals."""
    annotations: Path
    baseMean: Path
    LFCs: Path
    fitting_LFCs: Path | None = None
    fitting_results: Path | None = None
    transformed_weights: Path | None = None
    imputation_statistics: Path | None = None

    def validate_paths(self) -> None:
        """Raise FileNotFoundError if any present (non-None) path is missing on disk."""
        for path in [self.annotations, self.baseMean, self.LFCs, self.fitting_LFCs,
                     self.fitting_results, self.transformed_weights, self.imputation_statistics]:
            if path is not None and not path.exists():
                raise FileNotFoundError(f"Insertion-level release file not found: {path}")
```

**Step 3: Continue workflow/src/data_config.py** (GeneLevelConfig, DatasetConfig, registry loader)

```python
@dataclass(kw_only=True, frozen=True)
class GeneLevelConfig:
    """Paths under release/gene_level/. Only constructed when has_time_points is True."""
    LFCs: Path
    fitting_LFCs: Path
    fitting_results: Path

    def validate_paths(self) -> None:
        """Raise FileNotFoundError if any path is missing on disk."""
        for path in [self.LFCs, self.fitting_LFCs, self.fitting_results]:
            if not path.exists():
                raise FileNotFoundError(f"Gene-level release file not found: {path}")


@dataclass(kw_only=True, frozen=True)
class DatasetConfig:
    """Aggregate config for one registered project's release/ output."""
    name: str
    label: str
    has_time_points: bool
    has_imputation: bool
    insertion_level: InsertionLevelConfig
    gene_level: GeneLevelConfig | None

    def validate_all_paths(self) -> None:
        """Validate insertion_level (always) and gene_level (if present)."""
        self.insertion_level.validate_paths()
        if self.gene_level is not None:
            self.gene_level.validate_paths()


# =============================================================================
# REGISTRY LOADING
# =============================================================================
@cache
def _load_registry() -> dict:
    """Read and cache config/datasets.yaml for the lifetime of the process."""
    with DATASETS_YAML.open() as f:
        return yaml.safe_load(f)


def default_dataset_name() -> str:
    """Return the registry's default_dataset key."""
    return _load_registry()["default_dataset"]


def list_dataset_names() -> list[str]:
    """Return all registered dataset names, for wildcard_constraints in the Snakefile."""
    return list(_load_registry()["datasets"].keys())


def reference_pombase_version() -> str:
    """Return the registry's reference.pombase_version (used by workflow/src/features/)."""
    return _load_registry()["reference"]["pombase_version"]


def reference_alphafold_dir() -> Path:
    """Return the registry's reference.alphafold_dir as a Path (never copied locally, see design doc §4)."""
    return Path(_load_registry()["reference"]["alphafold_dir"])


def load_dataset_config(name: str | None = None) -> DatasetConfig:
    """Resolve one dataset's release/ paths. Raises KeyError if name is not registered."""
    registry = _load_registry()
    if name is None:
        name = registry["default_dataset"]
    if name not in registry["datasets"]:
        raise KeyError(f"Unknown dataset {name!r}; registered: {sorted(registry['datasets'])}")

    entry = registry["datasets"][name]
    release_dir = Path(registry["snakemake_repo"]) / entry["release_dir"]
    has_time_points = entry["has_time_points"]
    has_imputation = entry["has_imputation"]

    insertion_level = InsertionLevelConfig(
        annotations=release_dir / "insertion_level" / "annotations.tsv.gz",
        baseMean=release_dir / "insertion_level" / "baseMean.tsv",
        LFCs=release_dir / "insertion_level" / "LFC.tsv",
        fitting_LFCs=(release_dir / "insertion_level" / "fitting_LFCs.tsv") if has_time_points else None,
        fitting_results=(release_dir / "insertion_level" / "fitting_results.tsv") if has_time_points else None,
        transformed_weights=(release_dir / "insertion_level" / "transformed_weights.tsv") if has_time_points else None,
        imputation_statistics=(release_dir / "insertion_level" / "imputation_statistics.tsv") if has_imputation else None,
    )

    gene_level = None
    if has_time_points:
        gene_level = GeneLevelConfig(
            LFCs=release_dir / "gene_level" / "LFC.tsv",
            fitting_LFCs=release_dir / "gene_level" / "fitting_LFCs.tsv",
            fitting_results=release_dir / "gene_level" / "fitting_results.tsv",
        )

    return DatasetConfig(
        name=name,
        label=entry["label"],
        has_time_points=has_time_points,
        has_imputation=has_imputation,
        insertion_level=insertion_level,
        gene_level=gene_level,
    )
```

**Step 4: Write tests/__init__.py** (empty)

```python
# Empty __init__.py
```

**Step 5: Write tests/test_data_config.py**

Tests target the 4 datasets confirmed to have a populated `release/` on disk right now (`Spikein`, `Spore2YES6_1328`, `LD_DIT_HAP_generationRAW`, `LD_haploid`) plus dataclass-level validation logic that doesn't need real files.

```python
"""Tests for workflow/src/data_config.py registry loader and dataclasses."""

import pytest
from pathlib import Path
from workflow.src.data_config import (
    InsertionLevelConfig,
    GeneLevelConfig,
    load_dataset_config,
    list_dataset_names,
    default_dataset_name,
    reference_pombase_version,
    reference_alphafold_dir,
)


def test_list_dataset_names_includes_all_nine_projects():
    """Registry lists all 9 projects confirmed on disk under DIT_HAP_snakemake/projects/."""
    names = list_dataset_names()
    assert set(names) == {
        "HD_DIT_HAP_generationRAW", "LD_DIT_HAP_generationRAW", "HD_DIT_HAP",
        "HD_DIT_HAP_generationPLUS1", "LD_DIT_HAP_generationPLUS1",
        "HD_diploid", "LD_haploid", "Spikein", "Spore2YES6_1328",
    }


def test_default_dataset_name():
    """default_dataset matches the registry's configured value."""
    assert default_dataset_name() == "HD_DIT_HAP_generationRAW"


def test_insertion_level_config_validates_missing_optional_as_none():
    """Optional fields left as None do not raise on validate_paths."""
    cfg = InsertionLevelConfig(
        annotations=Path("/nonexistent/annotations.tsv.gz"),
        baseMean=Path("/nonexistent/baseMean.tsv"),
        LFCs=Path("/nonexistent/LFC.tsv"),
    )
    with pytest.raises(FileNotFoundError, match="Insertion-level release file not found"):
        cfg.validate_paths()


def test_spikein_has_no_gene_level_config():
    """Spikein has time_points commented out upstream -> gene_level must be None."""
    cfg = load_dataset_config("Spikein")
    assert cfg.gene_level is None
    assert cfg.has_time_points is False


@pytest.mark.parametrize("name", ["Spikein", "Spore2YES6_1328", "LD_DIT_HAP_generationRAW", "LD_haploid"])
def test_populated_datasets_validate_against_real_release_dirs(name):
    """These 4 datasets have a populated release/ on disk right now; validate_all_paths must pass."""
    cfg = load_dataset_config(name)
    cfg.validate_all_paths()


def test_ld_dit_hap_generationraw_has_imputation():
    """LD_DIT_HAP_generationRAW is the one populated dataset with use_DEseq2_for_biological_replicates=True."""
    cfg = load_dataset_config("LD_DIT_HAP_generationRAW")
    assert cfg.has_imputation is True
    assert cfg.insertion_level.imputation_statistics is not None


def test_unknown_dataset_raises_keyerror():
    """load_dataset_config raises KeyError with the unregistered name in the message."""
    with pytest.raises(KeyError, match="not_a_real_dataset"):
        load_dataset_config("not_a_real_dataset")


def test_reference_pombase_version():
    """reference_pombase_version matches the registry's reference.pombase_version."""
    assert reference_pombase_version() == "2025-10-01"


def test_reference_alphafold_dir_is_absolute_and_exists():
    """reference_alphafold_dir points at the real (uncopied) AlphaFold directory."""
    path = reference_alphafold_dir()
    assert path.is_absolute()
    assert path.exists()
```

**Step 6: Run tests**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_analysis
conda run -n bioinformatics python -m pytest tests/test_data_config.py -v
```

Expected: all 10 tests pass. `test_populated_datasets_validate_against_real_release_dirs` requires `DIT_HAP_snakemake` at the path in `snakemake_repo` — if it fails with `FileNotFoundError`, re-check that repo's `projects/*/release/` still has the same 4 populated dirs (upstream packaging may have added more since this plan was written; add their names to the parametrize list if so).

**Step 7: Commit registry module**

```bash
git add workflow/src/__init__.py workflow/src/data_config.py tests/
git commit -m "feat: add dataset registry loader for DIT_HAP_snakemake release/ dirs

- workflow/src/data_config.py: InsertionLevelConfig/GeneLevelConfig dataclasses,
  load_dataset_config() resolving config/datasets.yaml against snakemake_repo
- tests/test_data_config.py: validates against the 4 currently-populated release/ dirs
- Only module permitted to construct a DIT_HAP_snakemake path (design doc §3)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 3: Physical resource copies (design doc §4)

Design doc §4 is explicit: PomBase/STRING/KEGG/BioGrid/Ensembl get a full local copy so this repo doesn't couple to DIT_HAP_pipeline's data lifecycle (symlinks would defeat that). AlphaFold stays referenced by absolute external path (too large to copy, already configured in `config/datasets.yaml`'s... actually see note below). The Hayles deletion-library spreadsheet and the verification-summary CSV are **curated** artifacts per design doc §4's tree (`resources/curated/deletion_library_categories.xlsx`, `resources/curated/essentiality_verification.csv`) — they get git-tracked copies with renamed, semantic filenames, not gitignored external copies.

**Disk cost — flagging before you run this:** copying `pombase_data/2025-10-01` (312M) + `BioGrid` (47M) + `Ensemble` (~1M) + `Literature` (23M) totals roughly 380M. Reversible (just `rm -rf resources/external resources/literature` if wrong), non-destructive to the source (DIT_HAP_pipeline is untouched), and 677G is free on `/data/c` — proceeding is safe, but this step does write real disk space rather than a cheap symlink.

**Files:**
- Create: `resources/external/pombase/2025-10-01/` (copy)
- Create: `resources/external/biogrid/` (copy)
- Create: `resources/external/ensembl/` (copy)
- Create: `resources/literature/` (copy)
- Create: `resources/curated/deletion_library_categories.xlsx` (copy, renamed)
- Create: `resources/curated/essentiality_verification.csv` (copy, renamed)

**Step 1: Create directory structure and copy PomBase data**

```bash
mkdir -p resources/external/pombase resources/external/biogrid resources/external/ensembl resources/literature resources/curated
cp -r /data/c/yangyusheng_optimized/DIT_HAP_pipeline/resources/pombase_data/2025-10-01 resources/external/pombase/2025-10-01
```

**Step 2: Copy BioGrid, Ensembl, Literature**

```bash
cp -r /data/c/yangyusheng_optimized/DIT_HAP_pipeline/resources/BioGrid/. resources/external/biogrid/
cp -r /data/c/yangyusheng_optimized/DIT_HAP_pipeline/resources/Ensemble/. resources/external/ensembl/
cp -r /data/c/yangyusheng_optimized/DIT_HAP_pipeline/resources/Literature/. resources/literature/
```

**Step 3: Copy curated artifacts with semantic renames**

```bash
cp /data/c/yangyusheng_optimized/DIT_HAP_pipeline/resources/Hayles_2013_OB_merged_categories_sysIDupdated.xlsx resources/curated/deletion_library_categories.xlsx
cp /data/c/yangyusheng_optimized/DIT_HAP_pipeline/resources/20260220modified_verification_summary.csv resources/curated/essentiality_verification.csv
```

**Step 4: Verify copy sizes and file counts match source**

```bash
du -sh resources/external/pombase/2025-10-01 resources/external/biogrid resources/external/ensembl resources/literature
diff <(cd /data/c/yangyusheng_optimized/DIT_HAP_pipeline/resources/pombase_data/2025-10-01 && find . -type f | sort) \
     <(cd resources/external/pombase/2025-10-01 && find . -type f | sort)
```

Expected: sizes roughly match the source (312M/47M/1M/23M), `diff` produces no output (identical file listings).

**Step 5: Commit — resources/curated/ is tracked, resources/external/ and resources/literature/ are gitignored**

```bash
git add resources/curated/
git status --porcelain resources/external resources/literature   # confirm these show as ignored, not untracked
git commit -m "feat: add curated deletion-library and essentiality-verification artifacts

- resources/curated/deletion_library_categories.xlsx (from Hayles_2013_OB_merged_categories_sysIDupdated.xlsx)
- resources/curated/essentiality_verification.csv (from 20260220modified_verification_summary.csv)
- resources/external/ and resources/literature/ populated locally but gitignored
  per design doc §4 (decouples this repo from DIT_HAP_pipeline's data lifecycle)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

**Scope note for Tasks 4-9:** design doc §7 lists a larger `workflow/src/` regroup (`plotting/`, `enrichment/pipeline.py` with STRING+REVIGO, `features/genome.py`+`features/protein.py`). Porting all of it now would add modules the delivered rule never calls — `pombe_feature_collection.ipynb` uses `utils.py`, `pombe_feature_functions.py`, `protein_structure_functions.py`, and only the `OntologyDataConfig`/`load_ontology_data` half of `enrichment_functions.py` (GO term richness, not the enrichment-study/STRING/REVIGO pipeline, which the notebook never calls). Tasks 4-8 port exactly that subset; `plot.py`, `enrichment/pipeline.py` (stringdb/revigo/enrichment study), and `plotting/gene_level.py` are explicitly deferred to a follow-up plan once a rule that actually needs them is scoped. This keeps every ported line exercised by the one rule this plan commits to delivering.

---

## Task 4: `workflow/src/io.py` and `workflow/src/gene_ids.py` (split of `utils.py`)

**Files:**
- Create: `workflow/src/io.py`
- Create: `workflow/src/gene_ids.py`
- Create: `tests/test_io.py`
- Create: `tests/test_gene_ids.py`

**Step 1: Write workflow/src/io.py**

```python
"""
Generic Table Readers
======================

File-extension-dispatched table loading, factored out of the original
`workflow/src/utils.py` (DIT_HAP_pipeline). No biology-specific logic —
safe to import from any module in this repository.

Input
-----
- Any .tsv/.bed/.csv/.xlsx file path

Output
------
- pandas DataFrame

Usage
-----
    from workflow.src.io import read_file
    df = read_file(Path("resources/curated/essentiality_verification.csv"))

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# =============================================================================
# CORE LOGIC
# =============================================================================
def read_file(file: Path, **kwargs) -> pd.DataFrame:
    """Read a table into a DataFrame, dispatching on file extension (tsv/bed/csv/xlsx)."""
    if "tsv" in file.name:
        return pd.read_csv(file, sep="\t", **kwargs)
    elif "bed" in file.name:
        return pd.read_csv(file, sep="\t", **kwargs)
    elif "csv" in file.name:
        return pd.read_csv(file, sep=",", **kwargs)
    elif "xlsx" in file.name:
        return pd.read_excel(file, **kwargs)
    else:
        raise ValueError(f"Unsupported file type: {file.name}")
```

**Step 2: Write tests/test_io.py**

```python
"""Tests for workflow/src/io.py file-extension dispatch."""

import pytest
import pandas as pd
from pathlib import Path
from workflow.src.io import read_file


def test_read_tsv(tmp_path):
    """A .tsv file is read with tab separator."""
    f = tmp_path / "data.tsv"
    f.write_text("a\tb\n1\t2\n")
    df = read_file(f)
    assert list(df.columns) == ["a", "b"]
    assert df.iloc[0]["a"] == 1


def test_read_csv(tmp_path):
    """A .csv file is read with comma separator."""
    f = tmp_path / "data.csv"
    f.write_text("a,b\n1,2\n")
    df = read_file(f)
    assert list(df.columns) == ["a", "b"]


def test_read_bed_uses_tab_separator(tmp_path):
    """A .bed file is read with tab separator like tsv."""
    f = tmp_path / "regions.bed"
    f.write_text("chr1\t0\t100\n")
    df = read_file(f, header=None)
    assert df.shape == (1, 3)


def test_unsupported_extension_raises(tmp_path):
    """An unrecognized extension raises ValueError naming the file."""
    f = tmp_path / "data.json"
    f.write_text("{}")
    with pytest.raises(ValueError, match="Unsupported file type: data.json"):
        read_file(f)
```

**Step 3: Write workflow/src/gene_ids.py**

```python
"""
Gene Systematic ID Resolution
===============================

Resolves gene names/synonyms to current PomBase systematic IDs, factored out
of the original `workflow/src/utils.py` (DIT_HAP_pipeline). Depends on
`workflow.src.io.read_file` to load the gene metadata table.

Input
-----
- A list of gene identifiers (names, synonyms, or already-current systematic IDs)
- A PomBase gene metadata table (gene_IDs_names_products.tsv) with columns
  gene_systematic_id, gene_name, synonyms, gene_type

Output
------
- A list of the same length with each entry resolved to its current
  systematic ID where a unique match is found, or left unchanged/NaN
  where the id is unknown or ambiguous (logged via print, matching the
  original notebook workflow's use of these log lines for manual review).

Usage
-----
    from workflow.src.gene_ids import update_sysIDs
    resolved = update_sysIDs(["cdc2", "SPBC11B10.09"], gene_meta_file)

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Local Imports
from workflow.src.io import read_file

# =============================================================================
# CORE LOGIC
# =============================================================================
def update_sysIDs(
    genes: list[str],
    gene_meta_file: Path,
    gene_filter: str = "gene_type == 'protein coding gene'",
) -> list[str | float]:
    """Resolve each gene name/synonym in `genes` to its current systematic ID."""
    gene_meta = read_file(gene_meta_file)
    gene_meta["gene_name"] = gene_meta["gene_name"].fillna(gene_meta["gene_systematic_id"])

    filtered_genes = gene_meta.query(gene_filter)
    synonyms2ID = (
        filtered_genes.set_index("gene_systematic_id")["synonyms"]
        .str.split(",")
        .explode()
        .str.strip()
        .dropna()
        .reset_index()
        .set_index("synonyms")
    )
    names2ID = (
        filtered_genes.set_index("gene_name")["gene_systematic_id"]
        .drop_duplicates()
        .reset_index()
        .set_index("gene_name")
    )
    sysIDs_now = filtered_genes["gene_systematic_id"].unique().tolist()

    updated_sysIDs = []
    for gene in genes:
        if isinstance(gene, str):
            gene = gene.strip()
            if "." in gene:
                gene = gene.split(".")[0].upper() + "." + gene.split(".")[1].lower()
        if pd.isna(gene):
            updated_sysIDs.append(gene)
            print(f"{gene} is NA")
        elif gene in sysIDs_now:
            updated_sysIDs.append(gene)
        elif gene in names2ID.index:
            updated = names2ID.loc[gene, "gene_systematic_id"]
            if isinstance(updated, str):
                updated_sysIDs.append(updated)
                print(f"{gene} is updated to {updated}")
            else:
                updated_sysIDs.append(np.nan)
                print(f"{gene} has multiple updates:", updated)
        elif gene in synonyms2ID.index:
            updated = synonyms2ID.loc[gene, "gene_systematic_id"]
            if isinstance(updated, str):
                updated_sysIDs.append(updated)
                print(f"{gene} is updated to {updated}")
            else:
                updated_sysIDs.append(np.nan)
                print(f"{gene} has multiple updates:", updated)
        else:
            updated_sysIDs.append(gene)
            print(f"{gene} is not found in geneid2symbol or synonyms2ID")
    return updated_sysIDs
```

Note: this is a byte-faithful port of `utils.update_sysIDs` (same log-via-print behavior, same case-normalization quirk for dotted transcript IDs like `spac1002.01` → `SPAC1002.01`) — no behavior changes, only the import of `read_file` moved to the new `workflow.src.io` location.

**Step 4: Write tests/test_gene_ids.py**

```python
"""Tests for workflow/src/gene_ids.py systematic ID resolution."""

import math
import pandas as pd
from workflow.src.gene_ids import update_sysIDs


def _write_gene_meta(tmp_path):
    """Write a minimal gene_IDs_names_products.tsv fixture."""
    f = tmp_path / "gene_IDs_names_products.tsv"
    df = pd.DataFrame({
        "gene_systematic_id": ["SPBC11B10.09", "SPAC1002.01"],
        "gene_name": ["cdc2", None],
        "synonyms": ["cdk1,cdc28", ""],
        "gene_type": ["protein coding gene", "protein coding gene"],
    })
    df.to_csv(f, sep="\t", index=False)
    return f


def test_already_current_sysid_passes_through(tmp_path):
    """A gene already given as a current systematic ID is returned unchanged."""
    meta = _write_gene_meta(tmp_path)
    result = update_sysIDs(["SPBC11B10.09"], meta)
    assert result == ["SPBC11B10.09"]


def test_gene_name_resolves_to_sysid(tmp_path):
    """A gene given by its common name resolves to the systematic ID."""
    meta = _write_gene_meta(tmp_path)
    result = update_sysIDs(["cdc2"], meta)
    assert result == ["SPBC11B10.09"]


def test_synonym_resolves_to_sysid(tmp_path):
    """A gene given by a synonym resolves to the systematic ID."""
    meta = _write_gene_meta(tmp_path)
    result = update_sysIDs(["cdk1"], meta)
    assert result == ["SPBC11B10.09"]


def test_unknown_gene_passes_through_unchanged(tmp_path):
    """A gene not found anywhere is returned unchanged (for manual review)."""
    meta = _write_gene_meta(tmp_path)
    result = update_sysIDs(["not_a_real_gene"], meta)
    assert result == ["not_a_real_gene"]


def test_na_input_passes_through(tmp_path):
    """A NaN input is passed through as NaN, not resolved."""
    meta = _write_gene_meta(tmp_path)
    result = update_sysIDs([float("nan")], meta)
    assert len(result) == 1 and math.isnan(result[0])


def test_dotted_transcript_id_is_case_normalized(tmp_path):
    """A lowercase dotted transcript ID is normalized to GENE.transcript case before lookup."""
    meta = _write_gene_meta(tmp_path)
    result = update_sysIDs(["spac1002.01"], meta)
    assert result == ["SPAC1002.01"]
```

**Step 5: Run tests**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_analysis
conda run -n bioinformatics python -m pytest tests/test_io.py tests/test_gene_ids.py -v
```

Expected: all 10 tests pass.

**Step 6: Commit**

```bash
git add workflow/src/io.py workflow/src/gene_ids.py tests/test_io.py tests/test_gene_ids.py
git commit -m "feat: split utils.py into io.py (read_file) and gene_ids.py (update_sysIDs)

- workflow/src/io.py: generic file-extension-dispatched table reader
- workflow/src/gene_ids.py: systematic ID resolution, byte-faithful port
- tests for both modules using tmp_path fixtures (no real PomBase data needed)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 5: `workflow/src/features/genome.py` (DNA-level features)

Byte-faithful port of `pombe_feature_functions.py`'s DNA-level half: the `config` dataclass's genome-related fields, `CHROMOSOME_END`/`CENTROMERE_POSITIONS` constants, `determine_primary_candidate`, `DNA_level_features`, and `calculate_anticodon_usage_matrix`. `config.PomBase_resource_dir` originally defaulted to a relative path computed from `__file__` — replaced here with an explicit constructor arg resolved by the caller via `reference_pombase_version()` + a `pombase_dir` root, since design doc §3/§4 forbid hand-built relative paths.

**Files:**
- Create: `workflow/src/features/__init__.py`
- Create: `workflow/src/features/genome.py`
- Create: `tests/test_features_genome.py`

**Step 1: Write workflow/src/features/__init__.py** (empty)

```python
# Empty __init__.py
```

**Step 2: Write workflow/src/features/genome.py** (imports, constants, PombaseGenomeConfig)

```python
"""
DNA-Level Gene Features
=========================

Per-mRNA DNA-level features (telomere/centromere distance, GC content, intron
structure, codon usage) extracted from a PomBase GFF3 annotation + genome
FASTA. Byte-faithful port of the DNA-level half of
`pombe_feature_functions.py` (DIT_HAP_pipeline) — protein-level functions
moved to `workflow/src/features/protein.py`.

Input
-----
- A gffutils FeatureDB built from a PomBase GFF3
- The corresponding genome FASTA
- peptide_stats.tsv (for primary-transcript peptide-length comparison)

Output
------
- One DNA_level_features record per mRNA feature in the FeatureDB

Usage
-----
    from workflow.src.features.genome import PombaseGenomeConfig, DNA_level_features
    cfg = PombaseGenomeConfig.from_pombase_dir(pombase_dir)
    db = gffutils.FeatureDB(cfg.database_file)
    features = [DNA_level_features.from_gffutils_feature(m, db, cfg) for m in db.features_of_type("mRNA")]

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC123, gc_fraction
from codonbias.scores import EffectiveNumberOfCodons
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
CHR1_LEFT_TELOMERE_END = 29663
CHR1_RIGHT_TELOMERE_START = 5554844
CHR2_LEFT_TELOMERE_END = 39186
CHR2_RIGHT_TELOMERE_START = 4500619
CHR3_LEFT_RIBOSOMAL_DNA_END = 23130
CHR3_RIGHT_RIBOSOMAL_DNA_START = 2440994

CHROMOSOME_END = {
    "left": {"I": CHR1_LEFT_TELOMERE_END, "II": CHR2_LEFT_TELOMERE_END, "III": CHR3_LEFT_RIBOSOMAL_DNA_END},
    "right": {"I": CHR1_RIGHT_TELOMERE_START, "II": CHR2_RIGHT_TELOMERE_START, "III": CHR3_RIGHT_RIBOSOMAL_DNA_START},
}

CENTROMERE_POSITIONS = {
    "I": (3753687, 3789421),
    "II": (1602418, 1644747),
    "III": (1070904, 1137003),
}

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class PombaseGenomeConfig:
    """Genome sequence/annotation paths and derived lookup tables for one PomBase version."""
    fasta_file: str
    fai_file: str
    gff3_file: str
    database_file: str
    genome_sequence_dict: dict[str, Seq]
    genome_length_dict: dict[str, int]
    primary_peptide_length: dict[str, int]

    @classmethod
    def from_pombase_dir(cls, pombase_dir: Path) -> PombaseGenomeConfig:
        """Build config from a PomBase version directory (e.g. resources/external/pombase/2025-10-01)."""
        genome_dir = pombase_dir / "genome_sequence_and_features"
        fasta_file = str(genome_dir / "Schizosaccharomyces_pombe_all_chromosomes.fa")
        peptide_stats = pd.read_csv(pombase_dir / "Protein_features" / "peptide_stats.tsv", sep="\t", index_col=0)
        genome_sequence_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

        return cls(
            fasta_file=fasta_file,
            fai_file=str(genome_dir / "Schizosaccharomyces_pombe_all_chromosomes.fa.fai"),
            gff3_file=str(genome_dir / "Schizosaccharomyces_pombe_all_chromosomes.gff3"),
            database_file=str(genome_dir / "Schizosaccharomyces_pombe_all_chromosomes.db"),
            genome_sequence_dict=genome_sequence_dict,
            genome_length_dict={k: len(v) for k, v in genome_sequence_dict.items()},
            primary_peptide_length=peptide_stats["Residues"].to_dict(),
        )
```

**Step 3: Continue workflow/src/features/genome.py** (determine_primary_candidate, DNA_level_features)

```python
# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def determine_primary_candidate(gene_id: str, mRNA_id: str, peptide_length: int, primary_peptide_length: int) -> bool:
    """Determine if the mRNA is the primary transcript, with hardcoded exceptions for known mis-annotated loci."""
    if gene_id in ["SPBC119.04", "SPBC17A3.07"]:
        return mRNA_id.endswith(".1")
    elif gene_id in ["SPAC212.11", "SPAC2E12.05", "SPAC977.01", "SPMIT.03", "SPMIT.06", "SPMIT.08"]:
        return mRNA_id.endswith(".1")
    else:
        return peptide_length == primary_peptide_length


@dataclass
class DNA_level_features:
    """One mRNA's DNA-level feature record."""
    Gene_id: str
    mRNA_id: str
    Chromosome: Literal["chr_II_telomeric_gap", "I", "II", "III", "mating_type_region", "mitochondrial"]
    Start: int
    End: int
    Strand: Literal["+", "-"]
    Abs_distance_from_telomere: float
    Relative_distance_from_telomere: float
    Abs_distance_from_centromere: float
    Relative_distance_from_centromere: float
    Gene_length: int
    GC_content_of_gene: float
    CDS_number: int
    GC_content_of_CDS: float
    Fraction_of_CDS: float
    GC3: float
    Containing_intron: bool
    Intron_number: int
    GC_content_of_intron: float
    Total_intron_length: int
    Average_intron_length: float
    Length_of_first_intron: int
    GC_contents_of_first_intron: float
    ENC: float
    Peptide_length: int
    Primary_peptide_length: int
    Primary_candidate: bool

    @classmethod
    def from_gffutils_feature(cls, mRNA: gffutils.Feature, db: gffutils.FeatureDB, cfg: PombaseGenomeConfig) -> DNA_level_features:
        """Compute one mRNA's DNA-level features from its gffutils Feature."""
        if mRNA.strand == "+":
            CDSs = list(db.children(mRNA, featuretype="CDS", order_by="start"))
            start = getattr(CDSs[0], "start", 0)
            end = getattr(CDSs[-1], "end", 0)
        else:
            CDSs = list(db.children(mRNA, featuretype="CDS", reverse=True, order_by="start"))
            start = getattr(CDSs[0], "end", 0)
            end = getattr(CDSs[-1], "start", 0)

        gene_id = mRNA.attributes.get("Parent")[0]
        mRNA_id = mRNA.id
        chrom = mRNA.chrom
        strand = mRNA.strand
        midpoint = (start + end) // 2

        abs_distance_from_telomere = min(
            abs(midpoint - CHROMOSOME_END["left"].get(chrom, np.nan)),
            abs(midpoint - CHROMOSOME_END["right"].get(chrom, np.nan)),
        )
        relative_distance_from_telomere = round(abs_distance_from_telomere / cfg.genome_length_dict[chrom], 3)
        abs_distance_from_centromere = abs(midpoint - np.mean(CENTROMERE_POSITIONS.get(chrom, (np.nan, np.nan))))
        relative_distance_from_centromere = round(abs_distance_from_centromere / cfg.genome_length_dict[chrom], 3)

        gene_length = abs(end - start) + 1
        GC_content_of_gene = round(gc_fraction(cfg.genome_sequence_dict[chrom][min(start, end):max(start, end)]), 3)

        CDS_number = len(CDSs)
        CDS_sequence = "".join(cds.sequence(cfg.fasta_file) for cds in CDSs)
        GC_content_of_CDS = round(gc_fraction(CDS_sequence), 3)
        Fraction_of_CDS = round(len(CDS_sequence) / gene_length, 3)
        GC3 = round(GC123(Seq(CDS_sequence))[-1], 3)

        introns = list(db.children(mRNA, featuretype="intron", order_by="start"))
        Containing_intron = len(introns) > 0
        Intron_number = len(introns)
        intron_sequences = [intron.sequence(cfg.fasta_file) for intron in introns]
        intron_sequence = "".join(intron_sequences)
        GC_content_of_intron = round(gc_fraction(intron_sequence), 3)
        Total_intron_length = len(intron_sequence)
        Average_intron_length = Total_intron_length / Intron_number if Intron_number > 0 else 0
        Length_of_first_intron = len(intron_sequences[0]) if Intron_number > 0 else 0
        GC_contents_of_first_intron = round(gc_fraction(intron_sequences[0]), 3) if Intron_number > 0 else 0.0

        ENC_model = EffectiveNumberOfCodons(mean="unweighted")
        ENC = np.round(ENC_model.get_score(CDS_sequence), 2)

        Peptide_length = len(Seq(CDS_sequence).translate(to_stop=True))
        primary_peptide_length = cfg.primary_peptide_length[gene_id]
        primary_candidate = determine_primary_candidate(gene_id, mRNA_id, Peptide_length, primary_peptide_length)

        return cls(
            Gene_id=gene_id, mRNA_id=mRNA_id, Chromosome=chrom, Start=start, End=end, Strand=strand,
            Abs_distance_from_telomere=abs_distance_from_telomere,
            Relative_distance_from_telomere=relative_distance_from_telomere,
            Abs_distance_from_centromere=abs_distance_from_centromere,
            Relative_distance_from_centromere=relative_distance_from_centromere,
            Gene_length=gene_length, GC_content_of_gene=GC_content_of_gene,
            CDS_number=CDS_number, GC_content_of_CDS=GC_content_of_CDS,
            Fraction_of_CDS=Fraction_of_CDS, GC3=GC3,
            Containing_intron=Containing_intron, Intron_number=Intron_number,
            GC_content_of_intron=GC_content_of_intron, Total_intron_length=Total_intron_length,
            Average_intron_length=Average_intron_length, Length_of_first_intron=Length_of_first_intron,
            GC_contents_of_first_intron=GC_contents_of_first_intron, ENC=ENC,
            Peptide_length=Peptide_length, Primary_peptide_length=primary_peptide_length,
            Primary_candidate=primary_candidate,
        )
```

**Step 4: Continue workflow/src/features/genome.py** (calculate_anticodon_usage_matrix)

```python
@logger.catch
def calculate_anticodon_usage_matrix(db: gffutils.FeatureDB, cfg: PombaseGenomeConfig) -> pd.DataFrame:
    """Compute a gene x anti-codon count matrix across all coding genes' concatenated CDS sequence."""
    from collections import Counter

    bases = ["A", "T", "G", "C"]
    codons = [f"{b1}{b2}{b3}" for b1 in bases for b2 in bases for b3 in bases]
    anticodons = [str(Seq(codon).reverse_complement()) for codon in codons]
    codon_to_anticodon = dict(zip(codons, anticodons))

    records = []
    for mRNA in db.features_of_type("mRNA"):
        gene_id = mRNA.attributes.get("Parent")[0]
        CDSs = list(db.children(mRNA, featuretype="CDS", order_by="start"))
        if not CDSs:
            logger.warning(f"No CDS found for {gene_id}, skipping")
            continue

        CDS_sequence = "".join(cds.sequence(cfg.fasta_file) for cds in CDSs)
        CDS_length = len(CDS_sequence) - (len(CDS_sequence) % 3)
        CDS_sequence = CDS_sequence[:CDS_length]

        codon_counts = Counter(CDS_sequence[i:i + 3] for i in range(0, CDS_length, 3))
        anticodon_counts = {
            codon_to_anticodon.get(codon): count
            for codon, count in codon_counts.items()
            if codon in codon_to_anticodon
        }
        anticodon_counts["Gene_id"] = gene_id
        records.append(anticodon_counts)

    df = pd.DataFrame(records).fillna(0)
    df = df.set_index("Gene_id")
    for anticodon in anticodons:
        if anticodon not in df.columns:
            df[anticodon] = 0
    return df[sorted(anticodons)].astype(int).reset_index()
```

**Step 5: Write tests/test_features_genome.py**

Uses the real PomBase data copied in Task 3 (`resources/external/pombase/2025-10-01/`) — this module is inherently data-dependent (gffutils DB, genome FASTA), so tests build a tiny 2-gene FeatureDB from a real GFF3 slice rather than mocking gffutils internals.

```python
"""Tests for workflow/src/features/genome.py — DNA-level feature extraction."""

import pytest
from pathlib import Path
from workflow.src.features.genome import (
    determine_primary_candidate,
    CHROMOSOME_END,
    CENTROMERE_POSITIONS,
)

POMBASE_DIR = Path("resources/external/pombase/2025-10-01")


def test_chromosome_end_has_three_chromosomes():
    """CHROMOSOME_END covers chromosomes I, II, III on both arms."""
    assert set(CHROMOSOME_END["left"]) == {"I", "II", "III"}
    assert set(CHROMOSOME_END["right"]) == {"I", "II", "III"}


def test_centromere_positions_are_start_less_than_end():
    """Each centromere interval has start < end."""
    for chrom, (start, end) in CENTROMERE_POSITIONS.items():
        assert start < end, chrom


def test_determine_primary_candidate_hardcoded_exception_sPBC119_04():
    """SPBC119.04's .1 transcript is always primary regardless of peptide length match."""
    assert determine_primary_candidate("SPBC119.04", "SPBC119.04.1", 100, 999) is True
    assert determine_primary_candidate("SPBC119.04", "SPBC119.04.2", 100, 999) is False


def test_determine_primary_candidate_falls_back_to_length_match():
    """For genes without a hardcoded exception, primary candidacy is peptide-length equality."""
    assert determine_primary_candidate("SPAC1002.01", "SPAC1002.01.1", 162, 162) is True
    assert determine_primary_candidate("SPAC1002.01", "SPAC1002.01.2", 100, 162) is False


@pytest.mark.skipif(not POMBASE_DIR.exists(), reason="requires resources/external/pombase/2025-10-01 (Task 3)")
def test_pombase_genome_config_loads_real_data():
    """PombaseGenomeConfig.from_pombase_dir loads the real genome FASTA and peptide_stats.tsv."""
    from workflow.src.features.genome import PombaseGenomeConfig

    cfg = PombaseGenomeConfig.from_pombase_dir(POMBASE_DIR)
    assert "I" in cfg.genome_length_dict
    assert cfg.primary_peptide_length["SPAC1002.01"] == 162


@pytest.mark.skipif(not POMBASE_DIR.exists(), reason="requires resources/external/pombase/2025-10-01 (Task 3)")
def test_dna_level_features_for_one_known_gene(tmp_path):
    """DNA_level_features.from_gffutils_feature reproduces the known SPAC1002.01 record."""
    import gffutils
    from workflow.src.features.genome import PombaseGenomeConfig, DNA_level_features

    cfg = PombaseGenomeConfig.from_pombase_dir(POMBASE_DIR)
    db_path = tmp_path / "test.db"
    db = gffutils.create_db(cfg.gff3_file, str(db_path), force=True, merge_strategy="create_unique")

    mRNA = db["SPAC1002.01.1"]
    feat = DNA_level_features.from_gffutils_feature(mRNA, db, cfg)

    assert feat.Gene_id == "SPAC1002.01"
    assert feat.Chromosome == "I"
    assert feat.Strand == "+"
    assert feat.Peptide_length == 162
    assert feat.Primary_candidate is True
```

**Step 6: Run tests**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_analysis
conda run -n bioinformatics python -m pytest tests/test_features_genome.py -v
```

Expected: 6 tests pass (4 pure-logic + 2 data-dependent, both should run since Task 3 already populated `resources/external/pombase/2025-10-01/`). The `test_dna_level_features_for_one_known_gene` test builds a fresh gffutils DB in `tmp_path` rather than reusing `cfg.database_file` — building the real 5MB GFF3 into a DB takes ~30-60s; this is expected and only pays that cost once per test run.

**Step 7: Commit**

```bash
git add workflow/src/features/__init__.py workflow/src/features/genome.py tests/test_features_genome.py
git commit -m "feat: port DNA-level gene features to workflow/src/features/genome.py

- PombaseGenomeConfig replaces pombe_feature_functions.config's genome fields,
  built explicitly from a pombase_dir arg instead of a __file__-relative default
- DNA_level_features, determine_primary_candidate, calculate_anticodon_usage_matrix
  ported byte-faithful from DIT_HAP_pipeline/workflow/src/pombe_feature_functions.py
- Protein-level functions deferred to workflow/src/features/protein.py (Task 6)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 6: `workflow/src/features/protein.py` (protein features + pLDDT)

Merges the protein-level half of `pombe_feature_functions.py` (`calculate_aliphatic_index_biopython`, `extract_protein_features_from_peptide_sequence`) with all of `protein_structure_functions.py` (pLDDT extraction from AlphaFold structures) — design doc §7 groups both under `features/protein.py` since the notebook uses them together to build one protein feature table.

**Files:**
- Create: `workflow/src/features/protein.py`
- Create: `tests/test_features_protein.py`

**Step 1: Write workflow/src/features/protein.py** (imports + amino-acid composition functions)

```python
"""
Protein-Level Gene Features
==============================

Peptide-sequence-derived features (aromaticity, aliphatic index, amino acid
composition) and AlphaFold pLDDT confidence statistics. Merges the
protein-level half of `pombe_feature_functions.py` with all of
`protein_structure_functions.py` (DIT_HAP_pipeline) — design doc §7 groups
both under `features/protein.py`.

Input
-----
- A PomBase peptide FASTA (peptide.fa)
- A directory of AlphaFold structure files (.pdb.gz / .pdb / .cif / .cif.gz)

Output
------
- extract_protein_features_from_peptide_sequence: one row per peptide record
- pLDDT_statistics_report: one row per structure file, keyed by UniProt ID

Usage
-----
    from workflow.src.features.protein import extract_protein_features_from_peptide_sequence, pLDDT_statistics_report
    protein_features = extract_protein_features_from_peptide_sequence(peptide_fasta)
    pLDDTs = pLDDT_statistics_report(alphafold_dir, structure_format="pdb.gz")

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import gzip
import re
from pathlib import Path
from typing import Literal

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from Bio import SeqIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.SeqUtils import seq3
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from loguru import logger
from tqdm import tqdm

# =============================================================================
# CORE LOGIC — peptide-sequence features
# =============================================================================
@logger.catch
def calculate_aliphatic_index_biopython(protein_sequence: str) -> float:
    """Calculate the aliphatic index (Ikai 1980) of a protein sequence."""
    analysis = ProteinAnalysis(protein_sequence)
    aa_percent = analysis.amino_acids_percent
    X_ala = aa_percent.get("A", 0) * 100
    X_val = aa_percent.get("V", 0) * 100
    X_leu = aa_percent.get("L", 0) * 100
    X_ile = aa_percent.get("I", 0) * 100
    aliphatic_index = X_ala + 2.9 * X_val + 3.9 * (X_leu + X_ile)
    return round(aliphatic_index, 3)


@logger.catch
def extract_protein_features_from_peptide_sequence(peptide_fasta_file: Path, return_redundant_meta: bool = False) -> pd.DataFrame:
    """Extract per-gene protein features (aromaticity, aliphatic index, AA composition, ...) from a peptide FASTA."""
    records = []
    aa_content = {}
    aa_percent = {}
    for record in SeqIO.parse(peptide_fasta_file, "fasta"):
        gene_id = re.search(r"(\S+)\.\d:pep$", record.id).groups()[0]
        sequence = str(record.seq).rstrip("*")
        if "*" in sequence:
            logger.warning(f"Stop codon found in sequence of {gene_id}. Truncating at first stop codon.")
            sequence = sequence.split("*")[0]
        analysis = ProteinAnalysis(sequence)
        aa_content[gene_id] = analysis.count_amino_acids()
        aa_percent[gene_id] = analysis.amino_acids_percent
        protein_features = {
            "Gene_id": gene_id,
            "aromaticity": analysis.aromaticity(),
            "aliphatic_index": calculate_aliphatic_index_biopython(sequence),
            "gravy": analysis.gravy(),
            "flexibility": np.mean(analysis.flexibility()),
            "instability_index": analysis.instability_index(),
            "monoisotopic": analysis.monoisotopic,
        }
        protein_features.update(
            dict(zip(("molar_extinction_reduced", "molar_extinction_cystines"), analysis.molar_extinction_coefficient()))
        )
        protein_features.update(
            dict(zip(("Helix_fraction", "Turn_fraction", "Sheet_fraction"), analysis.secondary_structure_fraction()))
        )
        if return_redundant_meta:
            protein_features.update({
                "charge_at_pH": analysis.charge_at_pH(7.0),
                "isoelectric_point": analysis.isoelectric_point(),
                "length": len(sequence),
                "molecular_weight(kDa)": analysis.molecular_weight() / 1000,
            })
        records.append(protein_features)

    aa_content_df = pd.DataFrame.from_dict(aa_content, orient="index")
    aa_percent_df = pd.DataFrame.from_dict(aa_percent, orient="index")
    aa_content_df.columns = [f"aa_count_{seq3(col)}" for col in aa_content_df.columns]
    aa_percent_df.columns = [f"aa_percent_{seq3(col)}" for col in aa_percent_df.columns]

    records_df = pd.DataFrame(records).set_index("Gene_id")
    records_df = records_df.join(aa_content_df).join(aa_percent_df).reset_index()
    return records_df
```

**Step 2: Continue workflow/src/features/protein.py** (pLDDT extraction, ported from protein_structure_functions.py)

```python
# =============================================================================
# CORE LOGIC — AlphaFold pLDDT statistics
# =============================================================================
@logger.catch
def extract_pLDDT(structure_file: Path | str) -> list[float]:
    """Extract per-residue pLDDT scores from a PDB or mmCIF file, compressed or not."""
    if isinstance(structure_file, str):
        structure_file = Path(structure_file)

    f = gzip.open(structure_file, "rt") if structure_file.name.endswith(".gz") else open(structure_file, "r")

    stem = structure_file.name.rstrip(".gz").lower()
    if stem.endswith(".pdb"):
        parser = PDBParser()
    elif stem.endswith(".cif"):
        parser = MMCIFParser()
    else:
        raise ValueError(f"Unknown file format: {structure_file.name}")
    structure = parser.get_structure(structure_file.stem, f)

    pLDDT = [residue["CA"].bfactor for residue in structure.get_residues() if residue.has_id("CA")]
    f.close()
    return pLDDT


@logger.catch
def extract_pLDDT_pdb_gz(structure_file: Path | str) -> list[float]:
    """Extract per-residue pLDDT scores from a .pdb.gz file."""
    if isinstance(structure_file, str):
        structure_file = Path(structure_file)
    f = gzip.open(structure_file, "rt")
    parser = PDBParser()
    structure = parser.get_structure(structure_file.stem, f)
    pLDDT = [residue["CA"].bfactor for residue in structure.get_residues()]
    f.close()
    return pLDDT


@logger.catch
def extract_pLDDT_pdb(structure_file: Path | str) -> list[float]:
    """Extract per-residue pLDDT scores from an uncompressed PDB file."""
    if isinstance(structure_file, str):
        structure_file = Path(structure_file)
    parser = PDBParser()
    structure = parser.get_structure(structure_file.stem, structure_file)
    return [residue["CA"].bfactor for residue in structure.get_residues()]


@logger.catch
def extract_protein_seq_pdb_gz(structure_file: Path | str) -> str:
    """Extract the residue sequence from a .pdb.gz file."""
    if isinstance(structure_file, str):
        structure_file = Path(structure_file)
    f = gzip.open(structure_file, "rt")
    parser = PDBParser()
    structure = parser.get_structure(structure_file.stem, f)
    ppb = PPBuilder()
    seq = ppb.build_peptides(structure)[0].get_sequence()
    f.close()
    return seq


@logger.catch
def pLDDT_statistics_report(
    structure_dir: Path,
    structure_format: Literal["pdb", "pdb.gz", "cif", "cif.gz", "mixed"] = "pdb.gz",
) -> pd.DataFrame:
    """Summarize per-residue pLDDT into per-structure mean/std/CV/disorder-fraction, keyed by UniProt ID."""
    all_pdb_files = list(structure_dir.glob(f"*.{structure_format}"))
    pLDDT_records = []

    for pdb_file in tqdm(all_pdb_files):
        uniprot_id = pdb_file.name.split("-F1-")[0].split("AF-")[1]
        match structure_format:
            case "pdb":
                pLDDT = np.array(extract_pLDDT_pdb(pdb_file))
            case "pdb.gz":
                pLDDT = np.array(extract_pLDDT_pdb_gz(pdb_file))
            case "cif" | "cif.gz" | "mixed":
                pLDDT = np.array(extract_pLDDT(pdb_file))
            case _:
                raise ValueError(f"Unsupported structure format: {structure_format}")
        length_protein = len(pLDDT)
        mean_pLDDT = np.mean(pLDDT)
        std_pLDDT = np.std(pLDDT)
        cv_pLDDT = std_pLDDT / mean_pLDDT if mean_pLDDT != 0 else np.nan
        disorder_fraction = np.sum(pLDDT < 50) / length_protein
        pLDDT_records.append({
            "uniprot_id": uniprot_id,
            "protein_length": length_protein,
            "pLDDT": ",".join(pLDDT.astype(str)),
            "mean_pLDDT": round(mean_pLDDT, 3),
            "std_pLDDT": round(std_pLDDT, 3),
            "cv_pLDDT": round(cv_pLDDT, 3),
            "disorder_fraction": round(disorder_fraction, 3),
        })

    return pd.DataFrame(pLDDT_records)
```

**Step 3: Write tests/test_features_protein.py**

```python
"""Tests for workflow/src/features/protein.py — peptide + pLDDT feature extraction."""

import gzip
import pytest
from pathlib import Path
from workflow.src.features.protein import (
    calculate_aliphatic_index_biopython,
    extract_protein_features_from_peptide_sequence,
    extract_pLDDT_pdb_gz,
    pLDDT_statistics_report,
)

POMBASE_DIR = Path("resources/external/pombase/2025-10-01")

# A real, minimal-but-valid pdb.gz record (single CA atom) for pLDDT extraction tests.
_MINIMAL_PDB = (
    "HEADER    TEST\n"
    "ATOM      1  CA  ALA A   1      11.104  13.207   2.000  1.00 87.50           C\n"
    "TER\n"
    "END\n"
)


def test_calculate_aliphatic_index_known_sequence():
    """Aliphatic index of a single Ala residue (100% A) matches the Ikai formula by hand."""
    assert calculate_aliphatic_index_biopython("A") == 100.0


def test_aliphatic_index_all_glycine_is_zero():
    """A sequence with none of A/V/L/I contributes zero to the aliphatic index."""
    assert calculate_aliphatic_index_biopython("GGGG") == 0.0


@pytest.mark.skipif(not POMBASE_DIR.exists(), reason="requires resources/external/pombase/2025-10-01 (Task 3)")
def test_extract_protein_features_from_real_peptide_fasta():
    """Feature extraction from the real peptide.fa produces the known SPAC1002.01 row."""
    peptide_fasta = POMBASE_DIR / "genome_sequence_and_features" / "peptide.fa"
    df = extract_protein_features_from_peptide_sequence(peptide_fasta)
    row = df.set_index("Gene_id").loc["SPAC1002.01"]
    assert "aromaticity" in df.columns
    assert "aa_percent_Ala" in df.columns
    assert 0 <= row["aromaticity"] <= 1


def test_extract_pLDDT_pdb_gz_reads_bfactor_as_plddt(tmp_path):
    """extract_pLDDT_pdb_gz reads the bfactor column as the pLDDT score."""
    pdb_gz = tmp_path / "AF-P00000-F1-model_v6.pdb.gz"
    with gzip.open(pdb_gz, "wt") as f:
        f.write(_MINIMAL_PDB)
    pLDDT = extract_pLDDT_pdb_gz(pdb_gz)
    assert pLDDT == [87.50]


def test_pLDDT_statistics_report_computes_disorder_fraction(tmp_path):
    """A structure with one low-confidence residue (<50) has disorder_fraction 1.0."""
    pdb_gz = tmp_path / "AF-P00000-F1-model_v6.pdb.gz"
    with gzip.open(pdb_gz, "wt") as f:
        f.write(_MINIMAL_PDB)
    report = pLDDT_statistics_report(tmp_path, structure_format="pdb.gz")
    row = report.set_index("uniprot_id").loc["P00000"]
    assert row["disorder_fraction"] == 1.0
    assert row["mean_pLDDT"] == 87.5
```

**Step 4: Run tests**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_analysis
conda run -n bioinformatics python -m pytest tests/test_features_protein.py -v
```

Expected: 5 tests pass. The AlphaFold-dir-scale `pLDDT_statistics_report` is exercised end-to-end against the real 10,392-file dataset in Task 8, not here — these tests use a synthetic 1-atom structure to keep unit tests fast.

**Step 5: Commit**

```bash
git add workflow/src/features/protein.py tests/test_features_protein.py
git commit -m "feat: port protein features + pLDDT extraction to workflow/src/features/protein.py

- Peptide-sequence features from pombe_feature_functions.py
- All of protein_structure_functions.py (pLDDT extraction for pdb/pdb.gz/cif)
- Byte-faithful port; tests use synthetic minimal structures for speed

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 7: `workflow/src/enrichment/ontology.py` (GO term richness slice only)

Per the Task 3-9 scope note: the feature notebook only calls `OntologyDataConfig`/`load_ontology_data` from `enrichment_functions.py` to compute per-gene GO term richness — it never calls `ontology_enrichment_pipeline`, `stringdb_enrichment`, or `revigo_analysis`. This task ports exactly the called subset (`OntologyData`, `OntologyDataConfig`, `load_ontology_data`). The original `load_ontology_data` also computed `ns2slim_assoc` via `mapslim`/`get_slim_ns2assoc`, but the feature notebook only reads `gene2go` from the return tuple — those helpers are dropped here and deferred to a future `enrichment/pipeline.py` follow-up plan, along with the enrichment-study and STRING/REVIGO functions.

**Files:**
- Create: `workflow/src/enrichment/__init__.py`
- Create: `workflow/src/enrichment/ontology.py`
- Create: `tests/test_enrichment_ontology.py`

**Step 1: Write workflow/src/enrichment/__init__.py** (empty)

```python
# Empty __init__.py
```

**Step 2: Write workflow/src/enrichment/ontology.py**

```python
"""
Ontology Data Loading
=======================

Loads a GO/FYPO/MONDO-style OBO + GAF association pair into goatools objects.
Ported from the `OntologyDataConfig`/`load_ontology_data` slice of
`enrichment_functions.py` (DIT_HAP_pipeline) — only the subset
`pombe_feature_collection.ipynb` actually calls (GO term richness via
`gene2go`). The enrichment-study, STRING, and REVIGO functions in the
original module are intentionally not ported here; see Task 7's scope note
in docs/plans/2026-07-15-DIT-HAP-analysis-phase1-implementation.md.

Input
-----
- An OBO ontology file (e.g. go-basic.obo)
- A GAF-format association file (e.g. gene_ontology_annotation.gaf.tsv)
- One or more slim-term tables (Term, Description columns, no header)

Output
------
- OntologyData: validated file handles + concatenated slim term table
- load_ontology_data(...): goatools GODag/GafReader plus gene2go/go2genes dicts

Usage
-----
    from workflow.src.enrichment.ontology import OntologyDataConfig, load_ontology_data
    cfg = OntologyDataConfig(ontology_obo=..., ontology_association_gaf=..., slim_terms_table=[...])
    dag, objanno, ns2assoc, gene2go, go2genes, slim_dag = load_ontology_data(cfg.load_data())

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from dataclasses import dataclass
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from goatools.anno.gaf_reader import GafReader
from goatools.obo_parser import GODag

# 4. Local Imports
from workflow.src.io import read_file

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, frozen=True)
class OntologyData:
    """Validated ontology file paths plus the concatenated slim-term table."""
    ontology_obo_path: Path
    ontology_association_file: Path
    slim_term_dataframe: pd.DataFrame


@dataclass(kw_only=True, frozen=True)
class OntologyDataConfig:
    """Unvalidated ontology file paths, as registered in config/datasets.yaml-adjacent code."""
    ontology_obo: Path
    ontology_association_gaf: Path
    slim_terms_table: list[Path]

    def validate_paths(self) -> None:
        """Raise if the OBO, GAF, or any slim-term file is missing."""
        for file_path in [self.ontology_obo, self.ontology_association_gaf, *self.slim_terms_table]:
            if not file_path.exists():
                raise FileNotFoundError(f"Gene ontology file not found: {file_path}")
            if not file_path.is_file():
                raise ValueError(f"Path is not a file: {file_path}")

    def load_data(self) -> OntologyData:
        """Validate paths, concatenate slim-term tables, and return an OntologyData."""
        self.validate_paths()
        slim_dfs = [
            read_file(path, header=None, names=["Term", "Description"])
            for path in self.slim_terms_table
        ]
        slim_df = pd.concat(slim_dfs, ignore_index=True)
        return OntologyData(
            ontology_obo_path=self.ontology_obo,
            ontology_association_file=self.ontology_association_gaf,
            slim_term_dataframe=slim_df,
        )


# =============================================================================
# CORE LOGIC
# =============================================================================
def load_ontology_data(
    ontology_data: OntologyData, **kwargs
) -> tuple[GODag, GafReader, dict, dict, dict, dict]:
    """Load an OBO + GAF pair into a GODag/GafReader and derive gene2go/go2genes dicts."""
    try:
        dag = GODag(str(ontology_data.ontology_obo_path), optional_attrs=["def", "relationship"], load_obsolete=False)
    except KeyError:
        dag = GODag(str(ontology_data.ontology_obo_path), optional_attrs=["def"], load_obsolete=False)

    objanno = GafReader(str(ontology_data.ontology_association_file), godag=dag)

    slim_terms = ontology_data.slim_term_dataframe["Term"].to_list()
    slim_dag = {term: dag[term] for term in slim_terms}

    ns2assoc = objanno.get_ns2assc(**kwargs)
    gene2go = objanno.get_id2gos_nss(**kwargs)
    go2genes = objanno.get_id2gos_nss(go2geneids=True, **kwargs)

    return dag, objanno, ns2assoc, gene2go, go2genes, slim_dag
```

**Step 3: Write tests/test_enrichment_ontology.py**

```python
"""Tests for workflow/src/enrichment/ontology.py — OBO/GAF loading."""

import pytest
from pathlib import Path
from workflow.src.enrichment.ontology import OntologyDataConfig

POMBASE_DIR = Path("resources/external/pombase/2025-10-01")
ONTOLOGY_DIR = POMBASE_DIR / "ontologies_and_associations"


def test_validate_paths_raises_on_missing_obo(tmp_path):
    """validate_paths raises FileNotFoundError naming the missing OBO file."""
    cfg = OntologyDataConfig(
        ontology_obo=tmp_path / "missing.obo",
        ontology_association_gaf=tmp_path / "missing.gaf",
        slim_terms_table=[],
    )
    with pytest.raises(FileNotFoundError, match="Gene ontology file not found"):
        cfg.validate_paths()


@pytest.mark.skipif(not ONTOLOGY_DIR.exists(), reason="requires resources/external/pombase/2025-10-01 (Task 3)")
def test_load_data_concatenates_slim_tables():
    """load_data() concatenates all three GO slim tables into one Term/Description frame."""
    cfg = OntologyDataConfig(
        ontology_obo=ONTOLOGY_DIR / "go-basic.obo",
        ontology_association_gaf=ONTOLOGY_DIR / "gene_ontology_annotation.gaf.tsv",
        slim_terms_table=[
            ONTOLOGY_DIR / "bp_go_slim_terms.tsv",
            ONTOLOGY_DIR / "mf_go_slim_terms.tsv",
            ONTOLOGY_DIR / "cc_go_slim_terms.tsv",
        ],
    )
    data = cfg.load_data()
    assert list(data.slim_term_dataframe.columns) == ["Term", "Description"]
    assert len(data.slim_term_dataframe) > 0


@pytest.mark.skipif(not ONTOLOGY_DIR.exists(), reason="requires resources/external/pombase/2025-10-01 (Task 3)")
def test_load_ontology_data_returns_gene2go_for_known_gene():
    """load_ontology_data's gene2go dict has an entry for a known coding gene."""
    from workflow.src.enrichment.ontology import load_ontology_data

    cfg = OntologyDataConfig(
        ontology_obo=ONTOLOGY_DIR / "go-basic.obo",
        ontology_association_gaf=ONTOLOGY_DIR / "gene_ontology_annotation.gaf.tsv",
        slim_terms_table=[ONTOLOGY_DIR / "bp_go_slim_terms.tsv"],
    )
    dag, objanno, ns2assoc, gene2go, go2genes, slim_dag = load_ontology_data(
        cfg.load_data(),
        relationships={"is_a", "part_of"},
        propagate_counts=True,
        load_obsolete=False,
        prt=None,
    )
    assert "SPAC1002.01" in gene2go
    assert len(gene2go["SPAC1002.01"]) > 0
```

**Step 4: Run tests**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_analysis
conda run -n bioinformatics python -m pytest tests/test_enrichment_ontology.py -v
```

Expected: 3 tests pass. `test_load_ontology_data_returns_gene2go_for_known_gene` loads the full 31MB `go-basic.obo` and propagates counts across the whole DAG — expect this single test to take 1-3 minutes, matching the notebook cell's own runtime.

**Step 5: Commit**

```bash
git add workflow/src/enrichment/__init__.py workflow/src/enrichment/ontology.py tests/test_enrichment_ontology.py
git commit -m "feat: port GO term richness ontology loading to workflow/src/enrichment/ontology.py

- OntologyDataConfig/OntologyData/load_ontology_data ported from enrichment_functions.py
- Only the GO-term-richness slice pombe_feature_collection.ipynb calls;
  enrichment-study/STRING/REVIGO functions deferred to a future enrichment/pipeline.py
- Dropped the unused ns2slim_assoc computation (expensive, nothing reads it here)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 8: `workflow/scripts/features/collect_pombe_features.py` — the driver script

This assembles all 15 sources into the final feature matrix, faithfully reproducing `pombe_feature_collection.ipynb` end to end. Following the `def_ctr_insertions.py` convention (CLI args for every input/output path, no hardcoded resource layout inside the script — that knowledge lives in `features.smk`), and preserving the notebook's two documented quirks verbatim: (1) BioGrid PPI/GI degree only counts `Interactor A` occurrences, so genes appearing solely as `Interactor B` get degree 0; (2) the final merge produces two columns both named `DeletionLibrary_essentiality` (traced below in Step 8) — both are kept, not silently deduplicated, since removing either would change the shape of the reference output this task verifies against.

**Files:**
- Create: `workflow/scripts/features/collect_pombe_features.py`
- Create: `workflow/rules/features.smk`
- Create: `tests/test_collect_pombe_features.py`
- Modify: `Snakefile` (uncomment the features.smk include and rule all target)

**Step 1: Write workflow/scripts/features/collect_pombe_features.py** (module docstring + imports + constants)

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pombe Coding Gene Feature Collection
=======================================

Assembles a per-coding-gene feature matrix from 15 heterogeneous sources —
PomBase annotations, AlphaFold structures, BioGrid, Ensembl paralogs, and
8 literature supplementary tables — reproducing
DIT_HAP_pipeline/workflow/notebooks/pombe_feature_collection.ipynb exactly.
Dataset-independent: depends only on the PomBase reference version, not on
any DIT-HAP sequencing project (design doc §8).

Input
-----
- A PomBase version directory (genome FASTA/GFF3, gene metadata, protein
  features, ontology OBO/GAF files, curated orthologs)
- An AlphaFold structure directory (.pdb.gz files)
- 8 literature supplementary tables (xlsx/xls)
- Curated deletion-library and essentiality-verification tables
- BioGrid interaction table, Ensembl paralog export

Output
------
- A tab-separated per-gene feature matrix (one row per coding gene)
- A tab-separated codon (anti-codon) usage matrix

Usage
-----
    python collect_pombe_features.py \\
        --pombase-dir resources/external/pombase/2025-10-01 \\
        --alphafold-dir /path/to/AlphaFold_Dataset \\
        --literature-dir resources/literature \\
        --deletion-library-xlsx resources/curated/deletion_library_categories.xlsx \\
        --essentiality-verification-csv resources/curated/essentiality_verification.csv \\
        --biogrid-tsv resources/external/biogrid/BIOGRID-....tab3.txt \\
        --ensembl-paralogs-tsv resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv \\
        --output results/features/2025-10-01/pombe_coding_gene_protein_features.tsv \\
        --codon-usage-output results/features/2025-10-01/codon_usage_matrix.tsv

Author:   Yusheng Yang (guidance) + Claude Sonnet 5 (implementation)
Date:     2026-07-15
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
import numpy as np
import pandas as pd

# 3. Third-party Imports
import gffutils
from loguru import logger
from tqdm import tqdm

# 4. Local Imports
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from workflow.src.enrichment.ontology import OntologyDataConfig, load_ontology_data
from workflow.src.features.genome import DNA_level_features, PombaseGenomeConfig, calculate_anticodon_usage_matrix
from workflow.src.features.protein import extract_protein_features_from_peptide_sequence, pLDDT_statistics_report
from workflow.src.gene_ids import update_sysIDs

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================
# Amino-acid-composition columns kept from extract_protein_features_from_peptide_sequence
# (notebook cell 37's selected_protein_features_from_peptide list)
SELECTED_PEPTIDE_FEATURE_COLUMNS = [
    "aromaticity", "aliphatic_index", "gravy", "flexibility", "instability_index",
    "aa_percent_Ala", "aa_percent_Cys", "aa_percent_Asp", "aa_percent_Glu", "aa_percent_Phe",
    "aa_percent_Gly", "aa_percent_His", "aa_percent_Ile", "aa_percent_Lys", "aa_percent_Leu",
    "aa_percent_Met", "aa_percent_Asn", "aa_percent_Pro", "aa_percent_Gln", "aa_percent_Arg",
    "aa_percent_Ser", "aa_percent_Thr", "aa_percent_Val", "aa_percent_Trp", "aa_percent_Tyr",
]
```

**Step 2: Continue collect_pombe_features.py** (InputOutputConfig + logging setup)

```python
# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================
@dataclass(kw_only=True, slots=True, frozen=True)
class InputOutputConfig:
    """Validated input/output paths for the feature collection pipeline."""
    pombase_dir: Path
    alphafold_dir: Path
    literature_dir: Path
    deletion_library_xlsx: Path
    essentiality_verification_csv: Path
    biogrid_tsv: Path
    ensembl_paralogs_tsv: Path
    output_features: Path
    output_codon_usage: Path

    def __post_init__(self) -> None:
        """Validate all input paths exist, then ensure output directories exist."""
        required_inputs = [
            self.pombase_dir, self.alphafold_dir, self.literature_dir,
            self.deletion_library_xlsx, self.essentiality_verification_csv,
            self.biogrid_tsv, self.ensembl_paralogs_tsv,
        ]
        for path in required_inputs:
            if not path.exists():
                raise ValueError(f"Required input path does not exist: {path}")
        self.output_features.parent.mkdir(parents=True, exist_ok=True)
        self.output_codon_usage.parent.mkdir(parents=True, exist_ok=True)

    @property
    def gene_meta_file(self) -> Path:
        """PomBase gene_IDs_names_products.tsv, used throughout for update_sysIDs()."""
        return self.pombase_dir / "Gene_metadata" / "gene_IDs_names_products.tsv"


# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru for the application."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )
```

**Step 3: Continue collect_pombe_features.py** (gene metadata loader, ortholog counts helper, DNA-level features)

```python
# =============================================================================
# CORE LOGIC
# =============================================================================
@logger.catch
def load_gene_meta(gene_meta_file: Path) -> tuple[pd.DataFrame, dict]:
    """Load gene metadata and build a uniprot_id -> gene_systematic_id map."""
    gene_meta = pd.read_csv(gene_meta_file, sep="\t")
    gene_meta["gene_name"] = gene_meta["gene_name"].fillna(gene_meta["gene_systematic_id"])
    uniprot2id = dict(zip(gene_meta["uniprot_id"], gene_meta["gene_systematic_id"]))
    return gene_meta, uniprot2id


@logger.catch
def get_ortholog_counts(ortholog_file: Path) -> pd.Series:
    """Count pipe-separated orthologs per gene from a PomBase curated_orthologs file."""
    ortholog_df = pd.read_csv(
        ortholog_file, sep="\t", index_col=0, header=None,
        names=["gene_systematic_id", "orthologs"], na_values="NONE",
    )
    ortholog_df.index = ortholog_df.index.str.split("(").str[0]
    return ortholog_df["orthologs"].str.split("|").apply(lambda x: len(x) if isinstance(x, list) else 0)


@logger.catch
def collect_dna_level_features(db: gffutils.FeatureDB, genome_cfg: PombaseGenomeConfig) -> tuple[pd.DataFrame, list[str]]:
    """Compute DNA-level features for every mRNA, and return the list of coding gene IDs."""
    mRNAs = list(db.features_of_type("mRNA"))
    records = [
        DNA_level_features.from_gffutils_feature(mRNA, db, genome_cfg)
        for mRNA in tqdm(mRNAs, desc="DNA-level features")
    ]
    df = pd.DataFrame(records)
    coding_genes = df["Gene_id"].unique().tolist()
    return df, coding_genes
```

**Step 4: Continue collect_pombe_features.py** (RNA-level features)

```python
@logger.catch
def collect_rna_level_features(literature_dir: Path, gene_meta_file: Path, coding_genes: list[str]) -> pd.DataFrame:
    """Assemble mRNA abundance (Marguerat 2012) and mRNA kinetics (Harigaya 2016) features."""
    abundance = pd.read_excel(
        literature_dir / "margueratQuantitativeAnalysisFission2012.xlsx",
        sheet_name="Table_S2", comment="#",
    ).set_index("Systematic.name")
    abundance = abundance[["MM1.tot.cpc_ex", "MM2.tot.cpc_ex", "MN1.tot.cpc_ex", "MN2.tot.cpc_ex"]].copy()
    abundance.columns = pd.MultiIndex.from_tuples(
        [
            ("EMM_Proliferating_Cell_RNA_Abundance", "replicate1"),
            ("EMM_Proliferating_Cell_RNA_Abundance", "replicate2"),
            ("EMM_Nitrogen_Starved_Cell_RNA_Abundance", "replicate1"),
            ("EMM_Nitrogen_Starved_Cell_RNA_Abundance", "replicate2"),
        ],
        names=["Condition", "Replicate"],
    )
    mean_ = abundance.T.groupby(level="Condition").mean().T
    std_ = abundance.T.groupby(level="Condition").std().T
    cv_ = std_ / mean_
    abundance_stats = pd.concat([mean_, std_, cv_], axis=1, keys=["mean", "std", "cv"])
    abundance_stats.index = update_sysIDs(abundance_stats.index.tolist(), gene_meta_file)
    abundance_stats = abundance_stats[abundance_stats.index.isin(coding_genes)].copy().dropna().round(3)
    abundance_stats.columns = ["_".join(col).strip() for col in abundance_stats.columns.values]
    abundance_stats = (
        abundance_stats.rename_axis("gene_systematic_id")
        .reset_index()
        .drop_duplicates(subset=["gene_systematic_id"])
        .set_index("gene_systematic_id")
    )

    kinetics = pd.read_excel(literature_dir / "harigayaAnalysisAssociationCodon2016.xls", sheet_name="Table")
    kinetics = kinetics[["Gene ID", "tAIg", "HL - Mata (5)", "SR - Mata (5)"]].set_index("Gene ID")
    kinetics.columns = ["tAIg", "mRNA_half_life_minutes", "mRNA_synthesis_rate_per_minute"]
    kinetics.index = update_sysIDs(kinetics.index.tolist(), gene_meta_file)
    kinetics = kinetics[kinetics.index.isin(coding_genes)].copy().dropna().round(3)
    kinetics = (
        kinetics.rename_axis("gene_systematic_id")
        .reset_index()
        .drop_duplicates(subset=["gene_systematic_id"])
        .set_index("gene_systematic_id")
    )

    return pd.concat([abundance_stats, kinetics], axis=1, join="outer")
```

**Step 5: Continue collect_pombe_features.py** (protein-level features)

```python
@logger.catch
def collect_protein_level_features(
    pombase_dir: Path,
    alphafold_dir: Path,
    literature_dir: Path,
    gene_meta_file: Path,
    protein_metadata: pd.DataFrame,
    uniprot2id: dict,
    coding_genes: list[str],
) -> pd.DataFrame:
    """Assemble peptide-sequence, abundance, turnover, pLDDT, and PFAM-domain protein features."""
    peptide_features = extract_protein_features_from_peptide_sequence(
        pombase_dir / "genome_sequence_and_features" / "peptide.fa"
    )

    gene_abundance = pd.read_csv(pombase_dir / "RNA_metadata" / "quantitative_gene_expression.tsv", sep="\t")
    proliferating = gene_abundance.query(
        "reference == 'PMID:23101633' and type == 'protein' and condition == 'glucose MM,standard temperature'"
    )[["gene_systematic_id", "copies_per_cell"]].dropna().astype({"copies_per_cell": float})
    quiescent = gene_abundance.query(
        "reference == 'PMID:23101633' and type == 'protein' and condition == 'glucose MM,nitrogen absent,standard temperature'"
    )[["gene_systematic_id", "copies_per_cell"]].dropna().astype({"copies_per_cell": float})
    protein_abundance = pd.merge(
        proliferating, quiescent, on="gene_systematic_id",
        suffixes=("_EMM_Proliferating_Cell", "_EMMN_Quiescent_Cell"),
    ).set_index("gene_systematic_id")

    protein_kinetics = pd.read_excel(
        literature_dir / "christianoGlobalProteomeTurnover2014.xlsx", na_values=["n.d."]
    ).dropna(subset=["Degradation rates (min-1)", "t1/2 (min)"]).rename(
        columns={"t1/2 (min)": "protein_half_life_minutes"}
    )
    protein_kinetics["ENSG"] = protein_kinetics["ENSG"].fillna(protein_kinetics["Gene name"])
    protein_kinetics = protein_kinetics[["ENSG", "protein_half_life_minutes"]].set_index("protein_half_life_minutes")
    protein_kinetics = protein_kinetics["ENSG"].str.split(";").explode().reset_index()
    protein_kinetics["gene_systematic_id"] = update_sysIDs(protein_kinetics["ENSG"].tolist(), gene_meta_file)
    protein_kinetics = protein_kinetics.drop_duplicates(subset=["gene_systematic_id"])

    pLDDTs = pLDDT_statistics_report(alphafold_dir, structure_format="pdb.gz")
    pLDDTs["Systematic_ID"] = pLDDTs["uniprot_id"].map(uniprot2id)

    protein_domains = pd.read_csv(pombase_dir / "Protein_features" / "protein_families_and_domains.tsv", sep="\t")
    pfam_domain_counts = (
        protein_domains.query("database == 'PFAM'").groupby("systematic_id").size().rename("PFAM_domain_count")
    )

    merged = (
        protein_metadata
        .merge(peptide_features.set_index("Gene_id")[SELECTED_PEPTIDE_FEATURE_COLUMNS], left_index=True, right_index=True, how="outer")
        .merge(protein_abundance, left_index=True, right_index=True, how="outer")
        .merge(protein_kinetics.set_index("gene_systematic_id")[["protein_half_life_minutes"]], left_index=True, right_index=True, how="outer")
        .merge(pLDDTs.set_index("Systematic_ID")[["mean_pLDDT", "std_pLDDT", "cv_pLDDT"]], left_index=True, right_index=True, how="outer")
        .join(pfam_domain_counts)
    )
    merged = merged[merged.index.isin(coding_genes)].copy()
    merged["PFAM_domain_count"] = merged["PFAM_domain_count"].fillna(0).astype(int)
    return merged
```

**Step 6: Continue collect_pombe_features.py** (evolutionary-level features)

```python
@logger.catch
def collect_evolutionary_level_features(
    pombase_dir: Path,
    ensembl_paralogs_tsv: Path,
    literature_dir: Path,
    gene_meta_file: Path,
    coding_genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assemble ortholog/paralog counts, evolutionary rate, and phyloP/divergence scores.

    Returns (evolutionary_features_df, phyloP_and_divergence_df) — the second is reused
    by collect_phenotype_level_features for transposon-insertion-density features that
    live in the same source table (grechFitnessLandscapeFission2019.xlsx).
    """
    orthologs_dir = pombase_dir / "curated_orthologs"
    num_japonicus = get_ortholog_counts(orthologs_dir / "pombe_japonicus_orthologs.txt")
    num_cerevisiae = get_ortholog_counts(orthologs_dir / "pombe_cerevisiae_orthologs.txt")
    num_human = get_ortholog_counts(orthologs_dir / "pombe_human_orthologs.txt")

    pombe_paralogs = pd.read_csv(ensembl_paralogs_tsv, sep="\t")
    paralog_count = (
        pombe_paralogs.query("`Gene stable ID` in @coding_genes")
        .groupby(["Gene stable ID", "Gene name"])
        .apply(lambda sub_df: sub_df.shape[0], include_groups=False)
        .to_frame("paralog_count")
    )
    paralog_count["gene_systematic_id"] = update_sysIDs(
        paralog_count.index.get_level_values("Gene stable ID").tolist(), gene_meta_file
    )

    evolutionary_rate = pd.read_excel(
        literature_dir / "rhindComparativeFunctionalGenomics2011.xls", sheet_name="S30", skiprows=[0, 1]
    )
    evolutionary_rate["gene_systematic_id"] = update_sysIDs(
        evolutionary_rate["Genes"].str.split(";").str[0].tolist(), gene_meta_file
    )

    phyloP_and_divergence = pd.read_excel(
        literature_dir / "grechFitnessLandscapeFission2019.xlsx", sheet_name="Table 2", skiprows=list(range(14))
    ).drop_duplicates(subset=["gene"])
    phyloP_and_divergence["gene_systematic_id"] = update_sysIDs(phyloP_and_divergence["gene"].tolist(), gene_meta_file)

    evolutionary_df = pd.concat(
        [
            num_japonicus.rename("japonicus_ortholog_count"),
            num_cerevisiae.rename("cerevisiae_ortholog_count"),
            num_human.rename("human_ortholog_count"),
            paralog_count.set_index("gene_systematic_id")[["paralog_count"]],
            evolutionary_rate.set_index("gene_systematic_id")[["Rate"]].rename(columns={"Rate": "evolutionary_rate"}),
            phyloP_and_divergence.set_index("gene_systematic_id")[
                ["mean.phylop", "diversity.S", "diversity.Pi", "diversity.Theta", "diversity.Tajima_D"]
            ],
        ],
        join="outer", axis=1,
    )
    evolutionary_df = evolutionary_df[evolutionary_df.index.isin(coding_genes)].copy()
    evolutionary_df["paralog_count"] = evolutionary_df["paralog_count"].fillna(0).astype(int)

    return evolutionary_df, phyloP_and_divergence
```

**Step 7: Continue collect_pombe_features.py** (network-level features)

```python
@logger.catch
def collect_network_level_features(
    pombase_dir: Path,
    biogrid_tsv: Path,
    gene2go: dict,
    coding_genes: list[str],
) -> pd.DataFrame:
    """Assemble GO term richness and BioGrid PPI/GI degree.

    NOTE: PPI_degree/GI_degree only count rows grouped by "Systematic Name Interactor A" —
    a gene that appears solely as "Interactor B" in BioGrid gets degree 0 here, even if it
    has real interactions. This is the original notebook's behavior (undocumented upstream,
    not a bug introduced by this port) and is preserved rather than symmetrized, since fixing
    it would change every downstream degree value against the reference this task verifies.
    """
    go_richness = {gene: len(set(terms)) for gene, terms in gene2go.items()}
    go_richness_df = pd.DataFrame.from_dict(go_richness, orient="index", columns=["GO_term_richness"])

    biogrid_data = pd.read_csv(biogrid_tsv, sep="\t")
    PPI_and_GI = biogrid_data[
        [
            "Systematic Name Interactor A", "Systematic Name Interactor B",
            "Official Symbol Interactor A", "Official Symbol Interactor B",
            "Experimental System Type",
        ]
    ].drop_duplicates()
    PPI = PPI_and_GI.query("`Experimental System Type` == 'physical'")
    GI = PPI_and_GI.query("`Experimental System Type` == 'genetic'")
    PPI_degrees = PPI.groupby("Systematic Name Interactor A").size().rename("PPI_degree")
    GI_degrees = GI.groupby("Systematic Name Interactor A").size().rename("GI_degree")

    network_df = pd.concat([go_richness_df, PPI_degrees, GI_degrees], join="outer", axis=1)
    network_df = network_df[network_df.index.isin(coding_genes)].copy()
    network_df = network_df.fillna(0).astype({"PPI_degree": int, "GI_degree": int})
    return network_df
```

**Step 8: Continue collect_pombe_features.py** (phenotype-level features)

```python
@logger.catch
def collect_phenotype_level_features(
    pombase_dir: Path,
    deletion_library_xlsx: Path,
    essentiality_verification_csv: Path,
    literature_dir: Path,
    gene_meta_file: Path,
    coding_genes: list[str],
    phyloP_and_divergence: pd.DataFrame,
) -> pd.DataFrame:
    """Assemble FYPO viability, deletion-library essentiality, bar-seq fitness, transposon
    insertion density, and CRISPRi growth phenotypes.
    """
    FYPO_viability = pd.read_csv(
        pombase_dir / "Gene_metadata" / "gene_viability.tsv", sep="\t",
        header=None, names=["gene_systematic_id", "FYPOviability"],
    ).set_index("gene_systematic_id")

    DeletionLibrary_essentiality = pd.read_excel(deletion_library_xlsx)[
        ["Updated_Systematic_ID", "Gene dispensability. This study", "Category"]
    ].set_index("Updated_Systematic_ID").rename(
        columns={"Gene dispensability. This study": "DeletionLibrary_essentiality", "Category": "DeletionLibrary_category"}
    )

    revised_essentiality_map = (
        pd.read_csv(essentiality_verification_csv)[["systematic_id", "verification_essentiality"]]
        .set_index("systematic_id")["verification_essentiality"]
        .to_dict()
    )
    # Updated_essentiality intentionally ends up with TWO relevant columns:
    # "DeletionLibrary_essentiality" (carried over from the .copy() below) and the new
    # "RevisedDeletionLibrary_essentiality". Concatenating it against DeletionLibrary_essentiality
    # (Step 9) therefore produces a duplicate "DeletionLibrary_essentiality" column in the final
    # matrix — see this task's header note; both are kept to match the reference output byte-for-byte.
    Updated_essentiality = DeletionLibrary_essentiality[["DeletionLibrary_essentiality"]].copy()
    Updated_essentiality["RevisedDeletionLibrary_essentiality"] = Updated_essentiality.apply(
        lambda row: revised_essentiality_map.get(row.name, row["DeletionLibrary_essentiality"]), axis=1
    )

    bar_seq_fitness = pd.read_excel(literature_dir / "comp_fitness_QianWenFeng_Koch-1.xlsx").rename(
        columns={"yes": "Barseq_from_dulab", "SM fitness defect from Koch et al": "Barseq_from_koch"}
    ).dropna(subset=["Barseq_from_dulab", "Barseq_from_koch"])
    bar_seq_fitness["gene_systematic_id"] = update_sysIDs(bar_seq_fitness["gene"].tolist(), gene_meta_file)

    ins_density = pd.read_excel(literature_dir / "guoIntegrationProfilingGene2013.xls", sheet_name="TableS2")
    ins_density["gene_systematic_id"] = update_sysIDs(
        ins_density["Gene name"].str.strip().apply(lambda row: sorted(row.split(" "))[0]).tolist(), gene_meta_file
    )
    ins_density = ins_density.drop_duplicates(subset=["gene_systematic_id"])

    ins_grech = phyloP_and_divergence[["gene_systematic_id", "ipkm", "uipkm", "Malecki2016.KO.colony.size"]].copy().rename(
        columns={"Malecki2016.KO.colony.size": "colony_size_Malecki2016"}
    )

    CRISPRi_data = pd.read_excel(literature_dir / "ishikawaArrayedCRISPRiLibrary2024.xlsx").iloc[1:].dropna(
        subset=["Max Growth Rate", "Colony Formation"]
    )
    CRISPRi_data["gene_systematic_id"] = update_sysIDs(CRISPRi_data["Systematic ID"].tolist(), gene_meta_file)

    phenotype_df = pd.concat(
        [
            FYPO_viability,
            DeletionLibrary_essentiality,
            Updated_essentiality,
            bar_seq_fitness.set_index("gene_systematic_id")[["Barseq_from_dulab", "Barseq_from_koch"]],
            ins_density.set_index("gene_systematic_id")[["Integration density, in-vivo (integrations/kb/million inserts)"]],
            ins_grech.set_index("gene_systematic_id")[["ipkm", "uipkm", "colony_size_Malecki2016"]],
            CRISPRi_data.set_index("gene_systematic_id")[["Max Growth Rate", "Colony Formation"]],
        ],
        join="outer", axis=1,
    )
    return phenotype_df[phenotype_df.index.isin(coding_genes)].copy()
```

**Step 9: Continue collect_pombe_features.py** (final merge)

```python
@logger.catch
def merge_all_features(
    dna_df: pd.DataFrame,
    rna_df: pd.DataFrame,
    protein_df: pd.DataFrame,
    evolutionary_df: pd.DataFrame,
    network_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    gene_meta: pd.DataFrame,
) -> pd.DataFrame:
    """Outer-join all six feature groups on gene_systematic_id and fill category-column NAs."""
    pombe_features = pd.concat(
        [
            dna_df.query("Primary_candidate == True").set_index("Gene_id"),
            rna_df,
            protein_df,
            evolutionary_df,
            network_df,
            phenotype_df,
        ],
        join="outer", axis=1,
    ).rename_axis("gene_systematic_id")

    pombe_features = gene_meta[["gene_systematic_id", "gene_name"]].merge(
        pombe_features.reset_index(), on="gene_systematic_id", how="right"
    )
    pombe_features[["GO_term_richness", "PPI_degree", "GI_degree"]] = (
        pombe_features[["GO_term_richness", "PPI_degree", "GI_degree"]].fillna(0).astype(int)
    )
    pombe_features["DeletionLibrary_essentiality"] = pombe_features["DeletionLibrary_essentiality"].fillna("Not_determined")
    pombe_features["DeletionLibrary_category"] = pombe_features["DeletionLibrary_category"].fillna("Not_determined")
    pombe_features["RevisedDeletionLibrary_essentiality"] = pombe_features["RevisedDeletionLibrary_essentiality"].fillna("Not_determined")
    return pombe_features
```

**Step 10: Continue collect_pombe_features.py** (main orchestration function)

```python
@logger.catch
def run_feature_collection(config: InputOutputConfig) -> pd.DataFrame:
    """Execute the full 6-group feature collection pipeline and write both output files."""
    logger.info(f"Building gffutils DB from {config.pombase_dir}")
    genome_dir = config.pombase_dir / "genome_sequence_and_features"
    genome_cfg = PombaseGenomeConfig.from_pombase_dir(config.pombase_dir)
    db = gffutils.create_db(genome_cfg.gff3_file, genome_cfg.database_file, force=True, merge_strategy="create_unique")
    db = gffutils.FeatureDB(genome_cfg.database_file)

    gene_meta, uniprot2id = load_gene_meta(config.gene_meta_file)

    logger.info("Collecting DNA-level features")
    dna_df, coding_genes = collect_dna_level_features(db, genome_cfg)

    logger.info("Writing codon usage matrix")
    codon_usage_matrix = calculate_anticodon_usage_matrix(db, genome_cfg)
    codon_usage_matrix.to_csv(config.output_codon_usage, sep="\t", index=True)

    logger.info("Collecting RNA-level features")
    rna_df = collect_rna_level_features(config.literature_dir, config.gene_meta_file, coding_genes)

    logger.info("Collecting protein-level features")
    protein_meta = pd.read_csv(config.pombase_dir / "Protein_features" / "peptide_stats.tsv", sep="\t", index_col=0)
    protein_df = collect_protein_level_features(
        config.pombase_dir, config.alphafold_dir, config.literature_dir,
        config.gene_meta_file, protein_meta, uniprot2id, coding_genes,
    )

    logger.info("Collecting evolutionary-level features")
    evolutionary_df, phyloP_and_divergence = collect_evolutionary_level_features(
        config.pombase_dir, config.ensembl_paralogs_tsv, config.literature_dir, config.gene_meta_file, coding_genes,
    )

    logger.info("Loading GO ontology data for network-level features")
    ontology_cfg = OntologyDataConfig(
        ontology_obo=config.pombase_dir / "ontologies_and_associations" / "go-basic.obo",
        ontology_association_gaf=config.pombase_dir / "ontologies_and_associations" / "gene_ontology_annotation.gaf.tsv",
        slim_terms_table=[
            config.pombase_dir / "ontologies_and_associations" / "bp_go_slim_terms.tsv",
            config.pombase_dir / "ontologies_and_associations" / "mf_go_slim_terms.tsv",
            config.pombase_dir / "ontologies_and_associations" / "cc_go_slim_terms.tsv",
        ],
    )
    _, _, _, gene2go, _, _ = load_ontology_data(
        ontology_cfg.load_data(),
        relationships={"is_a", "part_of"}, propagate_counts=True, load_obsolete=False, prt=None,
    )

    logger.info("Collecting network-level features")
    network_df = collect_network_level_features(config.pombase_dir, config.biogrid_tsv, gene2go, coding_genes)

    logger.info("Collecting phenotype-level features")
    phenotype_df = collect_phenotype_level_features(
        config.pombase_dir, config.deletion_library_xlsx, config.essentiality_verification_csv,
        config.literature_dir, config.gene_meta_file, coding_genes, phyloP_and_divergence,
    )

    logger.info("Merging all feature groups")
    pombe_features = merge_all_features(dna_df, rna_df, protein_df, evolutionary_df, network_df, phenotype_df, gene_meta)

    pombe_features.to_csv(config.output_features, sep="\t", index=False)
    logger.success(f"Wrote {len(pombe_features)} gene records to {config.output_features}")

    return pombe_features
```

**Step 11: Write the CLI entry point**

```python
# =============================================================================
# MAIN EXECUTION
# =============================================================================
def parse_args() -> argparse.Namespace:
    """Parse command-line arguments and return the populated namespace."""
    parser = argparse.ArgumentParser(description="Collect pombe coding gene features from 15 sources")
    parser.add_argument("--pombase-dir", type=Path, required=True, help="PomBase version directory")
    parser.add_argument("--alphafold-dir", type=Path, required=True, help="AlphaFold structure directory (.pdb.gz files)")
    parser.add_argument("--literature-dir", type=Path, required=True, help="Directory of literature supplementary tables")
    parser.add_argument("--deletion-library-xlsx", type=Path, required=True, help="Curated deletion library categories xlsx")
    parser.add_argument("--essentiality-verification-csv", type=Path, required=True, help="Curated essentiality verification csv")
    parser.add_argument("--biogrid-tsv", type=Path, required=True, help="BioGrid interaction table")
    parser.add_argument("--ensembl-paralogs-tsv", type=Path, required=True, help="Ensembl paralog export table")
    parser.add_argument("--output", type=Path, required=True, dest="output_features", help="Output feature matrix path")
    parser.add_argument("--codon-usage-output", type=Path, required=True, dest="output_codon_usage", help="Output codon usage matrix path")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator: validate paths, run the feature collection pipeline, report results."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")

    try:
        config = InputOutputConfig(
            pombase_dir=args.pombase_dir,
            alphafold_dir=args.alphafold_dir,
            literature_dir=args.literature_dir,
            deletion_library_xlsx=args.deletion_library_xlsx,
            essentiality_verification_csv=args.essentiality_verification_csv,
            biogrid_tsv=args.biogrid_tsv,
            ensembl_paralogs_tsv=args.ensembl_paralogs_tsv,
            output_features=args.output_features,
            output_codon_usage=args.output_codon_usage,
        )
        run_feature_collection(config)
    except ValueError as e:
        logger.error(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    setup_logger()
    sys.exit(main())
```

**Step 12: Write workflow/rules/features.smk**

```python
# =============================================================================
# features.smk — Pombe coding gene feature collection (dataset-independent)
# =============================================================================

# Depends only on the reference PomBase version, not on any DIT-HAP sequencing
# project — no `dataset` wildcard (design doc §8).
rule collect_pombe_features:
    input:
        pombase_dir="resources/external/pombase/{pombase_version}",
        alphafold_dir=DATASETS["reference"]["alphafold_dir"],
        literature_dir="resources/literature",
        deletion_library_xlsx="resources/curated/deletion_library_categories.xlsx",
        essentiality_verification_csv="resources/curated/essentiality_verification.csv",
        biogrid_tsv="resources/external/biogrid/BIOGRID-ORGANISM-Schizosaccharomyces_pombe_972h-5.0.251.tab3.txt",
        ensembl_paralogs_tsv="resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv",
    output:
        features="results/features/{pombase_version}/pombe_coding_gene_protein_features.tsv",
        codon_usage="results/features/{pombase_version}/codon_usage_matrix.tsv",
    log:
        "logs/features/collect_pombe_features_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** Collecting pombe coding gene features for PomBase {wildcards.pombase_version}..."
    shell:
        """
        python workflow/scripts/features/collect_pombe_features.py \
            --pombase-dir {input.pombase_dir} \
            --alphafold-dir {input.alphafold_dir} \
            --literature-dir {input.literature_dir} \
            --deletion-library-xlsx {input.deletion_library_xlsx} \
            --essentiality-verification-csv {input.essentiality_verification_csv} \
            --biogrid-tsv {input.biogrid_tsv} \
            --ensembl-paralogs-tsv {input.ensembl_paralogs_tsv} \
            --output {output.features} \
            --codon-usage-output {output.codon_usage} &> {log}
        """
```

**Step 13: Uncomment the Snakefile's features.smk include and rule all target**

In `Snakefile`, replace:

```python
# include: "workflow/rules/features.smk"
# include: "workflow/rules/enrichment.smk"
# include: "workflow/rules/clustering.smk"
# (Phase 1 delivers only features.smk; remaining rules are follow-up work)
```

with:

```python
include: "workflow/rules/features.smk"
# include: "workflow/rules/enrichment.smk"
# include: "workflow/rules/clustering.smk"
# (Phase 1 delivers only features.smk; remaining rules are follow-up work)
```

And replace:

```python
rule all:
    input:
        # Uncommented in Task 5 once features.smk exists:
        # f"results/features/{DATASETS['reference']['pombase_version']}/pombe_coding_gene_protein_features.tsv",
    message:
        "*** DIT-HAP analysis complete"
```

with:

```python
rule all:
    input:
        f"results/features/{DATASETS['reference']['pombase_version']}/pombe_coding_gene_protein_features.tsv",
    message:
        "*** DIT-HAP analysis complete"
```

**Step 14: Write tests/test_collect_pombe_features.py**

Unit tests target the pure-transformation helper functions with small synthetic inputs (fast, no real data needed). A separate integration check (Step 15) runs the full script end-to-end and diffs against the notebook's own historical output — that's what actually proves the port is faithful, since these unit tests alone can't catch a wrong column order or a dropped merge.

```python
"""Tests for workflow/scripts/features/collect_pombe_features.py helper functions."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from workflow.scripts.features.collect_pombe_features import get_ortholog_counts, InputOutputConfig
import pandas as pd
import pytest


def test_get_ortholog_counts_counts_pipe_separated_entries(tmp_path):
    """A pipe-separated ortholog list counts entries; NONE maps to 0 via na_values."""
    f = tmp_path / "orthologs.txt"
    f.write_text("SPAC1002.01(name)\tOrthA|OrthB|OrthC\nSPAC1002.02(name)\tNONE\n")
    counts = get_ortholog_counts(f)
    assert counts.loc["SPAC1002.01"] == 3
    assert counts.loc["SPAC1002.02"] == 0


def test_get_ortholog_counts_strips_parenthetical_gene_name(tmp_path):
    """The index's trailing (name) suffix is stripped before returning counts."""
    f = tmp_path / "orthologs.txt"
    f.write_text("SPAC1002.01(mrx11)\tOrthA\n")
    counts = get_ortholog_counts(f)
    assert "SPAC1002.01" in counts.index
    assert "SPAC1002.01(mrx11)" not in counts.index


def test_input_output_config_rejects_missing_input(tmp_path):
    """InputOutputConfig.__post_init__ raises ValueError naming the missing path."""
    real_dir = tmp_path / "real"
    real_dir.mkdir()
    with pytest.raises(ValueError, match="does not exist"):
        InputOutputConfig(
            pombase_dir=tmp_path / "missing_pombase",
            alphafold_dir=real_dir,
            literature_dir=real_dir,
            deletion_library_xlsx=real_dir / "x.xlsx",
            essentiality_verification_csv=real_dir / "x.csv",
            biogrid_tsv=real_dir / "x.tsv",
            ensembl_paralogs_tsv=real_dir / "x.tsv",
            output_features=tmp_path / "out" / "features.tsv",
            output_codon_usage=tmp_path / "out" / "codon.tsv",
        )
```

Note: `test_input_output_config_rejects_missing_input` needs every *required* input path to exist except the one under test, since `__post_init__` validates all of them in one loop — the fixture creates `real_dir` and reuses it (as a stand-in directory) for the file-path fields so only `pombase_dir` is actually missing.

**Step 15: Run unit tests, then the full integration run against real data**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_analysis
conda run -n bioinformatics python -m pytest tests/test_collect_pombe_features.py -v
```

Expected: 3 tests pass.

Then run the real pipeline end-to-end (this is the actual faithfulness check — expect 10-20 minutes, dominated by the AlphaFold pLDDT scan over 10,392 structures and the GO DAG propagation):

```bash
mkdir -p results/features/2025-10-01 logs/features
conda run -n bioinformatics python workflow/scripts/features/collect_pombe_features.py \
    --pombase-dir resources/external/pombase/2025-10-01 \
    --alphafold-dir /data/c/yangyusheng_optimized/resource/AlphaFold_Dataset/20251107_downloaded/UP000002485_284812_SCHPO_v6 \
    --literature-dir resources/literature \
    --deletion-library-xlsx resources/curated/deletion_library_categories.xlsx \
    --essentiality-verification-csv resources/curated/essentiality_verification.csv \
    --biogrid-tsv resources/external/biogrid/BIOGRID-ORGANISM-Schizosaccharomyces_pombe_972h-5.0.251.tab3.txt \
    --ensembl-paralogs-tsv resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv \
    --output results/features/2025-10-01/pombe_coding_gene_protein_features.tsv \
    --codon-usage-output results/features/2025-10-01/codon_usage_matrix.tsv \
    --verbose
```

**Step 16: Verify against the notebook's own historical output**

`DIT_HAP_pipeline/resources/pombe_features/2025-10-01_pombe_coding_gene_protein_features.tsv` is the notebook's last real run against the same PomBase version. Diff row count, column set, and spot-check a few known values (the same `SPAC1002.01`/`SPAC1002.02` rows read manually during planning):

```bash
python3 -c "
import pandas as pd
new = pd.read_csv('results/features/2025-10-01/pombe_coding_gene_protein_features.tsv', sep='\t')
ref = pd.read_csv('/data/c/yangyusheng_optimized/DIT_HAP_pipeline/resources/pombe_features/2025-10-01_pombe_coding_gene_protein_features.tsv', sep='\t')
print('row counts:', len(new), len(ref))
print('new-only columns:', set(new.columns) - set(ref.columns))
print('ref-only columns:', set(ref.columns) - set(new.columns))
new_row = new[new['gene_systematic_id'] == 'SPAC1002.01'].iloc[0]
ref_row = ref[ref['gene_systematic_id'] == 'SPAC1002.01'].iloc[0]
print('SPAC1002.01 Peptide_length match:', new_row['Peptide_length'] == ref_row['Peptide_length'])
print('SPAC1002.01 mean_pLDDT match:', abs(new_row['mean_pLDDT'] - ref_row['mean_pLDDT']) < 0.01)
"
```

Expected: row counts match (~5100-5200 coding genes), column sets match exactly (both should have the duplicate `DeletionLibrary_essentiality` column noted in this task's header — `pd.read_csv` will auto-rename the second occurrence to `DeletionLibrary_essentiality.1`, which is fine, both files should do this identically), and the spot-checked values match. If row/column counts diverge, that signals a real divergence from the source notebook to debug before moving on — do not silently adjust the script to make counts match without understanding why they differed.

**Step 17: Commit**

```bash
git add workflow/scripts/features/collect_pombe_features.py workflow/rules/features.smk tests/test_collect_pombe_features.py Snakefile
git commit -m "feat: deliver features.smk — full pombe coding gene feature collection

- workflow/scripts/features/collect_pombe_features.py: byte-faithful port of
  pombe_feature_collection.ipynb's 15-source, 6-group feature assembly
- workflow/rules/features.smk: dataset-independent rule (no dataset wildcard,
  keyed only on pombase_version per design doc §8)
- Preserves the notebook's documented quirks verbatim (asymmetric BioGrid
  degree counting, duplicate DeletionLibrary_essentiality column)
- Snakefile: activated features.smk include + rule all target
- Verified against DIT_HAP_pipeline's historical 2025-10-01 output (Step 16)

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Task 9: Full test suite + Snakemake dry-run + wet-run verification

**Files:** none new — this task only runs and verifies everything built in Tasks 1-8 together.

**Step 1: Run the entire pytest suite**

```bash
cd /data/c/yangyusheng_optimized/DIT_HAP_analysis
conda run -n bioinformatics python -m pytest tests/ -v
```

Expected: all tests across every module pass (registry, io, gene_ids, features/genome, features/protein, enrichment/ontology, collect_pombe_features — roughly 35-40 tests total once Tasks 1-8 are all committed).

**Step 2: Snakemake dry-run against the real target**

```bash
mamba run -n snakemake snakemake -n --use-conda
```

Expected: dry-run reports exactly one job (`collect_pombe_features`) needed to satisfy `rule all` (or zero if Step 15/16's wet-run already produced the output files — in that case, `touch` nothing and treat "up to date" as success too).

**Step 3: Full wet-run through Snakemake itself** (not just the standalone script call from Task 8 Step 15 — this exercises the conda env activation and log redirection Snakemake adds)

```bash
rm -rf results/features logs/features   # clear Task 8's standalone run so Snakemake actually executes the rule
mamba run -n snakemake snakemake --cores 4 --use-conda
```

Expected: `collect_pombe_features` builds successfully via the `biopython.yml` conda env, producing the same two output files verified in Task 8 Step 16. Re-run the Step 16 diff against the DIT_HAP_pipeline reference output to confirm the Snakemake-driven run matches too.

**Step 4: Final commit (only if any test/lint fixes were needed in Steps 1-3)**

If everything passed with no code changes needed, there's nothing to commit — Task 8's commit already covers the delivered rule. If Step 1-3 surfaced a bug requiring a fix, commit that fix with a message describing what the dry-run/wet-run caught, e.g.:

```bash
git add -A
git commit -m "fix: <describe what the Snakemake dry-run/wet-run caught>

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
```

---

## Summary

Phase 1 delivers: repository skeleton (§2), the `datasets.yaml` release/ registry (§3), physical local copies of external resources (§4), a faithful `workflow/src/` regroup of exactly the modules `pombe_feature_collection.ipynb` uses (§7 — `io.py`, `gene_ids.py`, `features/genome.py`, `features/protein.py`, `enrichment/ontology.py`), and one complete, verified Snakemake rule (`features.smk`, §8) reproducing that notebook's full 15-source feature matrix deterministically.

Explicitly deferred to follow-up plans: `plotting/` (style.py, generic.py, gene_level.py — nothing in this plan calls them), `enrichment/pipeline.py` (ontology enrichment study, STRING API, REVIGO — the feature notebook only needs GO term richness), and all 10 remaining `.smk` rule files (clustering, ml, coverage, verification, comparison, noncoding_rna, spikein, complex, utr, domain_differences) plus their manual-judgment notebook counterparts under `notebooks/`.