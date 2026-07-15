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
