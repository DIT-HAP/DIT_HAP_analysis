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
