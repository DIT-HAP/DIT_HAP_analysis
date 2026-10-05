import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import pandas as pd
import pytest

from io_table import write_parquet


def _metrics(rows):
    """A minimal metrics frame: (source, group_id, z, q) rows."""
    return pd.DataFrame(
        [{"source": src, "group_id": gid, "group_name": f"name {gid}",
          "median_pairwise_distance_z": z, "q_value": q,
          "scored_member_names": ["ada1", "bms1"]}
         for src, gid, z, q in rows]
    )


def _config(tmp_path, **overrides):
    from export_cohorts import ExportConfig

    return ExportConfig(
        metrics=tmp_path / "metrics.parquet",
        output=tmp_path / "cohorts.xlsx",
        **overrides,
    )


def test_labels_apply_both_conditions_to_the_coherent_side(tmp_path):
    """Coherent needs q <= 0.05 AND z < -2; incoherent is z > 1 alone."""
    from export_cohorts import Cohort, label_cohorts

    table = _metrics([
        ("go_cc", "GO:1", -3.0, 0.01),   # coherent: both cuts pass
        ("go_cc", "GO:2", -3.0, 0.40),   # tight but not significant -> other
        ("go_cc", "GO:3", -1.0, 0.01),   # significant but not tight enough -> other
        ("go_cc", "GO:4", 2.0, 0.43),    # incoherent: z only, no q requirement
        ("go_cc", "GO:5", 0.5, 0.90),    # between the cuts -> other
    ])

    labelled = label_cohorts(table, _config(tmp_path))

    assert dict(zip(labelled["group_id"], labelled["cohort"])) == {
        "GO:1": Cohort.COHERENT.value,
        "GO:2": Cohort.OTHER.value,
        "GO:3": Cohort.OTHER.value,
        "GO:4": Cohort.INCOHERENT.value,
        "GO:5": Cohort.OTHER.value,
    }


def test_source_filter_keeps_one_source(tmp_path):
    """--source slices the table, so a per-source workbook is a cut of the combined one."""
    from export_cohorts import read_metrics

    metrics = tmp_path / "metrics.parquet"
    write_parquet(_metrics([("go_cc", "GO:1", -3.0, 0.01), ("go_bp", "GO:2", -3.0, 0.01)]), metrics)

    sliced = read_metrics(_config(tmp_path, source="go_cc"))

    assert list(sliced["group_id"]) == ["GO:1"]


def test_dedup_flag_joins_on_source_and_group_id(tmp_path):
    """group_id is not unique across sources, so the join must key on (source, group_id).

    A representatives-only table has no is_representative column: presence in it IS
    the flag, which is the original contract this keeps.
    """
    from export_cohorts import attach_dedup_columns

    metrics = _metrics([("go_cc", "GO:1", -3.0, 0.01), ("go_bp", "GO:1", -3.0, 0.01)])
    representatives = tmp_path / "reps.tsv"
    # Only go_bp's GO:1 survived de-duplication.
    pd.DataFrame({"source": ["go_bp"], "group_id": ["GO:1"]}).to_csv(representatives, sep="\t", index=False)

    flagged = attach_dedup_columns(metrics, representatives).set_index("source")

    assert flagged.loc["go_bp", "in_dedup_set"]
    assert not flagged.loc["go_cc", "in_dedup_set"]


def test_dedup_columns_carry_moonlighting_fraction_from_the_all_terms_table(tmp_path):
    """The all-terms table drives the flag from is_representative and adds the fraction."""
    from export_cohorts import attach_dedup_columns

    metrics = _metrics([("go_cc", "GO:1", -3.0, 0.01), ("go_cc", "GO:2", -3.0, 0.01)])
    dedup = tmp_path / "dedup.tsv"
    pd.DataFrame({
        "source": ["go_cc", "go_cc"],
        "group_id": ["GO:1", "GO:2"],
        "is_representative": [True, False],
        "moonlighting_fraction": [0.25, 1.0],
    }).to_csv(dedup, sep="\t", index=False)

    merged = attach_dedup_columns(metrics, dedup).set_index("group_id")

    assert merged.loc["GO:1", "in_dedup_set"]
    assert not merged.loc["GO:2", "in_dedup_set"]
    assert merged["moonlighting_fraction"].to_dict() == {"GO:1": 0.25, "GO:2": 1.0}


def test_workbook_holds_one_sheet_per_cohort(tmp_path):
    pytest.importorskip("openpyxl")
    from export_cohorts import build_sheets, label_cohorts, write_workbook

    config = _config(tmp_path)
    table = _metrics([
        ("go_cc", "GO:1", -3.0, 0.01),
        ("go_cc", "GO:2", 2.0, 0.43),
        ("go_cc", "GO:3", 0.5, 0.90),
    ])

    write_workbook(build_sheets(label_cohorts(table, config), config), config.output)

    sheets = pd.read_excel(config.output, sheet_name=None)
    assert list(sheets) == ["coherent", "incoherent", "all_terms", "thresholds"]
    assert list(sheets["coherent"]["group_id"]) == ["GO:1"]
    assert list(sheets["incoherent"]["group_id"]) == ["GO:2"]
    assert len(sheets["all_terms"]) == 3
    assert {"coherent_q_max", "coherent_z_threshold", "incoherent_z_threshold"} <= set(sheets["thresholds"]["cut"])
    # The list column is joined for display, not written as an Excel-hostile repr.
    assert sheets["coherent"]["scored_member_names"].iloc[0] == "ada1, bms1"
