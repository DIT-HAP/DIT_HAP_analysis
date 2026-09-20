import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import pandas as pd

from io_table import read_parquet, write_parquet


def _metrics(source, rows):
    """A minimal per-source metrics frame: (group_id, z-score) rows for `source`."""
    return pd.DataFrame(
        [{"source": source, "group_id": gid,
          "median_pairwise_distance_z": z, "median_pairwise_distance_p": 0.5, "q_value": 0.5}
         for gid, z in rows]
    )


def test_combine_concatenates_and_sorts_by_zscore(tmp_path):
    from combine_metrics import combine

    a = tmp_path / "a.parquet"
    b = tmp_path / "b.parquet"
    write_parquet(_metrics("go_cc", [("GO:1", -2.0), ("GO:2", 0.5)]), a)
    write_parquet(_metrics("go_bp", [("GO:9", -3.0), ("GO:8", 1.0)]), b)

    combined = combine([a, b])
    # All four rows present, sources preserved and distinguishable.
    assert len(combined) == 4
    assert set(combined["source"]) == {"go_cc", "go_bp"}
    # Sorted by median_pairwise_distance_z ascending (most coherent first).
    assert list(combined["median_pairwise_distance_z"]) == sorted(combined["median_pairwise_distance_z"])
    assert combined.iloc[0]["group_id"] == "GO:9"


def test_combine_tolerates_empty_source_tables(tmp_path):
    from combine_metrics import combine

    full = tmp_path / "full.parquet"
    empty = tmp_path / "empty.parquet"
    write_parquet(_metrics("go_cc", [("GO:1", -1.0)]), full)
    # An empty source table (no rows) — no group passed the size filter.
    write_parquet(_metrics("go_bp", []), empty)

    combined = combine([empty, full])
    assert len(combined) == 1
    assert combined.iloc[0]["source"] == "go_cc"


def test_combine_all_empty_returns_empty_frame(tmp_path):
    from combine_metrics import combine

    e1 = tmp_path / "e1.parquet"
    e2 = tmp_path / "e2.parquet"
    write_parquet(_metrics("go_cc", []), e1)
    write_parquet(_metrics("go_bp", []), e2)

    combined = combine([e1, e2])
    assert combined.empty


def test_zero_column_frame_round_trips_through_parquet(tmp_path):
    """The all-empty case writes a 0x0 frame, which must survive the round-trip.

    This is the contract that replaced the TSV path: on a fully empty source set
    the old combine wrote a header-less TSV that deduplicate_terms then failed to
    read at all (EmptyDataError -> exit 1), so an all-empty dataset stopped the
    pipeline instead of producing empty tables.
    """
    path = tmp_path / "empty.parquet"
    write_parquet(pd.DataFrame(), path)
    assert read_parquet(path).empty
