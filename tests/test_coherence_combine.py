import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import pandas as pd


def _metrics(source, rows):
    """A minimal per-source metrics frame: (group_id, z_score) rows for `source`."""
    return pd.DataFrame(
        [{"source": source, "group_id": gid, "z_score": z, "p_value": 0.5, "p_fdr": 0.5}
         for gid, z in rows]
    )


def test_combine_concatenates_and_sorts_by_zscore(tmp_path):
    from combine_metrics import combine

    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    _metrics("go_cc", [("GO:1", -2.0), ("GO:2", 0.5)]).to_csv(a, sep="\t", index=False)
    _metrics("go_bp", [("GO:9", -3.0), ("GO:8", 1.0)]).to_csv(b, sep="\t", index=False)

    combined = combine([a, b])
    # All four rows present, sources preserved and distinguishable.
    assert len(combined) == 4
    assert set(combined["source"]) == {"go_cc", "go_bp"}
    # Sorted by z_score ascending (most coherent first).
    assert list(combined["z_score"]) == sorted(combined["z_score"])
    assert combined.iloc[0]["group_id"] == "GO:9"


def test_combine_tolerates_empty_source_tables(tmp_path):
    from combine_metrics import combine

    full = tmp_path / "full.tsv"
    empty = tmp_path / "empty.tsv"
    _metrics("go_cc", [("GO:1", -1.0)]).to_csv(full, sep="\t", index=False)
    # An empty source table (header only) — no group passed the size filter.
    _metrics("go_bp", []).to_csv(empty, sep="\t", index=False)

    combined = combine([empty, full])
    assert len(combined) == 1
    assert combined.iloc[0]["source"] == "go_cc"


def test_combine_all_empty_returns_empty_frame(tmp_path):
    from combine_metrics import combine

    e1 = tmp_path / "e1.tsv"
    e2 = tmp_path / "e2.tsv"
    _metrics("go_cc", []).to_csv(e1, sep="\t", index=False)
    _metrics("go_bp", []).to_csv(e2, sep="\t", index=False)

    combined = combine([e1, e2])
    assert combined.empty
