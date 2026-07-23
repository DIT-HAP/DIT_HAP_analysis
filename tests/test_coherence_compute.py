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
    assert "GO:1" in groups
    assert "GO:2" not in groups


def test_build_groups_drops_broad_term_by_max_term_genes():
    from compute_coherence import build_groups
    bg = _background([f"g{i}" for i in range(10)])
    long = _long([("go_bp", "GO:9", "broad", f"g{i}", f"g{i}", 999) for i in range(3)])
    groups = build_groups(bg, long, min_group_size=3, max_group_size=300, max_term_genes=500)
    assert "GO:9" not in groups


def test_build_groups_drops_group_over_max_size():
    from compute_coherence import build_groups
    bg = _background([f"g{i}" for i in range(10)])
    long = _long([("go_cc", "GO:5", "toobig", f"g{i}", f"g{i}", 5) for i in range(5)])
    groups = build_groups(bg, long, min_group_size=3, max_group_size=4, max_term_genes=500)
    assert "GO:5" not in groups  # 5 DR-members > max_group_size 4


def test_build_groups_keeps_group_at_exact_boundaries():
    from compute_coherence import build_groups
    bg = _background([f"g{i}" for i in range(10)])
    # 3 DR-members == min_group_size (kept); n_group_genes == max_term_genes (kept, filter uses >)
    long = _long([("go_cc", "GO:6", "edge", f"g{i}", f"g{i}", 500) for i in range(3)])
    groups = build_groups(bg, long, min_group_size=3, max_group_size=300, max_term_genes=500)
    assert "GO:6" in groups


def test_load_fitting_results_drops_inf_rows(tmp_path):
    from compute_coherence import load_fitting_results
    df = pd.DataFrame(
        {"DR": [0.5, 0.8, np.inf], "DL": [1.0, 2.0, 3.0]},
        index=["SPAC1", "SPAC2", "SPINF"],
    )
    df.index.name = "Systematic ID"
    p = tmp_path / "fitting_results.tsv"
    df.to_csv(p, sep="\t")
    bg = load_fitting_results(p, dr_threshold=0.3)
    assert "SPINF" not in set(bg["Systematic ID"])
    assert np.isfinite(bg[["norm_DR", "norm_DL"]].to_numpy()).all()
