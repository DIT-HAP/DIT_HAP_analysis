"""Tests for the split gene-level clustering pipeline (shared module + driver configs)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.src.io import write_parquet, read_parquet

from workflow.src.clustering.candidates import (
    BEST_METHOD,
    DR_CAP,
    DL_DIVISOR,
    METHODS,
    cluster_one_method,
    finalize_auto_merge,
    finalize_direct,
    finalize_grid,
    renumber_by_dr,
    scale_features,
    score_labels,
)
from workflow.scripts.clustering.prepare_clustering_data import PrepareConfig
from workflow.scripts.clustering.cluster_one_method import MethodConfig
from workflow.scripts.clustering.select_candidate_clusters import SelectConfig, combine_metrics
from workflow.scripts.clustering.finalize_direct_clusters import FinalizeDirectConfig, run as run_finalize_direct
from workflow.scripts.clustering.finalize_auto_merge_clusters import FinalizeAutoMergeConfig, run as run_finalize_auto_merge
from workflow.scripts.clustering.finalize_grid_clusters import FinalizeGridConfig, run as run_finalize_grid


def test_best_method_is_pinned_to_kmeans():
    """The fragile set()[0] selection is replaced by an explicit pin (quirk #3)."""
    assert BEST_METHOD == "kmeans"
    assert METHODS[0] == "kmeans"


def test_scale_features_caps_DR_and_divides_DL():
    """DR above 1.3 is clamped to 1.3; DL is divided by 10 (byte-faithful quirk #1)."""
    df = pd.DataFrame(
        {"DR": [0.5, 1.3, 2.0, 1.29], "DL": [10.0, 20.0, 5.0, 100.0]},
        index=["g1", "g2", "g3", "g4"],
    )
    scaled = scale_features(df, ["DR", "DL"])
    # DR: 0.5 stays, 1.3 stays (strict < so equal is unchanged), 2.0 -> 1.3, 1.29 stays.
    assert list(scaled["DR"]) == [0.5, 1.3, DR_CAP, 1.29]
    # DL divided by 10.
    assert list(scaled["DL"]) == [1.0, 2.0, 0.5, 10.0]


def test_scale_features_respects_custom_dr_cap_and_dl_divisor():
    """Non-default dr_cap/dl_divisor (as would come from config/analysis.yaml) are honored."""
    df = pd.DataFrame(
        {"DR": [0.5, 1.0, 2.0], "DL": [10.0, 20.0, 5.0]},
        index=["g1", "g2", "g3"],
    )
    scaled = scale_features(df, ["DR", "DL"], dr_cap=1.0, dl_divisor=5)
    assert list(scaled["DR"]) == [0.5, 1.0, 1.0]
    assert list(scaled["DL"]) == [2.0, 4.0, 1.0]


def test_scale_features_drops_nan_rows():
    """dropna defines the clustered set — genes with NaN DR/DL are excluded."""
    df = pd.DataFrame(
        {"DR": [0.5, np.nan, 0.8], "DL": [10.0, 20.0, np.nan]},
        index=["g1", "g2", "g3"],
    )
    scaled = scale_features(df, ["DR", "DL"])
    assert list(scaled.index) == ["g1"]


def test_cluster_one_method_all_four_return_0based_labels():
    """Each of the four methods returns 0-based labels for the requested cluster count."""
    rng = np.random.default_rng(0)
    data = pd.DataFrame(rng.random((60, 2)), columns=["DR", "DL"])
    for method in METHODS:
        labels = cluster_one_method(method, data, n_clusters=4, random_state=42)
        assert labels.min() == 0, method
        assert len(np.unique(labels)) == 4, method


def test_cluster_one_method_rejects_unknown_method():
    """cluster_one_method raises ValueError on an unrecognized method name."""
    data = pd.DataFrame(np.random.default_rng(0).random((20, 2)), columns=["DR", "DL"])
    with pytest.raises(ValueError, match="Unknown clustering method"):
        cluster_one_method("dbscan", data, n_clusters=3, random_state=42)


def test_score_labels_reports_three_scores_and_count():
    """score_labels returns the three sklearn scores plus the cluster count."""
    rng = np.random.default_rng(1)
    data = pd.DataFrame(rng.random((60, 2)), columns=["DR", "DL"])
    labels = cluster_one_method("kmeans", data, n_clusters=3, random_state=42)
    scores = score_labels(data, labels)
    assert set(scores) == {"silhouette_score", "calinski_harabasz_score", "davies_bouldin_score", "n_clusters"}
    assert scores["n_clusters"] == 3


def test_kmeans_is_deterministic_given_seed():
    """Same seed + same data -> identical kmeans labels (reproducibility guarantee)."""
    rng = np.random.default_rng(2)
    data = pd.DataFrame(rng.random((80, 2)), columns=["DR", "DL"])
    r1 = cluster_one_method("kmeans", data, n_clusters=5, random_state=42)
    r2 = cluster_one_method("kmeans", data, n_clusters=5, random_state=42)
    assert np.array_equal(r1, r2)


def test_combine_metrics_stacks_ksweep_then_rounded_methods():
    """combine_metrics concatenates the unrounded k-sweep and 3dp-rounded method rows."""
    ksweep = pd.DataFrame({"k": [2], "inertia": [1.23456], "silhouette": [0.5]})
    methods = pd.DataFrame(
        {"method": ["kmeans"], "silhouette_score": [0.4531111], "n_clusters": [64]}
    )
    combined = combine_metrics(ksweep, methods)
    assert list(combined["table"]) == ["k_sweep", "method_comparison"]
    # k-sweep row keeps full precision; method row is rounded to 3 dp.
    assert combined.loc[0, "inertia"] == 1.23456
    assert combined.loc[1, "silhouette_score"] == 0.453


def test_prepare_config_rejects_missing_input(tmp_path):
    """PrepareConfig.validate raises ValueError when a required input file is absent."""
    cfg = PrepareConfig(
        fitting_results=tmp_path / "nope.tsv",
        essentiality_verification_csv=tmp_path / "also_nope.csv",
        output_annotated=tmp_path / "w" / "ann.parquet",
        output_scaled=tmp_path / "w" / "scaled.parquet",
        output_ksweep=tmp_path / "w" / "ksweep.parquet",
    )
    with pytest.raises(ValueError, match="Required input not found"):
        cfg.validate()


def test_finalize_direct_config_rejects_missing_input(tmp_path):
    """FinalizeDirectConfig.validate raises ValueError when a required input file is absent."""
    cfg = FinalizeDirectConfig(
        annotated_data=tmp_path / "nope.parquet",
        scaled_data=tmp_path / "also_nope.parquet",
        output=tmp_path / "out" / "final_clusters.tsv",
    )
    with pytest.raises(ValueError, match="Required input not found"):
        cfg.validate()


def test_method_config_rejects_unknown_method(tmp_path):
    """MethodConfig.validate raises ValueError on an unrecognized method name."""
    scaled = tmp_path / "scaled.parquet"
    scaled.write_bytes(b"")
    cfg = MethodConfig(
        method="spectral",
        scaled_data=scaled,
        output_labels=tmp_path / "w" / "l.parquet",
        output_metrics=tmp_path / "w" / "m.parquet",
    )
    with pytest.raises(ValueError, match="Unknown clustering method"):
        cfg.validate()


def test_select_config_rejects_missing_input(tmp_path):
    """SelectConfig.validate raises ValueError when a required input is absent."""
    cfg = SelectConfig(
        annotated_data=tmp_path / "ann.parquet",
        ksweep=tmp_path / "ksweep.parquet",
        best_labels=tmp_path / "labels.parquet",
        method_metrics=[tmp_path / "m.parquet"],
        output_clusters=tmp_path / "out" / "candidate_clusters.tsv",
        output_metrics=tmp_path / "out" / "clustering_metrics.tsv",
    )
    with pytest.raises(ValueError, match="Required input not found"):
        cfg.validate()


def _toy_annotated_scaled(n_per=20, seed=0):
    """Build a small annotated table + matching scaled matrix with 9 well-separated blobs in (DR, DL)."""
    rng = np.random.default_rng(seed)
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


def test_finalize_direct_produces_k_clusters_labelled_1_to_9():
    annotated, scaled = _toy_annotated_scaled()
    out = finalize_direct(annotated, scaled, method="kmeans", n_clusters=9, random_state=42, wt_cluster=9)
    assert "cluster" in out.columns
    assert "raw_cluster" not in out.columns          # direct has no pre-merge labels
    assert sorted(out["cluster"].unique()) == list(range(1, 10))
    assert out.index.name == "Systematic ID"
    assert {"DR", "DL", "A"}.issubset(out.columns)


def test_finalize_direct_assigns_lowest_DR_to_wt_cluster():
    annotated, scaled = _toy_annotated_scaled()
    out = finalize_direct(annotated, scaled, method="kmeans", n_clusters=9, random_state=42, wt_cluster=9)
    means = out.groupby("cluster")["DR"].mean()
    assert means.idxmin() == 9                       # WT = lowest DR
    non_wt = means.drop(index=9).sort_index()
    assert non_wt.is_monotonic_increasing            # ids 1..8 ascend in DR


def test_finalize_direct_is_deterministic():
    annotated, scaled = _toy_annotated_scaled()
    a = finalize_direct(annotated, scaled, method="kmeans", n_clusters=9, random_state=42, wt_cluster=9)
    b = finalize_direct(annotated, scaled, method="kmeans", n_clusters=9, random_state=42, wt_cluster=9)
    pd.testing.assert_series_equal(a["cluster"], b["cluster"])


def test_finalize_direct_only_labels_scaled_genes():
    """Genes dropped by scaling (NaN DR/DL) get no cluster (NaN)."""
    annotated, scaled = _toy_annotated_scaled(n_per=15)
    extra = pd.DataFrame({"DR": [np.nan], "DL": [np.nan], "A": [1.0]}, index=["ghost"])
    extra.index.name = "Systematic ID"
    annotated2 = pd.concat([annotated, extra])
    out = finalize_direct(annotated2, scaled, method="kmeans", n_clusters=9, random_state=42, wt_cluster=9)
    assert pd.isna(out.loc["ghost", "cluster"])
    assert out.loc["g0_0", "cluster"] in range(1, 10)


def test_finalize_direct_defaults_to_best_method():
    """Omitting `method` uses BEST_METHOD (kmeans) — byte-identical to the old auto path."""
    annotated, scaled = _toy_annotated_scaled()
    a = finalize_direct(annotated, scaled, n_clusters=9, random_state=42, wt_cluster=9)
    b = finalize_direct(annotated, scaled, method=BEST_METHOD, n_clusters=9, random_state=42, wt_cluster=9)
    pd.testing.assert_series_equal(a["cluster"], b["cluster"])


def test_renumber_by_dr_tiebreak_on_dl_then_raw_id():
    """Two groups with EQUAL mean DR are ordered by the mean_dl secondary key,
    deterministically and reproducibly (design doc §2 tie-break rule)."""
    # Geometry is hand-built so the tie is exact:
    #   A: DR=0.0 (clearly lowest -> WT)
    #   B: DR=1.0, DL=0.2   } identical mean DR, differing mean DL -> only the
    #   C: DR=1.0, DL=0.8   } mean_dl tiebreak can order these two.
    rows = []
    for i in range(10):
        rows.append((f"A{i}", 0.0, 0.5, "A"))
    for i in range(10):
        rows.append((f"B{i}", 1.0, 0.2, "B"))
    for i in range(10):
        rows.append((f"C{i}", 1.0, 0.8, "C"))
    idx = [r[0] for r in rows]
    annotated = pd.DataFrame(
        {"DR": [r[1] for r in rows], "DL": [r[2] for r in rows], "A": 1.0}, index=idx
    )
    annotated.index.name = "Systematic ID"
    raw = pd.Series([r[3] for r in rows], index=idx)

    final1 = renumber_by_dr(annotated, raw, n_clusters=3, wt_cluster=3)
    final2 = renumber_by_dr(annotated, raw, n_clusters=3, wt_cluster=3)
    pd.testing.assert_series_equal(final1, final2)          # reproducible

    assert final1.loc["A0"] == 3                             # lowest DR -> WT id 3
    b_id, c_id = final1.loc["B0"], final1.loc["C0"]
    assert b_id < c_id                                       # lower mean DL -> lower final id


def test_renumber_by_dr_raises_on_group_count_mismatch():
    """Fewer groups than n_clusters raises rather than silently truncating."""
    annotated = pd.DataFrame(
        {"DR": [0.0, 0.0, 1.0, 1.0], "DL": [0.1, 0.2, 0.3, 0.4], "A": 1.0},
        index=["a", "b", "c", "d"],
    )
    raw = pd.Series(["x", "x", "y", "y"], index=["a", "b", "c", "d"])
    with pytest.raises(ValueError, match="Expected 9 groups but got 2"):
        renumber_by_dr(annotated, raw, n_clusters=9, wt_cluster=9)


def test_finalize_auto_merge_recovers_known_groups():
    """Ward-merging 64 candidate centroids drawn from 9 well-separated blobs recovers
    the 9 groups, keeps raw_cluster, and DR-numbers with WT=9."""
    annotated, scaled = _toy_annotated_scaled(n_per=20)
    raw64 = pd.Series(
        cluster_one_method("kmeans", scaled, n_clusters=64, random_state=42),
        index=scaled.index,
    )
    out = finalize_auto_merge(annotated, scaled, raw64, n_clusters=9, wt_cluster=9)
    assert sorted(out["cluster"].dropna().unique()) == list(range(1, 10))
    assert "raw_cluster" in out.columns
    assert out["raw_cluster"].nunique() == 64
    means = out.groupby("cluster")["DR"].mean()
    assert means.idxmin() == 9
    assert means.drop(index=9).sort_index().is_monotonic_increasing


def test_finalize_auto_merge_raises_when_too_few_candidates():
    annotated, scaled = _toy_annotated_scaled(n_per=5)
    raw = pd.Series(
        cluster_one_method("kmeans", scaled, n_clusters=5, random_state=42),
        index=scaled.index,
    )
    with pytest.raises(ValueError, match="needs >= 9 candidate clusters"):
        finalize_auto_merge(annotated, scaled, raw, n_clusters=9, wt_cluster=9)


def _toy_grid_annotated_scaled(n_per=12, seed=0):
    """Toy data laid out ON a 3x3 (DR, DL) grid — one blob per cell — so axis cuts
    at dr=[0.4,0.8], dl=[0.4,0.8] fill all 9 cells (grid partitions rectangles, so
    the arbitrary 9-blob layout of _toy_annotated_scaled would leave cells empty)."""
    rng = np.random.default_rng(seed)
    drs, dls, idx = [], [], []
    for r, dr in enumerate((0.2, 0.6, 1.0)):
        for c, dl in enumerate((0.2, 0.6, 1.0)):
            for i in range(n_per):
                drs.append(dr + rng.normal(0, 0.01))
                dls.append(dl + rng.normal(0, 0.01))
                idx.append(f"r{r}c{c}_{i}")
    annotated = pd.DataFrame({"DR": drs, "DL": dls, "A": 1.0}, index=idx)
    annotated.index.name = "Systematic ID"
    scaled = annotated[["DR", "DL"]].copy()
    return annotated, scaled


def test_finalize_grid_assigns_cells_and_numbers_by_dr():
    """A 3x3 grid on grid-arranged data yields 9 clusters, DR-numbered, WT=9."""
    annotated, scaled = _toy_grid_annotated_scaled()
    out = finalize_grid(
        annotated, scaled, dr_cuts=[0.4, 0.8], dl_cuts=[0.4, 0.8],
        n_clusters=9, wt_cluster=9,
    )
    assert sorted(out["cluster"].dropna().unique()) == list(range(1, 10))
    assert out.groupby("cluster")["DR"].mean().idxmin() == 9


def test_finalize_grid_raises_on_cell_count_mismatch():
    annotated, scaled = _toy_annotated_scaled(n_per=10)
    with pytest.raises(ValueError, match="but final_n_clusters is 9"):
        finalize_grid(annotated, scaled, dr_cuts=[0.5], dl_cuts=[0.5], n_clusters=9, wt_cluster=9)


def test_finalize_direct_driver_writes_final_tsv(tmp_path):
    annotated, scaled = _toy_annotated_scaled()
    ap, sp = tmp_path / "annotated.parquet", tmp_path / "scaled.parquet"
    write_parquet(annotated, ap)
    write_parquet(scaled, sp)
    out = tmp_path / "final_clusters.tsv"
    cfg = FinalizeDirectConfig(
        annotated_data=ap, scaled_data=sp, output=out,
        method="kmeans", n_clusters=9, random_state=42, wt_cluster=9,
    )
    run_finalize_direct(cfg)
    df = pd.read_csv(out, sep="\t")
    assert "Systematic ID" in df.columns
    assert "cluster" in df.columns
    assert "raw_cluster" not in df.columns
    assert sorted(df["cluster"].dropna().unique()) == list(range(1, 10))


def test_finalize_auto_merge_driver_writes_final_tsv(tmp_path):
    annotated, scaled = _toy_annotated_scaled()
    raw64 = pd.Series(
        cluster_one_method("kmeans", scaled, n_clusters=64, random_state=42),
        index=scaled.index, name="cluster",
    )
    ap, sp, lp = tmp_path / "annotated.parquet", tmp_path / "scaled.parquet", tmp_path / "kmeans_labels.parquet"
    write_parquet(annotated, ap)
    write_parquet(scaled, sp)
    write_parquet(raw64, lp)
    out = tmp_path / "final_clusters.tsv"
    cfg = FinalizeAutoMergeConfig(
        annotated_data=ap, scaled_data=sp, candidate_labels=lp, output=out,
        n_clusters=9, wt_cluster=9,
    )
    run_finalize_auto_merge(cfg)
    df = pd.read_csv(out, sep="\t")
    assert "cluster" in df.columns
    assert "raw_cluster" in df.columns
    assert sorted(df["cluster"].dropna().unique()) == list(range(1, 10))


def test_finalize_grid_driver_writes_final_tsv(tmp_path):
    annotated, scaled = _toy_grid_annotated_scaled()
    ap, sp = tmp_path / "annotated.parquet", tmp_path / "scaled.parquet"
    write_parquet(annotated, ap)
    write_parquet(scaled, sp)
    out = tmp_path / "final_clusters.tsv"
    cfg = FinalizeGridConfig(
        annotated_data=ap, scaled_data=sp, output=out,
        dr_cuts=[0.4, 0.8], dl_cuts=[0.4, 0.8], n_clusters=9, wt_cluster=9,
    )
    run_finalize_grid(cfg)
    df = pd.read_csv(out, sep="\t")
    assert "cluster" in df.columns
    assert sorted(df["cluster"].dropna().unique()) == list(range(1, 10))
