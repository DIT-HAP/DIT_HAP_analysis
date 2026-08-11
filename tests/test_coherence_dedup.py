import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import itertools

import pandas as pd
import pytest

from deduplicate_terms import (
    UnionFind,
    overlap_coefficient,
    candidate_pairs,
    build_clusters,
    deduplicate,
    DedupConfig,
)


# --- primitives -------------------------------------------------------------
def test_overlap_coefficient_nested_is_one():
    """A subset nested in a superset has overlap coefficient 1.0 (the GO parent/child case)."""
    parent = {"a", "b", "c", "d"}
    child = {"a", "b"}
    assert overlap_coefficient(parent, child) == 1.0
    assert overlap_coefficient(set(), child) == 0.0


def test_overlap_coefficient_partial():
    assert overlap_coefficient({"a", "b", "c"}, {"b", "c", "d"}) == pytest.approx(2 / 3)


def test_candidate_pairs_matches_bruteforce():
    """The gene->terms inverted index yields exactly the member-sharing pairs."""
    sets = [{"a", "b"}, {"b", "c"}, {"x", "y"}, {"c"}]
    got = candidate_pairs(sets)
    brute = {
        (i, j) for i, j in itertools.combinations(range(len(sets)), 2)
        if sets[i] & sets[j]
    }
    assert got == brute
    assert (2, 3) not in got  # {x,y} shares nothing with {c}


def test_union_find_transitive():
    """A-B and B-C unions put A, B, C in one set even though A,C never met directly."""
    uf = UnionFind(4)
    uf.union(0, 1)
    uf.union(1, 2)
    groups = uf.groups()
    assert len(groups) == 2  # {0,1,2} and {3}
    roots = {uf.find(0), uf.find(1), uf.find(2)}
    assert len(roots) == 1
    assert uf.find(3) not in roots


# --- clustering -------------------------------------------------------------
def _sub(rows):
    """rows: list of (group_id, covered_genes-space-set). Minimal cols for build_clusters."""
    return pd.DataFrame(
        [{"group_id": gid, "covered_genes": ", ".join(sorted(genes))} for gid, genes in rows]
    )


def test_build_clusters_merges_high_overlap():
    """Two >=threshold-overlapping terms land in one cluster; a disjoint term stays alone."""
    sub = _sub([
        ("GO:1", {"a", "b", "c", "d"}),
        ("GO:2", {"a", "b", "c", "e"}),   # 3/4 overlap with GO:1 -> merge at 0.5
        ("GO:3", {"x", "y", "z"}),        # disjoint
    ])
    labels = build_clusters(sub, threshold=0.5, merge_dag_lineage=False, ancestors={})
    assert labels[0] == labels[1]
    assert labels[2] != labels[0]


def test_build_clusters_lineage_only_merges_member_sharing():
    """DAG lineage unites an ancestor/descendant pair ONLY when they share a member.

    GO:child shares one gene with GO:parent (overlap 1/3 < 0.5, so overlap alone
    would NOT merge), but the lineage rule merges them. A same-lineage but
    member-disjoint term must stay separate.
    """
    sub = _sub([
        ("GO:parent", {"a", "b", "c"}),
        ("GO:child", {"a", "m", "n"}),      # shares 'a' with parent; overlap 1/3
        ("GO:cousin", {"p", "q", "r"}),     # in lineage but no shared member
    ])
    ancestors = {"GO:child": {"GO:parent"}, "GO:cousin": {"GO:parent"}, "GO:parent": set()}
    labels = build_clusters(sub, threshold=0.5, merge_dag_lineage=True, ancestors=ancestors)
    assert labels[0] == labels[1]          # parent + child merged via lineage
    assert labels[2] != labels[0]          # cousin shares no member -> not merged


def test_build_clusters_lineage_off_keeps_low_overlap_separate():
    """With lineage off, the low-overlap parent/child pair stays in separate clusters."""
    sub = _sub([
        ("GO:parent", {"a", "b", "c"}),
        ("GO:child", {"a", "m", "n"}),
    ])
    ancestors = {"GO:child": {"GO:parent"}, "GO:parent": set()}
    labels = build_clusters(sub, threshold=0.5, merge_dag_lineage=False, ancestors=ancestors)
    assert labels[0] != labels[1]


# --- orchestration: representative selection --------------------------------
def _combined(rows):
    """rows: dicts with source, group_id, group_name, term_size, covered_genes, z_score, p_fdr."""
    return pd.DataFrame(rows)


def _cfg(tmp_path, **kw):
    return DedupConfig(
        combined=tmp_path / "c.tsv", obo=tmp_path / "o.obo",
        output_all=tmp_path / "all.tsv", output_representatives=tmp_path / "rep.tsv",
        **kw,
    )


def test_deduplicate_picks_best_qvalue_representative(tmp_path):
    """Within a redundant cluster the min-p_fdr term is the representative; all rows kept."""
    table = _combined([
        {"source": "go_bp", "group_id": "GO:1", "group_name": "big", "term_size": 200,
         "covered_genes": "a, b, c, d", "z_score": -5.0, "p_fdr": 0.05},
        {"source": "go_bp", "group_id": "GO:2", "group_name": "tight", "term_size": 20,
         "covered_genes": "a, b, c, e", "z_score": -7.0, "p_fdr": 0.01},  # best q
        {"source": "go_cc", "group_id": "GO:9", "group_name": "other", "term_size": 5,
         "covered_genes": "x, y, z", "z_score": -3.0, "p_fdr": 0.2},
    ])
    depth = {"GO:1": 3, "GO:2": 6, "GO:9": 4}
    ancestors = {"GO:1": set(), "GO:2": set(), "GO:9": set()}
    out = deduplicate(table, _cfg(tmp_path, scope="pooled", merge_dag_lineage=False), depth, ancestors)
    assert len(out) == 3  # nothing dropped
    reps = out[out["is_representative"]]
    # GO:1 & GO:2 are one cluster (overlap 3/4); GO:9 alone -> 2 clusters, 2 reps.
    assert set(reps["group_id"]) == {"GO:2", "GO:9"}
    cl = out[out["group_id"].isin(["GO:1", "GO:2"])]
    assert cl["representative_group_id"].nunique() == 1
    assert cl["representative_group_id"].iloc[0] == "GO:2"


def test_deduplicate_force_representative_overrides(tmp_path):
    """A forced group_id becomes its cluster's representative even with a worse q."""
    table = _combined([
        {"source": "go_bp", "group_id": "GO:1", "group_name": "big", "term_size": 200,
         "covered_genes": "a, b, c, d", "z_score": -5.0, "p_fdr": 0.05},
        {"source": "go_bp", "group_id": "GO:2", "group_name": "tight", "term_size": 20,
         "covered_genes": "a, b, c, e", "z_score": -7.0, "p_fdr": 0.01},
    ])
    depth = {"GO:1": 3, "GO:2": 6}
    ancestors = {"GO:1": set(), "GO:2": set()}
    cfg = _cfg(tmp_path, scope="pooled", merge_dag_lineage=False, force_representatives=["GO:1"])
    out = deduplicate(table, cfg, depth, ancestors)
    rep = out[out["is_representative"]]
    assert list(rep["group_id"]) == ["GO:1"]
    assert rep["representative_source"].iloc[0] == "forced"


def test_deduplicate_exactly_one_representative_per_cluster(tmp_path):
    """Every cluster has exactly one representative row."""
    table = _combined([
        {"source": "go_bp", "group_id": f"GO:{i}", "group_name": f"g{i}", "term_size": 10,
         "covered_genes": "a, b, c" if i < 3 else "x, y, z", "z_score": -float(i), "p_fdr": 0.01 * (i + 1)}
        for i in range(6)
    ])
    depth = {f"GO:{i}": 5 for i in range(6)}
    ancestors = {f"GO:{i}": set() for i in range(6)}
    out = deduplicate(table, _cfg(tmp_path, scope="pooled", merge_dag_lineage=False), depth, ancestors)
    per_cluster = out.groupby("redundancy_cluster")["is_representative"].sum()
    assert (per_cluster == 1).all()


def test_deduplicate_per_source_scope_keeps_sources_separate(tmp_path):
    """per_source scope never merges identical-member terms across sources."""
    table = _combined([
        {"source": "go_cc", "group_id": "GO:X", "group_name": "SSU", "term_size": 41,
         "covered_genes": "a, b, c", "z_score": -5.86, "p_fdr": 0.016},
        {"source": "go_macrocomplex", "group_id": "GO:X", "group_name": "SSU", "term_size": 41,
         "covered_genes": "a, b, c", "z_score": -5.86, "p_fdr": 0.031},
    ])
    depth = {"GO:X": 4}
    ancestors = {"GO:X": set()}
    out = deduplicate(table, _cfg(tmp_path, scope="per_source", merge_dag_lineage=False), depth, ancestors)
    # Same GO:X in two sources -> two clusters under per_source, both representatives.
    assert out["redundancy_cluster"].nunique() == 2
    assert out["is_representative"].sum() == 2
