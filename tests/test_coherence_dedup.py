import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import itertools

import pandas as pd
import pytest

from deduplicate_terms import (
    member_set,
    jaccard_index,
    candidate_pairs,
    build_clusters,
    deduplicate,
    group_member_table,
    DedupConfig,
)
from workflow.src.coherence.fractions import moonlighting_fraction


# --- primitives -------------------------------------------------------------
def test_jaccard_index_does_not_reward_mere_containment():
    """A subset nested in a superset scores |child|/|parent|, NOT 1.0.

    This is the property that motivated dropping the overlap coefficient, whose
    min() denominator scored ANY contained term 1.0. GO propagation makes
    containment ubiquitous, so that connected the whole redundancy graph: 2,097
    terms collapsed into 2 clusters once every source spelled genes with the same
    identifier. Pinned here so it cannot regress unnoticed.
    """
    parent = {"a", "b", "c", "d"}
    assert jaccard_index(parent, {"a", "b"}) == pytest.approx(0.5)
    assert jaccard_index(parent, {"a"}) == pytest.approx(0.25)
    assert jaccard_index(parent, parent) == 1.0
    assert jaccard_index(set(), {"a"}) == 0.0


def test_jaccard_index_partial():
    assert jaccard_index({"a", "b", "c"}, {"b", "c", "d"}) == pytest.approx(2 / 4)


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


def test_member_set_rejects_comma_joined_string():
    """scored_member_names is a LIST column now; a joined string would split into characters.

    The old TSV intermediate stored "a, b, c" and parsing it back was the reason
    for the join/parse round-trip. Now that the column is a real list, a string
    reaching here is a bug, and iterating it would silently produce {'a', 'b',
    ',', ' '} — so it raises instead.
    """
    assert member_set(["a", "b"]) == {"a", "b"}
    assert member_set(None) == set()
    assert member_set(float("nan")) == set()
    with pytest.raises(TypeError, match="list column"):
        member_set("a, b, c")


# --- clustering -------------------------------------------------------------
def _sub(rows):
    """rows: list of (group_id, scored_member_names list). Minimal cols for build_clusters."""
    return pd.DataFrame(
        [{"group_id": gid, "scored_member_names": sorted(genes)} for gid, genes in rows]
    )


def test_build_clusters_is_transitive():
    """A-B and B-C merges leave A, B, C in one cluster even though A and C never match."""
    sub = _sub([
        ("GO:1", {"a", "b", "c", "d"}),
        ("GO:2", {"a", "b", "c", "e"}),   # 3/4 with GO:1
        ("GO:3", {"a", "b", "c", "f"}),   # 3/4 with both, disjoint-ish otherwise
        ("GO:4", {"x", "y", "z"}),        # disjoint from all
    ])
    labels = build_clusters(sub, threshold=0.5, merge_dag_lineage=False, ancestors={})
    assert labels[0] == labels[1] == labels[2]
    assert labels[3] != labels[0]
    # dense labels numbered by first appearance, which the output sort keys on
    assert sorted(labels) == [0, 0, 0, 1]


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
    """rows: dicts with source, group_id, group_name, n_scored_members, scored_member_names,
    median_pairwise_distance_z, q_value."""
    return pd.DataFrame(rows)


def _cfg(tmp_path, **kw):
    return DedupConfig(
        combined=tmp_path / "c.tsv", obo=tmp_path / "o.obo",
        output_all=tmp_path / "all.tsv", output_representatives=tmp_path / "rep.tsv",
        output_group_members=tmp_path / "members.tsv",
        **kw,
    )


def test_deduplicate_picks_best_qvalue_representative(tmp_path):
    """Within a redundant cluster the min-q_value term is the representative; all rows kept."""
    table = _combined([
        {"source": "go_bp", "group_id": "GO:1", "group_name": "big", "n_scored_members": 200,
         "scored_member_names": ["a", "b", "c", "d"], "median_pairwise_distance_z": -5.0, "q_value": 0.05},
        {"source": "go_bp", "group_id": "GO:2", "group_name": "tight", "n_scored_members": 20,
         "scored_member_names": ["a", "b", "c", "e"], "median_pairwise_distance_z": -7.0, "q_value": 0.01},  # best q
        {"source": "go_cc", "group_id": "GO:9", "group_name": "other", "n_scored_members": 5,
         "scored_member_names": ["x", "y", "z"], "median_pairwise_distance_z": -3.0, "q_value": 0.2},
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
    # The collapsed terms ride on every row of the cluster, the representative's included.
    assert set(cl["non_representative_terms"]) == {"go_bp:big (GO:1)"}
    assert out.loc[out["group_id"] == "GO:9", "non_representative_terms"].iloc[0] == ""


def test_deduplicate_force_representative_overrides(tmp_path):
    """A forced group_id becomes its cluster's representative even with a worse q."""
    table = _combined([
        {"source": "go_bp", "group_id": "GO:1", "group_name": "big", "n_scored_members": 200,
         "scored_member_names": ["a", "b", "c", "d"], "median_pairwise_distance_z": -5.0, "q_value": 0.05},
        {"source": "go_bp", "group_id": "GO:2", "group_name": "tight", "n_scored_members": 20,
         "scored_member_names": ["a", "b", "c", "e"], "median_pairwise_distance_z": -7.0, "q_value": 0.01},
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
        {"source": "go_bp", "group_id": f"GO:{i}", "group_name": f"g{i}", "n_scored_members": 10,
         "scored_member_names": ["a", "b", "c"] if i < 3 else ["x", "y", "z"],
         "median_pairwise_distance_z": -float(i), "q_value": 0.01 * (i + 1)}
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
        {"source": "go_cc", "group_id": "GO:X", "group_name": "SSU", "n_scored_members": 41,
         "scored_member_names": ["a", "b", "c"], "median_pairwise_distance_z": -5.86, "q_value": 0.016},
        {"source": "go_macrocomplex", "group_id": "GO:X", "group_name": "SSU", "n_scored_members": 41,
         "scored_member_names": ["a", "b", "c"], "median_pairwise_distance_z": -5.86, "q_value": 0.031},
    ])
    depth = {"GO:X": 4}
    ancestors = {"GO:X": set()}
    out = deduplicate(table, _cfg(tmp_path, scope="per_source", merge_dag_lineage=False), depth, ancestors)
    # Same GO:X in two sources -> two clusters under per_source, both representatives.
    assert out["redundancy_cluster"].nunique() == 2
    assert out["is_representative"].sum() == 2


# --- linkage: single (connected components) vs complete (hierarchical cut) ---
def _chain_sub():
    """A~B and B~C clear 0.5; A~C is far below it. The classic chaining case."""
    return _sub([
        ("GO:A", {"a", "b", "c", "d", "e"}),              # J(A,B)=4/7=0.57
        ("GO:B", {"a", "b", "c", "d", "f", "g"}),         # J(B,C)=4/8=0.50
        ("GO:C", {"c", "d", "f", "g", "h", "i"}),         # J(A,C)=2/9=0.22
    ])


def test_single_linkage_chains_the_whole_line():
    """Default behaviour, pinned because it is what `complete` is being compared to."""
    labels = build_clusters(_chain_sub(), threshold=0.5, merge_dag_lineage=False, ancestors={})
    assert labels[0] == labels[1] == labels[2]


def test_complete_linkage_breaks_the_chain():
    """A and C are not redundant with each other, so no cluster may hold both."""
    labels = build_clusters(
        _chain_sub(), threshold=0.5, merge_dag_lineage=False, ancestors={},
        linkage_method="complete",
    )
    assert labels[0] == labels[1]          # A~B still merge
    assert labels[2] != labels[0]          # C does not ride in on B
    assert sorted(labels) == [0, 0, 1]     # first-appearance numbering preserved


def test_complete_linkage_holds_at_every_threshold():
    """The property that motivates it: no cluster ever contains a pair below threshold.

    Checked exhaustively rather than on one pair, because "the widest internal pair
    bounds the merge height" is the whole claim.
    """
    sub = _sub([
        ("GO:C1", {"a", "b", "c", "d", "e", "f"}),
        ("GO:C2", {"a", "b", "c", "d", "e", "g"}),
        ("GO:C3", {"d", "e", "f", "g", "h", "i"}),
        ("GO:C4", {"h", "i", "j", "k", "l", "m"}),
        ("GO:C5", {"a", "b", "h", "i", "j", "k"}),
        ("GO:C6", {"x", "y", "z", "w", "v", "u"}),      # disjoint from everything
    ])
    members = [set(row) for row in sub["scored_member_names"]]
    for threshold in (0.3, 0.5, 0.6, 0.8):
        labels = build_clusters(
            sub, threshold=threshold, merge_dag_lineage=False, ancestors={},
            linkage_method="complete",
        )
        for i, j in itertools.combinations(range(len(members)), 2):
            if labels[i] != labels[j] or not (members[i] | members[j]):
                continue
            overlap = len(members[i] & members[j]) / len(members[i] | members[j])
            assert overlap >= threshold, (
                f"thr={threshold}: terms {i},{j} share a cluster at Jaccard {overlap:.2f}"
            )


def test_complete_linkage_keeps_a_disjoint_term_alone():
    labels = build_clusters(
        _chain_sub(), threshold=0.5, merge_dag_lineage=False, ancestors={},
        linkage_method="complete",
    )
    assert len(set(labels)) == 2


def _existing_cfg(tmp_path, **kw):
    """A config whose declared inputs exist, so validate() reaches the param checks."""
    for name in ("c.tsv", "o.obo"):
        (tmp_path / name).touch()
    return _cfg(tmp_path, **kw)


@pytest.mark.parametrize("method", ["average", "complete"])
def test_lineage_rule_is_rejected_outside_single_linkage(tmp_path, method):
    """A DAG-ancestor edge is not a similarity, so it cannot ride on a linkage cut.

    Silently dropping it (or silently ignoring the linkage) would answer a different
    question than the config asked, so the config refuses the pair outright.
    """
    cfg = _existing_cfg(tmp_path, linkage=method, merge_dag_lineage=True)
    with pytest.raises(ValueError, match="merge_dag_lineage needs linkage='single'"):
        cfg.validate()


def test_unknown_linkage_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="linkage must be one of"):
        _existing_cfg(tmp_path, linkage="ward").validate()


def test_single_linkage_still_carries_the_lineage_edge(tmp_path):
    """The one combination that is allowed to set both must validate."""
    _existing_cfg(tmp_path, linkage="single", merge_dag_lineage=True).validate()


def test_defaults_are_self_consistent(tmp_path):
    """A bare DedupConfig() must validate, and its linkage must be `complete`.

    `complete` and the DAG-lineage rule are mutually exclusive (a similarity cannot
    express an ancestor edge), so a lineage default of True would make the
    out-of-the-box config — and therefore the CLI's own --help — advertise a
    combination that raises on the first call. The two constants have to move
    together; this pins the pair.
    """
    from deduplicate_terms import DEFAULT_LINKAGE, DEFAULT_MERGE_DAG_LINEAGE

    assert DEFAULT_LINKAGE == "complete"
    assert DEFAULT_MERGE_DAG_LINEAGE is False
    _existing_cfg(tmp_path).validate()   # no kwargs: exactly the defaults


# --- group-member long table + moonlighting ---------------------------------
def _annotated():
    """Two clusters + one singleton; `c` is in all three, `a` in two."""
    return pd.DataFrame({
        "redundancy_cluster": ["all:1", "all:1", "all:2", "all:3"],
        "is_representative": [True, False, True, True],
        "source": ["go_bp", "go_bp", "go_cc", "go_cc"],
        "group_id": ["GO:1", "GO:2", "GO:3", "GO:4"],
        "group_name": ["rep one", "alias", "second", "third"],
        "scored_member_names": [["a", "b"], ["c"], ["c", "d"], ["a", "c"]],
    })


def test_group_member_table_unions_the_cluster_and_names_it_by_the_representative():
    long, _mode = group_member_table(_annotated())

    assert list(long.columns) == ["source", "group_id", "group_name", "full_name",
                                  "gene", "redundancy_cluster", "n_groups", "is_moonlighting"]
    assert len(long) == 7  # {a,b,c} + {c,d} + {a,c}
    # The group's genes are the UNION over its terms: `c` comes from the alias and
    # still counts as a GO:1 gene, once.
    assert set(long[long["group_id"] == "GO:1"]["gene"]) == {"a", "b", "c"}
    rep = long[long["group_id"] == "GO:1"].iloc[0]
    assert rep["group_name"] == "rep one"
    assert rep["full_name"] == "rep one\nalias"  # representative first
    assert rep["source"] == "go_bp"


def test_moonlighting_flags_genes_above_the_views_cut_and_shares_are_per_term():
    long, cut = group_member_table(_annotated())
    per_gene = long.drop_duplicates("gene").set_index("gene")

    # Group counts: a in {all:1, all:3}, b in {all:1}, c in all three, d in {all:2}.
    assert per_gene["n_groups"].to_dict() == {"a": 2, "b": 1, "c": 3, "d": 1}
    # The cut is the 75th percentile of {2,1,3,1} = 2.25, so only c (3) is above it
    # — strictly above, which is what keeps a gene sitting ON the cut unflagged.
    assert cut == pytest.approx(2.25)
    assert per_gene["is_moonlighting"].to_dict() == {"a": False, "b": False, "c": True, "d": False}

    # Per term, over its own members, with `c` the only flagged gene: GO:1 {a,b} -> 0,
    # GO:2 {c} -> 1, GO:3 {c,d} -> 1/2, GO:4 {a,c} -> 1/2. The mapping run() does,
    # over the shared helper — the same one the per-view figures use.
    flagged = set(long.loc[long["is_moonlighting"], "gene"])
    shares = [moonlighting_fraction(members, flagged)
              for members in _annotated()["scored_member_names"]]
    assert shares == [0.0, 1.0, 0.5, 0.5]
