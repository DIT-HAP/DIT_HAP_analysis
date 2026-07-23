import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "coherence"))

import pandas as pd
import pytest


def _write_macrocomplex(tmp_path: Path) -> Path:
    df = pd.DataFrame({
        "complex_term_id": ["GO:0001", "GO:0001", "GO:0002"],
        "GO_term_name": ["alpha complex", "alpha complex", "beta complex"],
        "systematic_id": ["SPAC1", "SPAC2", "SPBC1"],
        "symbol": ["gene1", "gene2", "gene3"],
    })
    d = tmp_path / "ontologies_and_associations"
    d.mkdir(parents=True)
    df.to_csv(d / "macromolecular_complex_annotation.tsv", sep="\t", index=False)
    return tmp_path


def test_prepare_macrocomplex_returns_long_table(tmp_path):
    from prepare_annotation import prepare
    out = prepare("go_macrocomplex", _write_macrocomplex(tmp_path))
    assert list(out.columns) == ["source", "group_id", "group_name",
                                 "Systematic ID", "Name", "n_group_genes"]
    assert len(out) == 3


def test_prepare_unknown_source_raises(tmp_path):
    from prepare_annotation import prepare
    with pytest.raises(ValueError, match="unknown source"):
        prepare("kegg_nope", tmp_path)
