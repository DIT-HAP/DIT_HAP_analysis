# Remaining Notebooks Migration Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Port 11 remaining source notebooks from DIT_HAP_pipeline/workflow/notebooks/ into DIT_HAP_analysis as deterministic Snakemake rules + scripts (Batch A/B), split rules + human notebooks (Batch C), and pure human notebooks (Batch D).

**Architecture:** Each deterministic notebook becomes a workflow/rules/*.smk file + workflow/scripts/{stage}/*.py scripts, following the existing split-rule pattern (prepare → compute → finalize). Human-judgment steps remain as notebooks with explicit input/output contracts in their first cell. Design doc: docs/plans/2026-07-19-remaining-notebooks-migration-design.md.

**Tech Stack:** Python 3.12, Snakemake 9.13+, pandas, numpy, scipy, matplotlib, seaborn, loguru, goatools (for complex GO analysis), altair (for interactive notebooks only).

---

## Task 1: Spike-in analysis (spikein.smk)

**Files:**
- Create: `workflow/rules/spikein.smk`
- Create: `workflow/scripts/spikein/run_spikein_analysis.py`
- Create: `tests/test_spikein.py`

**Steps (TDD, 2-5 min each):**

1. Write the failing test in `tests/test_spikein.py`.
2. Run: `python -m pytest tests/test_spikein.py -v` → verify all tests fail (ImportError or AssertionError).
3. Write minimal implementation in `workflow/scripts/spikein/run_spikein_analysis.py` and `workflow/src/spikein/` if needed.
4. Run: `python -m pytest tests/test_spikein.py -v` → verify all tests pass.
5. Write `workflow/rules/spikein.smk`.
6. Run: `snakemake -n results/spikein/spike_in_stats.tsv` → verify dry-run succeeds.
7. Commit: `git commit -m "feat(spikein): add spike-in analysis rule + script + tests"`

**Test file `tests/test_spikein.py`:**

```python
"""Tests for spike-in analysis core computations."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.scripts.spikein.run_spikein_analysis import (
    assign_ratio_by_order,
    build_spike_sites_df,
    compute_linear_regression_stats,
    SPIKE_IN_RATIO,
)


SPIKE_IN_RATIO_EXPECTED = np.array([1.5, 4, 16, 64, 256, 1024]) / 100000


def test_spike_in_ratio_constant():
    """SPIKE_IN_RATIO constant matches known dilution series."""
    np.testing.assert_allclose(SPIKE_IN_RATIO, SPIKE_IN_RATIO_EXPECTED)


def test_assign_ratio_by_order_basic():
    """Reads rank ascending → lowest read = lowest ratio; relative values are log2-normalised."""
    sub = pd.DataFrame(
        {"Reads": [100.0, 200.0, 400.0, 800.0, 1600.0, 3200.0]},
        index=range(6),
    )
    result = assign_ratio_by_order(sub.copy(), SPIKE_IN_RATIO_EXPECTED)
    # Ratio assigned by rank (0-indexed ascending)
    np.testing.assert_allclose(result["Ratio"].values, SPIKE_IN_RATIO_EXPECTED)
    # Minimum read subtracted: lowest becomes 0
    assert result["Reads"].min() == 0.0
    # Relative_Dilution_Ratio: log2(ratio / max_ratio)
    expected_rel_dil = np.log2(SPIKE_IN_RATIO_EXPECTED / SPIKE_IN_RATIO_EXPECTED.max())
    np.testing.assert_allclose(result["Relative_Dilution_Ratio"].values, expected_rel_dil)


def test_assign_ratio_by_order_monotone_reads():
    """With perfectly ordered reads, the rank assignment is identity."""
    reads = np.array([10.0, 40.0, 160.0, 640.0, 2560.0, 10240.0])
    sub = pd.DataFrame({"Reads": reads}, index=range(6))
    result = assign_ratio_by_order(sub.copy(), SPIKE_IN_RATIO_EXPECTED)
    np.testing.assert_allclose(result["Ratio"].values, SPIKE_IN_RATIO_EXPECTED)


def test_build_spike_sites_df_shape():
    """build_spike_sites_df returns one row per spike-in site with expected columns."""
    mock_index = pd.MultiIndex.from_tuples(
        [("I", 3749394, "-"), ("II", 3344505, "-"), ("II", 185161, "-"),
         ("II", 1157130, "-"), ("II", 3065244, "-")],
        names=["Chr", "Coordinate", "Strand"],
    )
    reads = np.array([10.0, 20.0, 40.0, 80.0, 160.0, 320.0])
    mock_df = pd.DataFrame(
        {("Sample1", "0h"): reads[0], ("Sample2", "0h"): reads[1]},
        index=mock_index,
    )
    mock_df = pd.DataFrame(
        [[reads] for _ in range(5)],
        index=mock_index,
        columns=pd.MultiIndex.from_tuples([("S1",), ("S2",), ("S3",), ("S4",), ("S5",), ("S6",)]),
    )
    spike_in_sites = {
        "DY215": {"chr": "I", "coord": 3749394, "strand": "-"},
        "DY217": {"chr": "II", "coord": 3344505, "strand": "-"},
        "DY218": {"chr": "II", "coord": 185161, "strand": "-"},
        "DY339": {"chr": "II", "coord": 1157130, "strand": "-"},
        "DY348": {"chr": "II", "coord": 3065244, "strand": "-"},
    }
    df = build_spike_sites_df(mock_df, spike_in_sites)
    assert len(df) == 5
    assert "Strain" in df.columns


def test_compute_linear_regression_stats():
    """compute_linear_regression_stats returns slope, r_value, p_value, r2."""
    x = np.array([-10.0, -8.0, -6.0, -4.0, -2.0, 0.0])
    y = 0.95 * x + 0.1  # near-perfect linear
    stats = compute_linear_regression_stats(pd.Series(x), pd.Series(y))
    assert abs(stats["slope"] - 0.95) < 0.01
    assert stats["r2"] > 0.99
    assert "p_value" in stats
```

**Key logic to port from `spike_in.ipynb`:**
- 5 spike-in coordinates (DY215/DY217/DY218/DY339/DY348) → `config/analysis.yaml` under `spikein.coordinates`
- Extract reads by chr+coord+strand from `raw_reads.filtered.tsv` (MultiIndex `[Chr, Coordinate, Strand]`)
- `SPIKE_IN_RATIO = np.array([1.5, 4, 16, 64, 256, 1024]) / 100000`
- `assign_ratio_by_order`: rank reads ascending → assign ratios in order; subtract min read; compute `Relative_Read_Ratio = log2(reads / max_reads)`; compute `Relative_Dilution_Ratio = log2(ratio / max_ratio)` — applied per Strain group via `groupby("Strain").apply(...)`
- `linregress(Relative_Dilution_Ratio, Relative_Read_Ratio)` over all sites combined
- Output: `spike_in_stats.tsv` (long-form per site/sample) + `spike_in_correlation.pdf`

**Rule skeleton `workflow/rules/spikein.smk`:**

```python
# =============================================================================
# spikein.smk — Spike-in dilution linearity QC
# =============================================================================
#
# Standalone (no {dataset} wildcard): reads the Spikein project's filtered
# insertion table, extracts 5 known spike-in coordinates, fits a log-log linear
# regression of reads vs known dilution ratios, and emits a stats TSV + PDF.
#
# Single rule (no prepare/compute split — data is tiny and self-contained).

rule run_spikein_analysis:
    input:
        raw_reads=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets']['Spikein']['release_dir']}/insertion_level/raw_reads.filtered.tsv"
        ),
    output:
        stats="results/spikein/spike_in_stats.tsv",
        figure="results/spikein/spike_in_correlation.pdf",
    params:
        spike_in_sites=config.get("spikein", {}).get("coordinates", {}),
    log:
        "logs/spikein/run_spikein_analysis.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [spikein] Running spike-in linearity QC..."
    shell:
        """
        python workflow/scripts/spikein/run_spikein_analysis.py \
            --raw-reads {input.raw_reads} \
            --output-stats {output.stats} \
            --output-figure {output.figure} \
            --spike-in-sites '{params.spike_in_sites}' &> {log}
        """
```

**`config/analysis.yaml` additions for Task 1:**

```yaml
spikein:
  coordinates:
    DY215: {chr: "I",  coord: 3749394, strand: "-"}
    DY217: {chr: "II", coord: 3344505, strand: "-"}
    DY218: {chr: "II", coord: 185161,  strand: "-"}
    DY339: {chr: "II", coord: 1157130, strand: "-"}
    DY348: {chr: "II", coord: 3065244, strand: "-"}
```

---

## Task 2: Gene coverage analysis (coverage.smk)

**Files:**
- Create: `workflow/rules/coverage.smk`
- Create: `workflow/scripts/coverage/compute_coverage_stats.py`
- Create: `tests/test_coverage.py`

**Steps:**

1. Write failing tests in `tests/test_coverage.py`.
2. Run: `python -m pytest tests/test_coverage.py -v` → verify fail.
3. Write `workflow/scripts/coverage/compute_coverage_stats.py` (port logic from `gene_coverage_analysis.ipynb`).
4. Run: `python -m pytest tests/test_coverage.py -v` → verify pass.
5. Write `workflow/rules/coverage.smk`.
6. Run: `snakemake -n "results/coverage/{_DATASET}/coverage_stats.tsv"` → verify dry-run.
7. Commit: `git commit -m "feat(coverage): add gene coverage analysis rule + script + tests"`

**Test file `tests/test_coverage.py`:**

```python
"""Tests for gene coverage computation logic."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.scripts.coverage.compute_coverage_stats import (
    IN_GENE_FILTER,
    compute_insertion_coverage,
    compute_gene_coverage,
    compute_essentiality_coverage,
)


def _make_insertion_annotation(n_in_gene=30, n_intergenic=10):
    """Synthetic insertion annotation table with required columns."""
    rows = []
    for i in range(n_in_gene):
        rows.append({"Type": "Coding exon", "Distance_to_stop_codon": 10})
    for i in range(n_intergenic):
        rows.append({"Type": "Intergenic region", "Distance_to_stop_codon": 0})
    # Edge: in-gene but too close to stop codon
    rows.append({"Type": "Coding exon", "Distance_to_stop_codon": 3})
    idx = pd.MultiIndex.from_tuples(
        [(f"I", i * 100, "+", f"g{i}") for i in range(len(rows))],
        names=["Chr", "Coordinate", "Strand", "Gene"],
    )
    return pd.DataFrame(rows, index=idx)


def test_in_gene_filter_constant():
    """Exact filter string is preserved from source notebook (quirk)."""
    assert IN_GENE_FILTER == "Type != 'Intergenic region' and Distance_to_stop_codon > 4"


def test_compute_insertion_coverage_counts():
    """In-gene count = rows passing filter; intergenic = complement."""
    annotation = _make_insertion_annotation(n_in_gene=30, n_intergenic=10)
    # 30 in-gene with Distance_to_stop_codon=10, 1 edge with Distance=3 (fails), 10 intergenic
    result = compute_insertion_coverage(annotation)
    assert result["total"] == 41
    assert result["in_gene"] == 30  # edge case excluded
    assert result["intergenic"] == 11


def test_compute_gene_coverage_counts():
    """covered = DR not NaN; not_covered = DR is NaN."""
    gene_result = pd.DataFrame({
        "Systematic ID": ["g1", "g2", "g3", "g4"],
        "DR": [0.5, None, 0.8, None],
        "DeletionLibrary_essentiality": ["E", "V", "E", "V"],
    })
    result = compute_gene_coverage(gene_result)
    assert result["total"] == 4
    assert result["covered"] == 2
    assert result["not_covered"] == 2


def test_compute_essentiality_coverage_essential():
    """Essential (E) gene coverage split is correct."""
    gene_result = pd.DataFrame({
        "Systematic ID": [f"g{i}" for i in range(6)],
        "DR": [0.5, None, 0.8, None, 0.7, None],
        "DeletionLibrary_essentiality": ["E", "E", "E", "V", "V", "V"],
    })
    result = compute_essentiality_coverage(gene_result)
    assert result["essential"]["total"] == 3
    assert result["essential"]["covered"] == 2
    assert result["non_essential"]["total"] == 3
    assert result["non_essential"]["covered"] == 2
```

**Key logic from `gene_coverage_analysis.ipynb`:**
- `IN_GENE_FILTER = "Type != 'Intergenic region' and Distance_to_stop_codon > 4"` (exact string, verbatim quirk)
- Inputs: `insertion_level/fitting_results.tsv` (MultiIndex) + `insertion_level/annotations.tsv` (same MultiIndex)
- Count total vs in-gene insertions using the filter; intergenic = total - in-gene
- Per-chromosome fractions from `coverage_meta` (insertion density analysis intermediate — recompute from annotations directly rather than loading the upstream report file)
- Gene coverage: `DR.notna()` as covered; split by `DeletionLibrary_essentiality == 'E'` / `'V'`
- Donut chart using `workflow/src/plotting/generic.py::donut_chart`
- DR/DL histograms split by essentiality category (3 rows × 2 cols subplot)
- Output: `coverage_stats.tsv` + `coverage_figures.pdf`

**Rule skeleton `workflow/rules/coverage.smk`:**

```python
# =============================================================================
# coverage.smk — Gene insertion coverage statistics
# =============================================================================
#
# Per-dataset: computes insertion coverage fractions (in-gene vs intergenic)
# and gene coverage (covered vs not covered) per essentiality class.
# Uses the exact IN_GENE_FILTER string from the source notebook (quirk).

rule compute_coverage_stats:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/insertion_level/fitting_results.tsv"
        ),
        annotations=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/insertion_level/annotations.tsv"
        ),
        gene_level=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        essentiality_verification_csv="resources/curated/essentiality_verification.csv",
    output:
        stats="results/coverage/{dataset}/coverage_stats.tsv",
        figures="results/coverage/{dataset}/coverage_figures.pdf",
    log:
        "logs/coverage/compute_coverage_stats_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [coverage] Computing insertion + gene coverage for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/coverage/compute_coverage_stats.py \
            --fitting-results {input.fitting_results} \
            --annotations {input.annotations} \
            --gene-level {input.gene_level} \
            --essentiality-verification-csv {input.essentiality_verification_csv} \
            --output-stats {output.stats} \
            --output-figures {output.figures} &> {log}
        """
```

---

## Task 3: Deletion library verification (verification.smk)

**Files:**
- Create: `workflow/rules/verification.smk`
- Create: `workflow/scripts/verification/compare_deletion_library.py`
- Create: `tests/test_verification.py`

**Steps:**

1. Write failing tests in `tests/test_verification.py`.
2. Run: `python -m pytest tests/test_verification.py -v` → verify fail.
3. Write `workflow/scripts/verification/compare_deletion_library.py` and add color maps to `workflow/src/plotting/style.py`.
4. Run: `python -m pytest tests/test_verification.py -v` → verify pass.
5. Write `workflow/rules/verification.smk`.
6. Run: `snakemake -n "results/verification/{_DATASET}/verification_stats.tsv"` → verify dry-run.
7. Commit: `git commit -m "feat(verification): add deletion library verification rule + script + tests"`

**Test file `tests/test_verification.py`:**

```python
"""Tests for deletion library verification logic."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.scripts.verification.compare_deletion_library import (
    merge_deletion_library,
    compute_category_stats,
    CATEGORY_COLOR_MAP,
    DONUT_COLOR_MAP,
)


def _make_gene_results():
    return pd.DataFrame({
        "Systematic ID": ["g1", "g2", "g3", "g4", "g5"],
        "DR": [0.9, 0.1, 0.5, 0.8, 0.2],
        "DL": [5.0, 1.0, 3.0, 4.0, 2.0],
        "DeletionLibrary_essentiality": ["E", "V", "E", "V", "E"],
    })


def _make_deletion_library():
    return pd.DataFrame({
        "Updated_Systematic_ID": ["g1", "g2", "g3", "g4", "g5"],
        "Category": ["spores", "WT", "germinated", "small colonies", "E"],
    })


def test_category_color_map_has_required_keys():
    """CATEGORY_COLOR_MAP must cover all expected deletion phenotype labels."""
    required = {"WT", "small colonies", "very small colonies", "E",
                "E (tiny colonies)", "microcolonies", "germinated", "spores", "Not verified"}
    assert required.issubset(set(CATEGORY_COLOR_MAP.keys()))


def test_donut_color_map_has_required_keys():
    """DONUT_COLOR_MAP must cover all expected donut chart categories."""
    required = {"spores", "germinated", "microcolonies", "E",
                "E (tiny colonies)", "very small colonies", "small colonies", "WT"}
    assert required.issubset(set(DONUT_COLOR_MAP.keys()))


def test_merge_deletion_library_joins_on_systematic_id():
    """merge_deletion_library joins on Systematic ID / Updated_Systematic_ID."""
    gene = _make_gene_results()
    dl = _make_deletion_library()
    merged = merge_deletion_library(gene, dl)
    assert "Category" in merged.columns
    assert len(merged) == 5


def test_compute_category_stats_returns_counts():
    """compute_category_stats returns count per category in the merged frame."""
    gene = _make_gene_results()
    dl = _make_deletion_library()
    merged = merge_deletion_library(gene, dl)
    stats = compute_category_stats(merged)
    assert "category" in stats.columns
    assert "count" in stats.columns
    assert stats["count"].sum() == 5


def test_category_with_essentiality_flag():
    """small colonies + E essentiality → 'small colonies (E)' label."""
    from workflow.scripts.verification.compare_deletion_library import apply_category_with_essentiality
    row_sc_e = pd.Series({"Category": "small colonies", "DeletionLibrary_essentiality": "E"})
    row_sc_v = pd.Series({"Category": "small colonies", "DeletionLibrary_essentiality": "V"})
    assert apply_category_with_essentiality(row_sc_e) == "small colonies (E)"
    assert apply_category_with_essentiality(row_sc_v) == "small colonies"
```

**Key logic from `compare_with_deletion_library.ipynb`:**
- `CATEGORY_COLOR_MAP` and `DONUT_COLOR_MAP` use `COLORS` indices 0, 2, 4, 5, 7, -1, -4 → add both dicts to `workflow/src/plotting/style.py` as module-level constants (import from there in script)
- Merge gene-level fitting results with `deletion_library_categories.xlsx` (`Updated_Systematic_ID` → `Systematic ID`) via `DeletionLibrary_essentiality` column
- `apply_category_with_essentiality`: if `Category == 'small colonies'` and `DeletionLibrary_essentiality == 'E'` → `"small colonies (E)"` else `Category`
- Merge with `essentiality_verification.csv` (simplified: map `"E,small colonies"` → `"E"`, `"E,WT"` → `"E"`, `"Leu-condition"` → `"E"`)
- Produce donut chart per category + DR scatter plot
- Output: `verification_stats.tsv` + `deletion_library_comparison.pdf`

**Rule skeleton `workflow/rules/verification.smk`:**

```python
# =============================================================================
# verification.smk — Deletion library phenotype verification
# =============================================================================
#
# Per-dataset: merges gene-level DIT-HAP results with Hayles-2013 deletion
# library categories and the curated essentiality verification table.
# Produces donut charts per phenotype category + DR scatter.

rule compare_deletion_library:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
        deletion_library="resources/curated/deletion_library_categories.xlsx",
        essentiality_verification="resources/curated/essentiality_verification.csv",
    output:
        stats="results/verification/{dataset}/verification_stats.tsv",
        figures="results/verification/{dataset}/deletion_library_comparison.pdf",
    log:
        "logs/verification/compare_deletion_library_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [verification] Comparing deletion library for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/verification/compare_deletion_library.py \
            --fitting-results {input.fitting_results} \
            --deletion-library {input.deletion_library} \
            --essentiality-verification {input.essentiality_verification} \
            --output-stats {output.stats} \
            --output-figures {output.figures} &> {log}
        """
```

---

## Task 4: Non-coding RNA analysis (noncoding_rna.smk)

**Files:**
- Create: `workflow/rules/noncoding_rna.smk`
- Create: `workflow/scripts/noncoding_rna/analyze_noncoding_rna.py`
- Create: `tests/test_noncoding_rna.py`

**Steps:**

1. Write failing tests in `tests/test_noncoding_rna.py`.
2. Run: `python -m pytest tests/test_noncoding_rna.py -v` → verify fail.
3. Write `workflow/scripts/noncoding_rna/analyze_noncoding_rna.py`.
4. Run: `python -m pytest tests/test_noncoding_rna.py -v` → verify pass.
5. Write `workflow/rules/noncoding_rna.smk`.
6. Run: `snakemake -n "results/noncoding_rna/{_DATASET}/ncrna_stats.tsv"` → verify dry-run.
7. Commit: `git commit -m "feat(noncoding_rna): add ncRNA analysis rule + script + tests"`

**Test file `tests/test_noncoding_rna.py`:**

```python
"""Tests for non-coding RNA analysis core computations."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import numpy as np
import pytest

from workflow.scripts.noncoding_rna.analyze_noncoding_rna import (
    normalize_chromosome_names,
    merge_gtrnadb_by_position,
    extract_tRNA_amino_acid_and_anticodon,
    compute_tRNA_copy_number,
)


def test_normalize_chromosome_names_replaces_chr_prefix():
    """chrI/II/III → I/II/III (exact quirk from source notebook)."""
    df = pd.DataFrame({"#Chr": ["chrI", "chrII", "chrIII", "mitochondrial"]})
    result = normalize_chromosome_names(df)
    assert list(result["#Chr"]) == ["I", "II", "III", "mitochondrial"]


def test_merge_gtrnadb_by_position_uses_chr_start_end():
    """GtRNAdb merge is on #Chr+Start+End (not Name) — key quirk."""
    ncrna = pd.DataFrame({
        "#Chr": ["I", "II"],
        "Start": [100, 200],
        "End": [200, 300],
        "Systematic ID": ["g1", "g2"],
    })
    gtrnadb = pd.DataFrame({
        "#Chr": ["I", "III"],
        "Start": [100, 500],
        "End": [200, 600],
        "GtRNAdb_Name": ["tRNA-Ala-AGC-1-1", "tRNA-Gly-GCC-2-1"],
    })
    merged = merge_gtrnadb_by_position(ncrna, gtrnadb)
    assert len(merged) == 2
    assert merged.loc[merged["Systematic ID"] == "g1", "GtRNAdb_Name"].iloc[0] == "tRNA-Ala-AGC-1-1"
    assert pd.isna(merged.loc[merged["Systematic ID"] == "g2", "GtRNAdb_Name"].iloc[0])


def test_extract_tRNA_amino_acid_and_anticodon():
    """Amino acid and anticodon parsed from GtRNAdb_Name field."""
    row_with = pd.Series({
        "Systematic ID": "SPATRNA.Ala1",
        "GtRNAdb_Name": "Schpo_chr1.trna1-AlaAGC",
    })
    # Amino acid from sysID (TRNA<AA>.)
    row_sys = pd.Series({
        "Systematic ID": "SPTRNAALA.01",
        "GtRNAdb_Name": "schiPomb_972H-tRNA-Ala-AGC-1-1",
    })
    result = extract_tRNA_amino_acid_and_anticodon(row_sys)
    assert result["Anticodon"] == "AGC"


def test_compute_tRNA_copy_number():
    """tRNA_copy_number = count of tRNAs sharing the same Amino_Acid+Anticodon."""
    df = pd.DataFrame({
        "Amino_Acid": ["Ala", "Ala", "Gly", "Ala"],
        "Anticodon": ["AGC", "AGC", "GCC", "AGC"],
        "Systematic ID": ["t1", "t2", "t3", "t4"],
    })
    result = compute_tRNA_copy_number(df)
    assert list(result.loc[result["Systematic ID"].isin(["t1", "t2", "t4"]), "tRNA_copy_number"]) == [3, 3, 3]
    assert result.loc[result["Systematic ID"] == "t3", "tRNA_copy_number"].iloc[0] == 1
```

**Key logic from `non_coding_RNA_analysis.ipynb`:**
- `normalize_chromosome_names`: `.replace({"chrI": "I", "chrII": "II", "chrIII": "III"})` on `#Chr` column in GtRNAdb
- `merge_gtrnadb_by_position`: merge on `["#Chr", "Start", "End"]` (not Name) — left join ncRNA meta onto GtRNAdb
- Input ncRNA fitting results path: confirm via `data_config`; expected upstream path pattern `{release}/gene_level/noncoding_rna_fitting_results.tsv`
- Marguerat 2012 Excel: `resources/literature/margueratQuantitativeAnalysisFission2012.xlsx`, sheet `Table_S2`, columns `MM1.tot.cpc_ex / MM2.tot.cpc_ex / MN1.tot.cpc_ex / MN2.tot.cpc_ex`; index on `Systematic.name`; produce `mean` per condition (EMM_Proliferating / EMM_Nitrogen_Starved)
- `extract_tRNA_amino_acid_and_anticodon`: amino acid from `re.search(r"TRNA(\w+)\.", sysID)`, anticodon from `GtRNAdb_Name.split("-")[2]`
- Filter nuclear tRNAs: `Feature == 'tRNA' and #Chr != 'mitochondrial'`
- Output: `ncrna_stats.tsv` + `ncrna_analysis.pdf`

**Rule skeleton `workflow/rules/noncoding_rna.smk`:**

```python
# =============================================================================
# noncoding_rna.smk — Non-coding RNA depletion analysis
# =============================================================================
#
# Per-dataset: merges ncRNA gene-level stats with GtRNAdb tRNA annotations
# (matched by chr+start+end, not name) and Marguerat 2012 mRNA abundance.
# Analyzes tRNA vs other ncRNA depletion patterns.

rule analyze_noncoding_rna:
    input:
        ncrna_fitting=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/noncoding_rna_fitting_results.tsv"
        ),
        ncrna_bed="resources/pombase/{pombase_version}/genome_region/non_coding_rna.bed",
        gtrnadb_bed="resources/pombase/schiPomb_972H-tRNAs.bed",
        marguerat_excel="resources/literature/margueratQuantitativeAnalysisFission2012.xlsx",
    output:
        stats="results/noncoding_rna/{dataset}/ncrna_stats.tsv",
        figures="results/noncoding_rna/{dataset}/ncrna_analysis.pdf",
    params:
        pombase_version=lambda wc: DATASETS["reference"]["pombase_version"],
    log:
        "logs/noncoding_rna/analyze_noncoding_rna_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [noncoding_rna] Analyzing ncRNA depletion for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/noncoding_rna/analyze_noncoding_rna.py \
            --ncrna-fitting {input.ncrna_fitting} \
            --ncrna-bed {input.ncrna_bed} \
            --gtrnadb-bed {input.gtrnadb_bed} \
            --marguerat-excel {input.marguerat_excel} \
            --output-stats {output.stats} \
            --output-figures {output.figures} &> {log}
        """
```

---

## Task 5: Large-scale study comparison (comparison.smk)

**Files:**
- Create: `workflow/rules/comparison.smk`
- Create: `workflow/scripts/comparison/compare_large_scale_studies.py`
- Create: `tests/test_comparison.py`

**Steps:**

1. Write failing tests in `tests/test_comparison.py`.
2. Run: `python -m pytest tests/test_comparison.py -v` → verify fail.
3. Write `workflow/scripts/comparison/compare_large_scale_studies.py`.
4. Run: `python -m pytest tests/test_comparison.py -v` → verify pass.
5. Write `workflow/rules/comparison.smk`.
6. Add `comparison.gRNA_data_file` to `config/analysis.yaml`.
7. Run: `snakemake -n "results/comparison/{_DATASET}/fitness_correlation_stats.tsv"` → verify dry-run.
8. Commit: `git commit -m "feat(comparison): add large-scale study comparison rule + script + tests"`

**Test file `tests/test_comparison.py`:**

```python
"""Tests for large-scale study comparison core computations."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.scripts.comparison.compare_large_scale_studies import (
    clip_density_columns,
    compute_pearson_r,
    CLIP_UPPER,
    DENSITY_COLUMNS,
)


def test_clip_upper_constant():
    """clip(upper=200) is the exact value from source notebook."""
    assert CLIP_UPPER == 200


def test_density_columns_include_required_names():
    """Integration density, ipkm, uipkm columns must be clipped."""
    required = {"Integration density, in-vivo (integrations/kb/million inserts)", "ipkm", "uipkm"}
    assert required.issubset(set(DENSITY_COLUMNS))


def test_clip_density_columns_caps_at_200():
    """Values above 200 are clipped to exactly 200."""
    df = pd.DataFrame({
        "Integration density, in-vivo (integrations/kb/million inserts)": [50.0, 250.0, 200.0],
        "ipkm": [100.0, 300.0, 199.0],
        "uipkm": [10.0, 201.0, 5.0],
        "other_col": [1000.0, 2000.0, 3000.0],
    })
    result = clip_density_columns(df)
    assert result["Integration density, in-vivo (integrations/kb/million inserts)"].max() == 200.0
    assert result["ipkm"].max() == 200.0
    assert result["uipkm"].max() == 200.0
    # other_col untouched
    assert result["other_col"].max() == 3000.0


def test_compute_pearson_r_returns_r_and_pvalue():
    """compute_pearson_r returns (r, p_value) for two numeric series."""
    x = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0])
    y = pd.Series([1.1, 1.9, 3.1, 3.9, 5.1])
    r, p = compute_pearson_r(x, y)
    assert abs(r - 1.0) < 0.05
    assert p < 0.05


def test_compute_pearson_r_ignores_nan_pairs():
    """NaN in either column → drop pair before correlation."""
    x = pd.Series([1.0, 2.0, np.nan, 4.0])
    y = pd.Series([1.0, 2.0, 3.0, np.nan])
    r, p = compute_pearson_r(x, y)
    # Only (1,1) and (2,2) survive — perfect positive correlation
    assert abs(r - 1.0) < 0.01
```

**Key logic from `compare_with_other_large_scale_studies.ipynb`:**
- `CLIP_UPPER = 200`; apply `.clip(upper=200)` to `"Integration density, in-vivo (integrations/kb/million inserts)"`, `"ipkm"`, `"uipkm"` before any plotting
- Merge DIT-HAP cluster data with gRNA fitted params (on `Systematic ID`); merge with `pombe_coding_gene_protein_features.tsv`
- Pairwise scatter: `um_DIT_HAP` vs `um_gRNA` vs `Barseq_from_dulab` vs `Barseq_from_koch` vs `Integration density` vs `colony_size_Malecki2016` etc.
- KDE overlay: `scipy.stats.gaussian_kde` on each pair; contour lines on scatter
- Pearson r + p-value per pair; annotate on scatter subplot
- `gRNA_data_file` registered in `config/analysis.yaml` → script receives path via arg
- Output: `fitness_correlation_stats.tsv` + `pairwise_fitness_comparison.pdf`

**`config/analysis.yaml` additions for Task 5:**

```yaml
comparison:
  gRNA_data_file: "resources/curated/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv"
  clip_upper: 200
```

**Rule skeleton `workflow/rules/comparison.smk`:**

```python
# =============================================================================
# comparison.smk — Pairwise fitness comparison with other large-scale studies
# =============================================================================
#
# Batch B (requires resources/curated/final_clusters.tsv).
# Per-dataset: merges DIT-HAP data with gRNA, Barseq, integration density,
# colony size. Pairwise scatter with KDE overlay + Pearson r statistics.

rule compare_large_scale_studies:
    input:
        final_clusters="resources/curated/final_clusters.tsv",
        protein_features=lambda wc: (
            f"results/features/{DATASETS['reference']['pombase_version']}/"
            "pombe_coding_gene_protein_features.tsv"
        ),
        gRNA_data=config.get("comparison", {}).get(
            "gRNA_data_file",
            "resources/curated/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv"
        ),
    output:
        stats="results/comparison/{dataset}/fitness_correlation_stats.tsv",
        figures="results/comparison/{dataset}/pairwise_fitness_comparison.pdf",
    params:
        clip_upper=config.get("comparison", {}).get("clip_upper", 200),
    log:
        "logs/comparison/compare_large_scale_studies_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [comparison] Running pairwise fitness comparison for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/comparison/compare_large_scale_studies.py \
            --final-clusters {input.final_clusters} \
            --protein-features {input.protein_features} \
            --grna-data {input.gRNA_data} \
            --clip-upper {params.clip_upper} \
            --output-stats {output.stats} \
            --output-figures {output.figures} &> {log}
        """
```

---

## Task 6: Complex coherence analysis (complex.smk)

**Files:**
- Create: `workflow/rules/complex.smk`
- Create: `workflow/scripts/complex/analyze_complex_modules.py`
- Create: `workflow/scripts/complex/compute_complex_coherence.py`
- Create: `workflow/src/complex/__init__.py`
- Create: `workflow/src/complex/coherence.py`
- Create: `tests/test_complex_coherence.py`

**Steps:**

1. Write failing tests in `tests/test_complex_coherence.py`.
2. Run: `python -m pytest tests/test_complex_coherence.py -v` → verify fail.
3. Write `workflow/src/complex/coherence.py` (geometric_median, coherence_metrics, compute_distance_zscore).
4. Write `workflow/scripts/complex/compute_complex_coherence.py` and `analyze_complex_modules.py`.
5. Run: `python -m pytest tests/test_complex_coherence.py -v` → verify pass.
6. Write `workflow/rules/complex.smk`.
7. Add complex config section to `config/analysis.yaml`.
8. Run: `snakemake -n "results/complex/{_DATASET}/complex_coherence_metrics.tsv"` → verify dry-run.
9. Commit: `git commit -m "feat(complex): add complex coherence analysis rule + src + tests"`

**Test file `tests/test_complex_coherence.py`:**

```python
"""Tests for complex coherence algorithm (Weiszfeld geometric median + permutation test)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.src.complex.coherence import (
    geometric_median,
    coherence_metrics,
    compute_distance_zscore,
    EPSILON,
)


def test_epsilon_guard_value():
    """Zero-distance epsilon is exactly 1e-5 (source notebook quirk)."""
    assert EPSILON == 1e-5


def test_geometric_median_single_point():
    """geometric_median of a single point is that point itself."""
    points = np.array([[1.0, 2.0]])
    gm = geometric_median(points)
    np.testing.assert_allclose(gm, [1.0, 2.0], atol=1e-6)


def test_geometric_median_collinear_symmetric():
    """geometric_median of symmetric points converges to centroid."""
    points = np.array([[-1.0, 0.0], [1.0, 0.0], [0.0, 0.0]])
    gm = geometric_median(points)
    np.testing.assert_allclose(gm, [0.0, 0.0], atol=1e-4)


def test_geometric_median_uses_component_wise_median_init():
    """Initialization from component-wise median (not mean)."""
    # With points that have outliers, median init is more robust than mean init.
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 0.0], [100.0, 0.0]])
    gm = geometric_median(points)
    # Should converge near the cluster of 3, not pulled to outlier
    assert gm[0] < 2.0


def test_coherence_metrics_returns_required_keys():
    """coherence_metrics returns centroid_x, centroid_y + 6 distance stats."""
    points = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 0.866], [0.5, -0.866]])
    result = coherence_metrics(points)
    required = {"centroid_x", "centroid_y", "median_distance", "mean_distance",
                "std_distance", "min_distance", "max_distance", "mpd"}
    assert required.issubset(set(result.keys()))


def test_compute_distance_zscore_returns_zscore_and_pvalue():
    """compute_distance_zscore returns observed_mpd, z_score, p_value, n_permutations."""
    rng = np.random.default_rng(42)
    all_points = rng.standard_normal((100, 2))
    complex_indices = list(range(5))  # tight cluster
    result = compute_distance_zscore(all_points, complex_indices, n_permutations=100, random_state=42)
    assert "observed_mpd" in result
    assert "z_score" in result
    assert "p_value" in result
    assert result["n_permutations"] == 100


def test_compute_distance_zscore_tight_cluster_has_low_zscore():
    """A tight cluster should have a lower MPD than random (negative z-score)."""
    rng = np.random.default_rng(42)
    # Background: spread out
    all_points = rng.standard_normal((200, 2)) * 5
    # Tight complex: 8 points near origin
    all_points[:8] = rng.standard_normal((8, 2)) * 0.1
    result = compute_distance_zscore(all_points, list(range(8)), n_permutations=500, random_state=42)
    assert result["z_score"] < 0  # observed MPD < permutation mean
```

**Key logic from `complex_analysis.ipynb` sections 4-5:**
- **Section 4 (module visualization)**: named module dict (cytoplasmic translation, kinetochore, mitochondria, vesicle, vacuolar ATPase) → `config/analysis.yaml` under `complex.modules`; use `workflow/src/plotting/gene_level.py` visualization functions
- **Section 5 (Weiszfeld geometric median)**:
  - Init: `median = np.median(points, axis=0)` (component-wise)
  - Iteration: `distances = np.linalg.norm(points - median, axis=1)` → clip: `distances = np.where(distances < EPSILON, EPSILON, distances)`; `weights = 1.0 / distances`; `median = np.sum(weights[:, None] * points, axis=0) / weights.sum()`; converge when `np.linalg.norm(new - old) < 1e-7`
  - `coherence_metrics()`: computes centroid (= geometric median) x/y, all pairwise L2 distances, MPD = `np.median(pairwise_distances)`
  - `compute_distance_zscore()`: for `n_permutations=1000`, sample `size` random genes from background, compute MPD; z-score = `(observed - mean(null)) / std(null)`
  - Group filter: complexes with `3 <= size <= 300` and `median DR of members > 0.3`
- Output: `complex_coherence_metrics.tsv` + `complex_module_visualization.pdf` + `coherence_analysis.pdf`

**`config/analysis.yaml` additions for Task 6:**

```yaml
complex:
  min_complex_size: 3
  max_complex_size: 300
  dr_threshold: 0.3
  n_permutations: 1000
  random_state: 42
  modules:
    cytoplasmic_translation: []  # filled from PomBase macromolecular complex annotation
    kinetochore: []
    mitochondria: []
    vesicle: []
    vacuolar_ATPase: []
```

**Rule skeleton `workflow/rules/complex.smk`:**

```python
# =============================================================================
# complex.smk — Macromolecular complex coherence analysis
# =============================================================================
#
# Batch B (requires resources/curated/final_clusters.tsv).
# Two-part: module visualization (named complexes) + coherence permutation test
# (all complexes 3<=size<=300 with DR>0.3).
# Weiszfeld geometric median in workflow/src/complex/coherence.py.

_CPLX_WORK = "results/complex/{dataset}/_work"

rule analyze_complex_modules:
    input:
        final_clusters="resources/curated/final_clusters.tsv",
        complex_annotation="resources/pombase/{pombase_version}/macromolecular_complex.tsv",
    output:
        module_viz=f"{_CPLX_WORK}/module_visualization_done.flag",
        module_figure="results/complex/{dataset}/complex_module_visualization.pdf",
    params:
        pombase_version=lambda wc: DATASETS["reference"]["pombase_version"],
        modules=config.get("complex", {}).get("modules", {}),
    log:
        "logs/complex/analyze_complex_modules_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [complex] Visualizing named complex modules for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/complex/analyze_complex_modules.py \
            --final-clusters {input.final_clusters} \
            --complex-annotation {input.complex_annotation} \
            --modules '{params.modules}' \
            --output-flag {output.module_viz} \
            --output-figure {output.module_figure} &> {log}
        """

rule compute_complex_coherence:
    input:
        final_clusters="resources/curated/final_clusters.tsv",
        complex_annotation="resources/pombase/{pombase_version}/macromolecular_complex.tsv",
    output:
        metrics="results/complex/{dataset}/complex_coherence_metrics.tsv",
        coherence_figure="results/complex/{dataset}/coherence_analysis.pdf",
    params:
        pombase_version=lambda wc: DATASETS["reference"]["pombase_version"],
        min_size=config.get("complex", {}).get("min_complex_size", 3),
        max_size=config.get("complex", {}).get("max_complex_size", 300),
        dr_threshold=config.get("complex", ).get("dr_threshold", 0.3),
        n_permutations=config.get("complex", {}).get("n_permutations", 1000),
        random_state=config.get("complex", {}).get("random_state", 42),
    log:
        "logs/complex/compute_complex_coherence_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [complex] Computing coherence metrics for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/complex/compute_complex_coherence.py \
            --final-clusters {input.final_clusters} \
            --complex-annotation {input.complex_annotation} \
            --min-size {params.min_size} \
            --max-size {params.max_size} \
            --dr-threshold {params.dr_threshold} \
            --n-permutations {params.n_permutations} \
            --random-state {params.random_state} \
            --output-metrics {output.metrics} \
            --output-figure {output.coherence_figure} &> {log}
        """
```

---

## Task 7: UTR insertion analysis — deterministic part (utr.smk)

**Files:**
- Create: `workflow/rules/utr.smk`
- Create: `workflow/scripts/utr/classify_utr_insertions.py`
- Create: `tests/test_utr.py`

**Steps:**

1. Write failing tests in `tests/test_utr.py`.
2. Run: `python -m pytest tests/test_utr.py -v` → verify fail.
3. Write `workflow/scripts/utr/classify_utr_insertions.py`.
4. Run: `python -m pytest tests/test_utr.py -v` → verify pass.
5. Write `workflow/rules/utr.smk`.
6. Run: `snakemake -n "results/utr/{_DATASET}/utr_insertion_stats.tsv"` → verify dry-run.
7. Commit: `git commit -m "feat(utr): add UTR insertion classification rule + script + tests"`

**Test file `tests/test_utr.py`:**

```python
"""Tests for UTR insertion classification (assign_UTR_type)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import pytest

from workflow.scripts.utr.classify_utr_insertions import (
    assign_UTR_type,
    UTR_DISTANCE_THRESHOLD,
    filter_intergenic_near_gene,
)


def test_utr_distance_threshold_constant():
    """distance_threshold=400 bp is the exact value from the source notebook."""
    assert UTR_DISTANCE_THRESHOLD == 400


def _row(strand, dist_start, dist_end, gene_strand="+"):
    return pd.Series({
        "Strand": strand,
        "gene_strand": gene_strand,
        "Distance_to_region_start": dist_start,
        "Distance_to_region_end": dist_end,
    })


def test_assign_UTR_type_5utr_forward_strand():
    """Insertion near gene start (forward strand) → 5UTR."""
    row = _row(strand="+", dist_start=50, dist_end=600, gene_strand="+")
    assert assign_UTR_type(row) == "5UTR"


def test_assign_UTR_type_3utr_forward_strand():
    """Insertion near gene end (forward strand) → 3UTR."""
    row = _row(strand="+", dist_start=600, dist_end=80, gene_strand="+")
    assert assign_UTR_type(row) == "3UTR"


def test_assign_UTR_type_5utr_reverse_strand():
    """Insertion near gene END (reverse strand) = 5' end → 5UTR."""
    row = _row(strand="-", dist_start=600, dist_end=100, gene_strand="-")
    assert assign_UTR_type(row) == "5UTR"


def test_assign_UTR_type_3utr_reverse_strand():
    """Insertion near gene START (reverse strand) = 3' end → 3UTR."""
    row = _row(strand="-", dist_start=50, dist_end=600, gene_strand="-")
    assert assign_UTR_type(row) == "3UTR"


def test_assign_UTR_type_ambiguous_both_close():
    """Insertion within threshold of both boundaries → '5UTR' (start takes priority)."""
    row = _row(strand="+", dist_start=100, dist_end=100, gene_strand="+")
    result = assign_UTR_type(row)
    assert result in ("5UTR", "3UTR", "ambiguous")  # implementation choice


def test_filter_intergenic_near_gene():
    """Filter: Type == Intergenic region AND (dist_start < 400 OR dist_end < 400)."""
    df = pd.DataFrame({
        "Type": ["Intergenic region", "Intergenic region", "Coding exon", "Intergenic region"],
        "Distance_to_region_start": [100, 500, 50, 500],
        "Distance_to_region_end": [600, 600, 600, 350],
    })
    result = filter_intergenic_near_gene(df)
    # Row 0: intergenic, dist_start < 400 → pass
    # Row 1: intergenic, both >= 400 → fail
    # Row 2: not intergenic → fail
    # Row 3: intergenic, dist_end < 400 → pass
    assert len(result) == 2
    assert 0 in result.index
    assert 3 in result.index
```

**Key logic from `upstream_and_downstream_analysis.ipynb`:**
- `UTR_DISTANCE_THRESHOLD = 400` (bp)
- `filter_intergenic_near_gene`: `"Type == 'Intergenic region' and (Distance_to_region_start < 400 or Distance_to_region_end < 400)"`
- `assign_UTR_type(row)`:
  - If `gene_strand == "+"`: `Distance_to_region_start < threshold` → `"5UTR"`; `Distance_to_region_end < threshold` → `"3UTR"`
  - If `gene_strand == "-"`: `Distance_to_region_end < threshold` → `"5UTR"`; `Distance_to_region_start < threshold` → `"3UTR"`
  - (strand-aware because region_start/end are genomic coordinates, not transcript-relative)
- Merge insertion stats + gene-level stats; compute `um_ratio` (insertion DR / gene median DR) and `A_ratio` (insertion A / gene median A)
- Output: `utr_insertion_stats.tsv` (per-insertion, with UTR type column)

**Rule skeleton `workflow/rules/utr.smk`:**

```python
# =============================================================================
# utr.smk — UTR insertion classification (deterministic part of Batch C)
# =============================================================================
#
# Per-dataset: classifies intergenic insertions near gene boundaries as 5UTR or
# 3UTR (strand-aware, distance_threshold=400bp). Merges with insertion + gene
# stats, computes um_ratio and A_ratio.
# The human-review notebook is notebooks/domain_analysis/review_utr_insertions.ipynb.

rule classify_utr_insertions:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/insertion_level/fitting_results.tsv"
        ),
        annotations=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/insertion_level/annotations.tsv"
        ),
        gene_level=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
    output:
        stats="results/utr/{dataset}/utr_insertion_stats.tsv",
    params:
        distance_threshold=400,
    log:
        "logs/utr/classify_utr_insertions_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [utr] Classifying UTR insertions for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/utr/classify_utr_insertions.py \
            --fitting-results {input.fitting_results} \
            --annotations {input.annotations} \
            --gene-level {input.gene_level} \
            --distance-threshold {params.distance_threshold} \
            --output-stats {output.stats} &> {log}
        """
```

---

## Task 8: Domain differences — deterministic part (domain_differences.smk)

**Files:**
- Create: `workflow/rules/domain_differences.smk`
- Create: `workflow/scripts/domain_differences/compute_domain_stats.py`
- Create: `tests/test_domain_differences.py`

**Steps:**

1. Write failing tests in `tests/test_domain_differences.py`.
2. Run: `python -m pytest tests/test_domain_differences.py -v` → verify fail.
3. Write `workflow/scripts/domain_differences/compute_domain_stats.py`.
4. Run: `python -m pytest tests/test_domain_differences.py -v` → verify pass.
5. Write `workflow/rules/domain_differences.smk`.
6. Run: `snakemake -n "results/domain_differences/{_DATASET}/domain_candidate_stats.tsv"` → verify dry-run.
7. Commit: `git commit -m "feat(domain_differences): add domain stats rule + script + tests"`

**Test file `tests/test_domain_differences.py`:**

```python
"""Tests for domain difference (intra-gene DR heterogeneity) computation."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import pytest

from workflow.scripts.domain_differences.compute_domain_stats import (
    DR_THRESHOLD,
    compute_insertion_fraction,
    filter_high_dr_genes,
    compute_domain_candidate_stats,
)


def test_dr_threshold_constant():
    """DR > 0.15 threshold is the exact value from source notebook."""
    assert DR_THRESHOLD == 0.15


def _make_insertion_data():
    """Synthetic in-gene insertions for 2 genes."""
    rows = []
    # Gene g1: 10 insertions distributed across the gene
    for i in range(10):
        rows.append({
            "Systematic ID": "g1",
            "Coordinate": i * 100,
            "Distance_to_start_codon": i * 50,
            "Distance_to_stop_codon": (9 - i) * 50,
            "DR": 0.6 + i * 0.01,
        })
    # Gene g2: 5 insertions (DR below threshold)
    for i in range(5):
        rows.append({
            "Systematic ID": "g2",
            "Coordinate": i * 200,
            "Distance_to_start_codon": i * 100,
            "Distance_to_stop_codon": (4 - i) * 100,
            "DR": 0.05,
        })
    return pd.DataFrame(rows)


def _make_gene_stats():
    return pd.DataFrame({
        "Systematic ID": ["g1", "g2", "g3"],
        "DR": [0.6, 0.05, 0.3],
        "DL": [5.0, 1.0, 3.0],
    })


def test_filter_high_dr_genes_applies_threshold():
    """filter_high_dr_genes keeps only genes with gene-level DR > 0.15."""
    gene_stats = _make_gene_stats()
    result = filter_high_dr_genes(gene_stats)
    assert set(result["Systematic ID"]) == {"g1", "g3"}
    assert "g2" not in result["Systematic ID"].values


def test_compute_insertion_fraction_range():
    """Insertion fraction is in [0, 1] for all insertions."""
    insertions = _make_insertion_data()
    gene_stats = _make_gene_stats()
    # Only g1 passes DR threshold
    high_dr = filter_high_dr_genes(gene_stats)
    result = compute_insertion_fraction(insertions, high_dr)
    assert (result["insertion_fraction"] >= 0).all()
    assert (result["insertion_fraction"] <= 1).all()


def test_compute_domain_candidate_stats_output_columns():
    """Output has required columns for downstream domain review notebook."""
    insertions = _make_insertion_data()
    gene_stats = _make_gene_stats()
    result = compute_domain_candidate_stats(insertions, gene_stats)
    required = {"Systematic ID", "n_insertions", "mean_insertion_fraction",
                "std_insertion_fraction", "gene_DR"}
    assert required.issubset(set(result.columns))
```

**Key logic from `genes_with_domain_differences.ipynb`:**
- `DR_THRESHOLD = 0.15`: select in-gene insertions for genes with gene-level `DR > 0.15`
- `compute_insertion_fraction`: position each insertion relative to gene length as `insertion_fraction = Distance_to_start_codon / (Distance_to_start_codon + Distance_to_stop_codon)`; clamp to [0, 1]
- `compute_domain_candidate_stats`: per gene, compute `n_insertions`, `mean_insertion_fraction`, `std_insertion_fraction`, `gene_DR`; sort by `std_insertion_fraction` descending (candidates with most heterogeneous insertion distribution)
- Inputs: insertion-level `fitting_results.tsv` (MultiIndex) + `annotations.tsv` + gene-level `fitting_results.tsv`
- Apply `IN_GENE_FILTER` (same as coverage.smk) to select in-gene insertions
- Output: `domain_candidate_stats.tsv`

**Rule skeleton `workflow/rules/domain_differences.smk`:**

```python
# =============================================================================
# domain_differences.smk — Intra-gene DR heterogeneity candidates (Batch C)
# =============================================================================
#
# Per-dataset: selects in-gene insertions for genes with DR>0.15, computes
# insertion fraction statistics (position relative to start/stop codon).
# Human review in notebooks/domain_analysis/review_domain_differences.ipynb.

rule compute_domain_stats:
    input:
        fitting_results=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/insertion_level/fitting_results.tsv"
        ),
        annotations=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/insertion_level/annotations.tsv"
        ),
        gene_level=lambda wc: (
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets'][wc.dataset]['release_dir']}/gene_level/fitting_results.tsv"
        ),
    output:
        stats="results/domain_differences/{dataset}/domain_candidate_stats.tsv",
    params:
        dr_threshold=0.15,
    log:
        "logs/domain_differences/compute_domain_stats_{dataset}.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [domain_differences] Computing domain candidate stats for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/domain_differences/compute_domain_stats.py \
            --fitting-results {input.fitting_results} \
            --annotations {input.annotations} \
            --gene-level {input.gene_level} \
            --dr-threshold {params.dr_threshold} \
            --output-stats {output.stats} &> {log}
        """
```

---

## Task 9: Human notebooks (Batch C + D)

**Files:**
- Create: `notebooks/domain_analysis/review_utr_insertions.ipynb`
- Create: `notebooks/domain_analysis/review_domain_differences.ipynb`
- Create: `notebooks/enrichment/further_analysis.ipynb`
- Create: `notebooks/thesis_figures/figures.ipynb`
- Create: `notebooks/reference/visualize_dit_hap_style.ipynb`

**Steps:**

1. Create each notebook with the mandatory first and second cells.
2. Adapt logic from the source notebook (replace hardcoded paths with config dataclass).
3. Verify each notebook opens and the config cell executes without error: `jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=60 notebooks/<path> --output /tmp/test_nb.ipynb` (run on a machine with data available, or skip with `@pytest.mark.integration`).
4. Commit: `git commit -m "feat(notebooks): add human review + analysis notebooks (Batch C+D)"`

**Notebook structure (all five notebooks follow this template):**

Cell 1 (Markdown):
```markdown
## Inputs

- `results/{stage}/{dataset}/...` — (list actual files)
- `resources/curated/...` — (curated inputs if any)

## Outputs

- `reports/{stage}/...` — (PDFs / Excel files produced)
```

Cell 2 (Code — Config dataclass):
```python
from dataclasses import dataclass
from pathlib import Path
import pandas as pd

DATASET = "LD_DIT_HAP"  # change as needed

@dataclass
class Config:
    dataset: str = DATASET
    results_dir: Path = Path(f"../../results")
    resources_dir: Path = Path("../../resources")

    def __post_init__(self):
        # Paths resolved here — fail fast if files are missing
        self.stats_file = self.results_dir / f"<stage>/{self.dataset}/<file>.tsv"
        # ... other paths

cfg = Config()
```

**`notebooks/domain_analysis/review_utr_insertions.ipynb` specifics:**
- Inputs: `results/utr/{dataset}/utr_insertion_stats.tsv` + `resources/curated/non_essential_domain_candidates.xlsx`
- Key logic: load UTR stats, filter by gene DR threshold, plot intragenic insertion visualizations using `workflow/src/plotting/gene_level.py::intragenic_insertion_visualization` for candidate genes
- Config: add `workflow/src/` to `sys.path` via `sys.path.insert(0, str(Path("../../workflow/src")))`

**`notebooks/domain_analysis/review_domain_differences.ipynb` specifics:**
- Inputs: `results/domain_differences/{dataset}/domain_candidate_stats.tsv` + `resources/curated/non_essential_domain_candidates.xlsx`
- Key logic: rank-order plot of `std_insertion_fraction`; histogram of insertion positions for curated gene list; intragenic insertion visualization for top candidates
- Sort candidates by `std_insertion_fraction` descending

**`notebooks/enrichment/further_analysis.ipynb` specifics:**
- Inputs: `results/enrichment/raw/{dataset}/{pombase_version}/go_enrichment_full_filtered.tsv` + `resources/curated/enrichment_term_categorization.xlsx` + `resources/curated/final_clusters.tsv`
- Key logic (from `further_analysis_based_on_enrichment.ipynb`):
  - `get_paralog_cluster(x, genes2cluster)`: split paralog string by comma, map each gene to cluster id (leave as gene name if not in dict), join with comma
  - `distribution_bar_for_given_genes()` from `workflow/src/plotting/gene_level.py::distribution_bar_for_given_genes`
  - `donut_chart()` from `workflow/src/plotting/generic.py`
  - Cluster-level enrichment bar charts; paralog cluster analysis using `gene_features["paralogs_of_genes"]`
- Config property `go_enrichment_results` loads the categorized GO enrichment Excel

**`notebooks/thesis_figures/figures.ipynb` specifics:**
- Inputs: all `results/{stage}/{dataset}/` outputs (clustering, enrichment, comparison, complex, verification, coverage)
- Output: `reports/thesis/*.pdf` (manuscript-quality figures)
- Config: `output_dir = Path("../../reports/thesis"); output_dir.mkdir(parents=True, exist_ok=True)`
- Note: this is a summary notebook — it assembles figures from results of upstream stages; no new computation

**`notebooks/reference/visualize_dit_hap_style.ipynb` specifics:**
- Copy from `DIT_HAP_pipeline/workflow/notebooks/visualize_dit_hap_style.ipynb` (copy as-is)
- Update imports: replace `sys.path.append("../src")` with `sys.path.insert(0, str(Path("../../workflow/src")))`
- Replace any direct path references with `workflow/src/plotting/style.py`
- First cell: add standard Inputs/Outputs contract markdown

---

## Task 10: Snakefile integration + config updates

**Files:**
- Modify: `Snakefile` (add 8 new includes)
- Modify: `config/analysis.yaml` (add spikein, comparison, complex sections)

**Steps:**

1. Add 8 includes to `Snakefile` (after existing includes block).
2. Add commented targets to `rule all` for each new stage.
3. Add config sections to `config/analysis.yaml`: `spikein.coordinates`, `comparison.gRNA_data_file`, `complex.*`
4. Run: `snakemake -n results/spikein/spike_in_stats.tsv` → verify dry-run (Batch A target).
5. Run: `snakemake -n "results/coverage/{_DATASET}/coverage_stats.tsv"` → verify dry-run.
6. Run: `snakemake -n "results/verification/{_DATASET}/verification_stats.tsv"` → verify dry-run.
7. Run: `snakemake -n "results/noncoding_rna/{_DATASET}/ncrna_stats.tsv"` → verify dry-run.
8. Run: `snakemake -n "results/comparison/{_DATASET}/fitness_correlation_stats.tsv"` → verify dry-run (Batch B).
9. Run: `snakemake -n "results/complex/{_DATASET}/complex_coherence_metrics.tsv"` → verify dry-run.
10. Run: `snakemake -n "results/utr/{_DATASET}/utr_insertion_stats.tsv"` → verify dry-run (Batch C).
11. Run: `snakemake -n "results/domain_differences/{_DATASET}/domain_candidate_stats.tsv"` → verify dry-run.
12. Commit: `git commit -m "feat(snakefile): integrate 8 new analysis stages into Snakefile + config"`

**Snakefile additions (after line 35, after `include: "workflow/rules/pcr_qc.smk"`):**

```python
# --- Batch A: no final_clusters.tsv dependency ---
include: "workflow/rules/spikein.smk"
include: "workflow/rules/coverage.smk"
include: "workflow/rules/verification.smk"
include: "workflow/rules/noncoding_rna.smk"
# --- Batch B: requires resources/curated/final_clusters.tsv ---
include: "workflow/rules/comparison.smk"
include: "workflow/rules/complex.smk"
# --- Batch C: split deterministic + human notebooks ---
include: "workflow/rules/utr.smk"
include: "workflow/rules/domain_differences.smk"
```

**`rule all` additions (commented, after existing commented targets):**

```python
        # Batch A (no final_clusters.tsv dependency):
        # "results/spikein/spike_in_stats.tsv",
        # f"results/coverage/{_DATASET}/coverage_stats.tsv",
        # f"results/verification/{_DATASET}/verification_stats.tsv",
        # f"results/noncoding_rna/{_DATASET}/ncrna_stats.tsv",
        # Batch B (requires resources/curated/final_clusters.tsv):
        # f"results/comparison/{_DATASET}/fitness_correlation_stats.tsv",
        # f"results/complex/{_DATASET}/complex_coherence_metrics.tsv",
        # Batch C (requires insertion-level results):
        # f"results/utr/{_DATASET}/utr_insertion_stats.tsv",
        # f"results/domain_differences/{_DATASET}/domain_candidate_stats.tsv",
```

**Complete `config/analysis.yaml` additions:**

```yaml
# --- Spike-in QC (spikein.smk) ---
spikein:
  coordinates:
    DY215: {chr: "I",  coord: 3749394, strand: "-"}
    DY217: {chr: "II", coord: 3344505, strand: "-"}
    DY218: {chr: "II", coord: 185161,  strand: "-"}
    DY339: {chr: "II", coord: 1157130, strand: "-"}
    DY348: {chr: "II", coord: 3065244, strand: "-"}

# --- Large-scale comparison (comparison.smk) ---
comparison:
  gRNA_data_file: "resources/curated/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv"
  clip_upper: 200

# --- Complex coherence (complex.smk) ---
complex:
  min_complex_size: 3
  max_complex_size: 300
  dr_threshold: 0.3
  n_permutations: 1000
  random_state: 42
  modules:
    cytoplasmic_translation: []
    kinetochore: []
    mitochondria: []
    vesicle: []
    vacuolar_ATPase: []
```

**Expected dry-run output for each target (example for spikein):**

```
$ snakemake -n results/spikein/spike_in_stats.tsv
Building DAG of jobs...
Job stats:
job                     count
--------------------  -------
run_spikein_analysis        1
total                       1

1 of 1 steps (100%) done
This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

---

## Implementation order summary

```
Task 1  → spikein.smk + test_spikein.py            (Batch A, standalone)
Task 2  → coverage.smk + test_coverage.py           (Batch A, per-dataset)
Task 3  → verification.smk + test_verification.py  (Batch A, per-dataset)
Task 4  → noncoding_rna.smk + test_noncoding_rna.py (Batch A, per-dataset)
Task 5  → comparison.smk + test_comparison.py       (Batch B, needs final_clusters)
Task 6  → complex.smk + src/complex/ + tests        (Batch B, needs final_clusters)
Task 7  → utr.smk + test_utr.py                    (Batch C, deterministic part)
Task 8  → domain_differences.smk + tests           (Batch C, deterministic part)
Task 9  → 5 human notebooks                        (Batch C+D, no new tests)
Task 10 → Snakefile + config/analysis.yaml          (integration)
```

All Tasks 1-8 follow strict TDD: write failing test → implement → verify pass → dry-run → commit.
Tasks 9-10 do not have unit tests (notebooks) or are integration-level (Snakefile wiring).
