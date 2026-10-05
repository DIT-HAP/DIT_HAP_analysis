# =============================================================================
# Snakefile — DIT-HAP analysis entry point
# =============================================================================

from snakemake.utils import min_version
from pathlib import Path
import yaml

min_version("9.0")

workdir: "/data/c/yangyusheng_optimized/DIT_HAP_analysis"

# This project's analysis parameters (clustering k, enrichment thresholds, ml
# splits, ...). Unlike datasets.yaml (a data registry read directly below), these
# ARE experiment parameters, so they flow through Snakemake's `config` object.
configfile: "config/analysis.yaml"

# ---------------------------------------------------------------------------
# Dataset registry (not a Snakemake configfile — see design doc §8)
# ---------------------------------------------------------------------------
with open("config/datasets.yaml") as f:
    DATASETS = yaml.safe_load(f)

wildcard_constraints:
    dataset="|".join(DATASETS["datasets"].keys()),

# ---------------------------------------------------------------------------
# Includes
# ---------------------------------------------------------------------------
# Rule filenames carry a two-level prefix: <chapter><section>, e.g. 1a_pcr_qc.
# The chapter is the analysis block, the section the position inside it — so
# `ls workflow/rules/` reads in analysis order. Chapters are narrative order,
# NOT dependency depth: same-chapter files have no order between them.
#
#   1  data plausibility + reference layer          (1a pcr_qc, 1b features, 1c annotate)
#   2  characterisation that does not need clusters (2a coverage, 2b verification)
#   3  group-level analysis in fitness space        (3a coherence)
#
# The clustering block and everything downstream of it (enrichment, ml, ...) are
# still unprefixed in both places; they take the next free chapter when they do.
#
# Every prefixed stage carries the same prefix into results/, so `ls results/`
# reads in the same analysis order as `ls workflow/rules/`: results/1a_pcr_qc/,
# results/1b_features/, results/1c_annotation/ (the stage noun — the rule file is
# the verb), results/2a_coverage/, results/2b_verification/, results/3a_coherence/.
# The prefix stops at the rule file and the results/logs directories: the script
# and library directories keep the unprefixed stage name (workflow/scripts/
# coherence/, workflow/src/coherence/).
#
# Included in prefix order. clustering.smk must stay ahead of enrichment /
# ml / comparison — they call its selected_variant() / final_clusters_path();
# 1c_annotate.smk must stay ahead of 2a_coverage.smk — coverage reads _ANNOT_REF.
include: "workflow/rules/1a_pcr_qc.smk"
include: "workflow/rules/1b_features.smk"
include: "workflow/rules/1c_annotate.smk"
include: "workflow/rules/2a_coverage.smk"
include: "workflow/rules/2b_verification.smk"
include: "workflow/rules/3a_coherence.smk"
include: "workflow/rules/clustering.smk"
include: "workflow/rules/enrichment.smk"
include: "workflow/rules/enrichment_network.smk"
include: "workflow/rules/ml.smk"
include: "workflow/rules/noncoding_rna.smk"
include: "workflow/rules/comparison.smk"
include: "workflow/rules/utr.smk"
include: "workflow/rules/domain_differences.smk"

# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
# Following the repo convention, per-stage targets are listed but commented —
# uncomment (or pass on the CLI) to run a specific stage. The core chain is:
#   clustering spine -> finalize VARIANT -> enrichment / ml.
# Finalize has named variants (config.clustering.variants), each clustering for
# itself (no fixed candidate stage). Every variant produces
# results/clustering/{dataset}/{variant}/final_clusters.tsv (+ metrics.tsv):
# direct/auto_merge/grid via deterministic scripts, manual_merge via
# finalize_manual_merge (executes notebooks/clustering/finalize_gene_clusters.ipynb
# headlessly). enrichment fans out per variant; ml uses config.clustering.selected_variant.
_REF = DATASETS["reference"]["pombase_version"]
_DATASET = DATASETS["default_dataset"]
_SELECTED_VARIANT = config["clustering"]["selected_variant"]
# Coherence fans out by grouping database; its figures are one set PER source, so
# rule all has to expand over this rather than name one folder.
_COHERENCE_SOURCES = config["coherence"].get("sources", ["go_macrocomplex", "go_cc", "go_bp"])

rule all:
    input:
        # f"results/1b_features/{_REF}/pombe_coding_gene_protein_features.tsv",
        # Selected finalize variant's clusters (buildable variants only):
        # f"results/clustering/{_DATASET}/{_SELECTED_VARIANT}/final_clusters.tsv",
        # Compare ALL buildable variants (builds every variant + a metrics table):
        # f"results/clustering/{_DATASET}/variant_metrics_comparison.tsv",
        # f"results/clustering/{_DATASET}/all_variants_cluster_scatter.pdf",
        # Enrichment (per variant):
        # f"results/enrichment/raw/{_DATASET}/{_SELECTED_VARIANT}/{_REF}/go_enrichment_full_filtered.tsv",
        # Network enrichment (optional, hits STRING/REVIGO — run explicitly):
        # f"results/enrichment/network/{_DATASET}/{_SELECTED_VARIANT}/{_REF}/go_enrichment_full_revigo.tsv",
        # ML AutoML (target x mode; uses selected_variant):
        # f"results/ml/models/{_DATASET}/{_REF}/DR_Explain/metrics.tsv",
        # f"results/ml/models/{_DATASET}/{_REF}/DL_Explain/metrics.tsv",
        # Library-prep QC (no dataset wildcard; spike-in stats feed the figure's
        # panel d, so both targets come from 1a_pcr_qc.smk):
        # "results/1a_pcr_qc/PCR_quality_control.pdf",
        # "results/1a_pcr_qc/spike_in_stats.tsv",
        f"results/2a_coverage/{_DATASET}/coverage_stats.tsv",
        f"results/2a_coverage/{_DATASET}/detailed_genes.xlsx",
        # Ten single-figure PDFs (cnsplots) — one coverage composition + a DR and
        # a DL distribution per annotation dimension, + insertion placement; see
        # 2a_coverage.smk for the full list.
        f"results/2a_coverage/{_DATASET}/coverage_by_deletion_essentiality.pdf",
        # f"results/2a_coverage/{_DATASET}/coverage_by_characterisation_status.pdf",
        # f"results/2a_coverage/{_DATASET}/coverage_by_FYPOviability.pdf",
        # f"results/2a_coverage/{_DATASET}/coverage_insertion_placement.pdf",
        # f"results/2a_coverage/{_DATASET}/coverage_dr_by_deletion_essentiality.pdf",
        # f"results/2a_coverage/{_DATASET}/coverage_dl_by_deletion_essentiality.pdf",
        # f"results/2a_coverage/{_DATASET}/coverage_dr_by_characterisation_status.pdf",
        # f"results/2b_verification/{_DATASET}/verification_stats.tsv",
        # f"results/2b_verification/{_DATASET}/deletion_library_comparison.pdf",
        # f"results/2b_verification/{_DATASET}/verification_category_boxplot.pdf",
        # f"results/2b_verification/{_DATASET}/critical_genes",
        # f"results/2b_verification/{_DATASET}/depletion_curves",
        # f"results/noncoding_rna/{_DATASET}/ncrna_stats.tsv",
        # Batch B (requires resources/curated/final_clusters.tsv):
        # f"results/comparison/{_DATASET}/fitness_correlation_stats.tsv",
        # 3a_coherence is organised as one folder per view — {source}/, combined/
        # and dedup/ — with the FOLDER name (never the file name) saying which view
        # a file belongs to, so each view exposes the same names where the artifact
        # exists for it. The view-neutral inputs live at the dataset root.
        # dedup/ — the de-duplicated representative set as a first-class result: its
        # own terms table, coherence figure, attribution run and cohort workbook.
        f"results/3a_coherence/{_DATASET}/dedup/coherence_terms_representatives.tsv",
        f"results/3a_coherence/{_DATASET}/dedup/coherence_group_members_long.tsv",
        f"results/3a_coherence/{_DATASET}/dedup/coherence.pdf",
        f"results/3a_coherence/{_DATASET}/dedup/incoherence_attribution.tsv",
        f"results/3a_coherence/{_DATASET}/dedup/incoherence_attribution.pdf",
        f"results/3a_coherence/{_DATASET}/dedup/coherence_cohorts.xlsx",
        # The de-duplication audit pages (Altair + pyvis, ~1-3 MB), and the
        # interactive explorer over both pooled views. The per-source explorer is
        # deliberately NOT here: it is the same page over one source's rows, and
        # pre-building five of them for the one you are reading is not worth it
        # (build it by naming the path, e.g. .../go_cc/interactive_scatter.html).
        f"results/3a_coherence/{_DATASET}/dedup/redundancy_network.html",
        f"results/3a_coherence/{_DATASET}/dedup/redundancy_overview.html",
        f"results/3a_coherence/{_DATASET}/combined/interactive_scatter.html",
        f"results/3a_coherence/{_DATASET}/dedup/interactive_scatter.html",
        # combined/ — one figure with every source on it, coloured by source, + the
        # representative set as a sixth series, and the threshold-filtered cohort
        # workbook over the pooled metrics.
        f"results/3a_coherence/{_DATASET}/combined/coherence.pdf",
        f"results/3a_coherence/{_DATASET}/combined/coherence_cohorts.xlsx",
        # Per-source coherence figures — one coherence overview + one named-group
        # scatter grid per registered source (3 sources = 6 figures + their
        # .review.png siblings), each with its own cohort workbook sliced out of the
        # combined table so all views share one q scale.
        expand(f"results/3a_coherence/{_DATASET}/{{source}}/coherence.pdf",
               source=_COHERENCE_SOURCES),
        expand(f"results/3a_coherence/{_DATASET}/{{source}}/group_scatter.pdf",
               source=_COHERENCE_SOURCES),
        expand(f"results/3a_coherence/{_DATASET}/{{source}}/coherence_cohorts.xlsx",
               source=_COHERENCE_SOURCES),
        # Incoherence attribution (why a group is dispersed), per configured source:
        expand(f"results/3a_coherence/{_DATASET}/{{source}}/incoherence_attribution.tsv",
               source=config["coherence"].get("attribution_sources", ["go_macrocomplex"])),
        expand(f"results/3a_coherence/{_DATASET}/{{source}}/incoherence_attribution.pdf",
               source=config["coherence"].get("attribution_sources", ["go_macrocomplex"])),
        # Batch C (requires insertion-level results):
        # f"results/utr/{_DATASET}/utr_insertion_stats.tsv",
        # f"results/domain_differences/{_DATASET}/domain_candidate_stats.tsv",
    message:
        "*** DIT-HAP analysis complete"
