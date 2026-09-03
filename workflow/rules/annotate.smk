# =============================================================================
# annotate.smk — Per-gene annotation reference + annotated tables
# =============================================================================
#
# Parses PomBase + SGD sources once into a single wide annotation table (24
# columns) keyed by pombe systematic ID, then joins it onto tables that need it.
# The reference is dataset-independent (like features.smk): it depends only on a
# PomBase version and an SGD snapshot, not on any DIT-HAP sequencing project.
#
# The split between "build the reference" and "annotate a table" exists because
# parsing the GO OBO/GAF and the 200k-row SGD phenotype table costs far more than
# the join; build_annotation_reference pays that cost once and caches it as a
# parquet (per repo convention: intermediates are parquet, not TSV).
#
# Four rules, in dependency order:
#   fetch_sgd_data             -> resources/external/sgd/{sgd_version}/  (NETWORK)
#   build_annotation_reference -> gene_annotation_reference.parquet
#   annotate_table             -> {name}.annotated.tsv   (per config.annotate.tables)
#   build_annotated_workbook   -> {dataset}_annotated.xlsx (coverage + verification)
#
# The SGD snapshot is pinned by config.annotate.sgd_version rather than
# auto-detected from whatever directories exist, so re-fetching SGD cannot
# silently change an existing annotation; _SGD_VERSION below is the single place
# the other rules read it from.
_ANNOT_CFG = config.get("annotate", {})
_SGD_VERSION = _ANNOT_CFG.get("sgd_version", "2026-08-11")
_ANNOT_TABLES = _ANNOT_CFG.get("tables", {})
_ANNOT_REF = (
    f"results/annotation/{DATASETS['reference']['pombase_version']}/"
    f"{_SGD_VERSION}/gene_annotation_reference.parquet"
)

# `name` is a config key, not a path component of the input — constrain it to the
# registered names so a typo fails at DAG build instead of matching greedily.
wildcard_constraints:
    name="|".join(_ANNOT_TABLES.keys()) if _ANNOT_TABLES else "$^",


# -----------------------------------------------------------------------------
# Stage 0: Fetch the SGD snapshot (network)
# -----------------------------------------------------------------------------
# The ONLY rule in this repo that makes a network call, and it is not reachable
# from `rule all` — resources/external/{pombase,biogrid,sgd}/ are hand-filled,
# git-ignored local caches. It is a rule rather than a bare script so
# build_annotation_reference can declare the SGD directory as an input and have
# it materialise on demand.
#
# Marked `retries: 2` (phenotype_data.tab is ~39 MB and SGD can be slow) on top
# of curl's own --retry. Downloads land in .part files and are renamed, so an
# interrupted run cannot leave a truncated table looking complete.
#
# The directory is NOT `temp()`: it is a versioned snapshot the reference's
# provenance points at, and re-downloading 39 MB to rebuild a parquet would be
# wasteful. Delete it by hand to force a re-fetch.
rule fetch_sgd_data:
    output:
        sgd_dir=directory("resources/external/sgd/{sgd_version}"),
    retries: 2
    log:
        "logs/annotate/fetch_sgd_data_{sgd_version}.log",
    message:
        "*** [annotate] Fetching SGD tables for snapshot {wildcards.sgd_version} (network)..."
    shell:
        """
        bash workflow/scripts/annotate/fetch_sgd_data.sh {wildcards.sgd_version} &> {log}
        """


# -----------------------------------------------------------------------------
# Stage 1: Build the annotation reference
# -----------------------------------------------------------------------------
rule build_annotation_reference:
    input:
        pombase_dir="resources/external/pombase/{pombase_version}",
        sgd_dir="resources/external/sgd/{sgd_version}",
        deletion_library_xlsx="resources/curated/deletion_library_categories.xlsx",
        grna_parameters_tsv="resources/curated/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv",
    output:
        reference="results/annotation/{pombase_version}/{sgd_version}/gene_annotation_reference.parquet",
    log:
        "logs/annotate/build_annotation_reference_{pombase_version}_{sgd_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [annotate] Building gene annotation reference (PomBase {wildcards.pombase_version}, SGD {wildcards.sgd_version})..."
    shell:
        """
        python workflow/scripts/annotate/build_annotation_reference.py \
            --pombase-dir {input.pombase_dir} \
            --sgd-dir {input.sgd_dir} \
            --deletion-library-xlsx {input.deletion_library_xlsx} \
            --grna-parameters-tsv {input.grna_parameters_tsv} \
            --output {output.reference} &> {log}
        """


# -----------------------------------------------------------------------------
# Stage 2: Annotate a registered table
# -----------------------------------------------------------------------------
# annotate_pombe_genes.py joins the reference onto ANY table with a column of
# pombe systematic IDs, so it has no single natural DAG position. Rather than
# leaving it as a hand-typed command, the tables worth annotating are registered
# in config.annotate.tables (name -> input path + gene column), and {name}
# selects one. A one-off join is then a config entry, which is reproducible and
# reviewable; for a genuinely throwaway table, calling the script directly still
# works exactly as before.
#
# Unmatched IDs are reported in the log, never dropped (unless the entry sets
# drop_unmatched) — a wrong gene_column looks exactly like a table with no
# annotations, so that summary is the only signal.
rule annotate_table:
    input:
        table=lambda wc: _ANNOT_TABLES[wc.name]["input"].format(
            dataset=wc.dataset, variant=wc.variant
        ),
        reference=_ANNOT_REF,
    output:
        annotated="results/annotation/{dataset}/{variant}/tables/{name}.annotated.tsv",
    params:
        gene_column=lambda wc: _ANNOT_TABLES[wc.name]["gene_column"],
        columns_flag=lambda wc: (
            "--columns " + " ".join(_ANNOT_TABLES[wc.name]["columns"])
            if _ANNOT_TABLES[wc.name].get("columns") else ""
        ),
        drop_flag=lambda wc: (
            "--drop-unmatched" if _ANNOT_TABLES[wc.name].get("drop_unmatched", False) else ""
        ),
    log:
        "logs/annotate/annotate_table_{dataset}_{variant}_{name}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [annotate] Annotating table '{wildcards.name}' for {wildcards.dataset} x {wildcards.variant}..."
    shell:
        """
        python workflow/scripts/annotate/annotate_pombe_genes.py \
            --input {input.table} \
            --gene-column '{params.gene_column}' \
            --annotation-reference {input.reference} \
            {params.columns_flag} \
            {params.drop_flag} \
            --output {output.annotated} &> {log}
        """


# -----------------------------------------------------------------------------
# Stage 3: Consolidated annotated workbook
# -----------------------------------------------------------------------------
# Consolidates coverage detailed_genes.xlsx (16 sheets: All genes + 15 categories)
# and verification critical_genes/*.tsv (4 groups) into one 20-sheet annotated
# workbook, each sheet name carrying its gene count. Both inputs are fixed
# Snakemake paths, so unlike annotate_table this needs no registry — naming the
# output builds coverage and verification first.
#
# critical_genes is a directory() output of verification_boxplots; the script
# globs *.tsv inside it, so the directory is the input and the glob stays in the
# script (Snakemake cannot enumerate files it did not declare).
rule build_annotated_workbook:
    input:
        detailed_xlsx="results/coverage/{dataset}/detailed_genes.xlsx",
        critical_dir="results/verification/{dataset}/critical_genes",
        reference=_ANNOT_REF,
    output:
        workbook="results/annotation/{dataset}/{dataset}_annotated.xlsx",
    log:
        "logs/annotate/build_annotated_workbook_{dataset}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [annotate] Building consolidated annotated workbook for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/annotate/build_annotated_workbook.py \
            --detailed-xlsx {input.detailed_xlsx} \
            --critical-dir {input.critical_dir} \
            --annotation-reference {input.reference} \
            --output {output.workbook} &> {log}
        """
