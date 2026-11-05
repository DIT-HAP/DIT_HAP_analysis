# =============================================================================
# 3c_enrichment.smk — GO/FYPO/MONDO cluster enrichment + STRING/REVIGO network
# =============================================================================
#
# Over-representation analysis for each merged cluster in a finalize variant's
# final_clusters.tsv, then the two optional network passes over its results.
#
# Deterministic half (locally computed, goatools + PomBase OBO/GAF):
#   prepare  -> gene sets (background/nonWT/per-cluster) + id->name + gene lists
#   enrich   -> one job per ontology (GO / FYPO / MONDO), each writing workbooks
#   finalize -> concat TSVs + the filtered GO table the figures read
#
# Non-deterministic half (external APIs, NOT in `rule all`):
#   STRING-db functional enrichment and REVIGO semantic-similarity annotation.
#   Both are cached under resources/external/enrichment_cache/{dataset}, so a
#   second run is deterministic and offline.
#
# Fans out over clustering finalize VARIANTS: every variant's clusters get their
# own enrichment under {dataset}/{variant}/{pombase_version}/. final_clusters.tsv
# is sourced via final_clusters_path(dataset, variant) (3b_clustering.smk) and is
# BUILT by the matching finalize rule: direct/auto_merge/grid via deterministic
# scripts, manual_merge via finalize_manual_merge.
#
# Figure rules live here too and read only the finalized TSVs, never the core:
#   plot_enrichment                      one dotplot per variant (opt-in)
#   plot_variant_enrichment_comparison   heatmap + counts across every variant

_ENRICH = config["enrichment"]
_ENRICH_VARIANTS = list(config["clustering"]["variants"])
_ENRICH_ONTOLOGIES = ["GO", "FYPO", "MONDO"]

_ERAW = "results/3c_enrichment/{dataset}/{variant}/{pombase_version}"
_EWORK = f"{_ERAW}/_work"
_ENETWORK = f"{_ERAW}/network"
# Same path with the non-expanded wildcards escaped, for expand() calls that
# enumerate variants and leave dataset/pombase_version for Snakemake to fill.
_ERAW_PARTIAL = "results/3c_enrichment/{{dataset}}/{{variant}}/{{pombase_version}}"

wildcard_constraints:
    ontology="|".join(_ENRICH_ONTOLOGIES),


# --- Preprocessing spine (gene sets + id->name + gene lists) ---
rule prepare_genesets:
    input:
        final_clusters=lambda wc: final_clusters_path(wc.dataset, wc.variant),
        pombase_dir="resources/external/pombase/{pombase_version}",
        deletion_library_xlsx="resources/curated/deletion_library_categories.xlsx",
    output:
        genesets=f"{_EWORK}/genesets.pkl",
        id2name=f"{_EWORK}/id2name.pkl",
        all_genes=f"{_ERAW}/DIT_HAP_all_genes.txt",
    params:
        output_dir=_ERAW,
        work_dir=_EWORK,
        wt_cluster=_ENRICH["wt_cluster"],
    log:
        "logs/3c_enrichment/prepare_genesets_{dataset}_{variant}_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [enrichment] Preparing gene sets for {wildcards.dataset}/{wildcards.variant}..."
    shell:
        """
        python workflow/scripts/enrichment/prepare_genesets.py \
            --final-clusters {input.final_clusters} \
            --pombase-dir {input.pombase_dir} \
            --deletion-library-xlsx {input.deletion_library_xlsx} \
            --output-dir {params.output_dir} \
            --work-dir {params.work_dir} \
            --wt-cluster {params.wt_cluster} &> {log}
        """


# --- One ontology (fanned out by the `ontology` wildcard) ---
rule enrich_one_ontology:
    input:
        genesets=f"{_EWORK}/genesets.pkl",
        id2name=f"{_EWORK}/id2name.pkl",
        pombase_dir="resources/external/pombase/{pombase_version}",
    output:
        frames=f"{_EWORK}/{{ontology}}_frames.pkl",
    params:
        output_dir=_ERAW,
        work_dir=_EWORK,
        wt_cluster=_ENRICH["wt_cluster"],
        fdr_threshold=_ENRICH["fdr_threshold"],
    log:
        "logs/3c_enrichment/enrich_{ontology}_{dataset}_{variant}_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [enrichment] Running {wildcards.ontology} enrichment for {wildcards.dataset}/{wildcards.variant}..."
    shell:
        """
        python workflow/scripts/enrichment/enrich_one_ontology.py \
            --ontology {wildcards.ontology} \
            --genesets {input.genesets} \
            --id2name {input.id2name} \
            --pombase-dir {input.pombase_dir} \
            --output-dir {params.output_dir} \
            --work-dir {params.work_dir} \
            --wt-cluster {params.wt_cluster} \
            --fdr-threshold {params.fdr_threshold} &> {log}
        """


# --- Finalize (concat TSVs + filtered GO target) ---
rule finalize_enrichment:
    input:
        frames=expand(f"{_EWORK}/{{ontology}}_frames.pkl", ontology=_ENRICH_ONTOLOGIES, allow_missing=True),
    output:
        go_filtered=f"{_ERAW}/go_enrichment_full_filtered.tsv",
        go_full=f"{_ERAW}/go_enrichment_full.tsv",
    params:
        work_dir=_EWORK,
        output_dir=_ERAW,
        pop_count_max=_ENRICH["pop_count_max"],
    log:
        "logs/3c_enrichment/finalize_enrichment_{dataset}_{variant}_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [enrichment] Finalizing enrichment TSVs for {wildcards.dataset}/{wildcards.variant}..."
    shell:
        """
        python workflow/scripts/enrichment/finalize_enrichment.py \
            --work-dir {params.work_dir} \
            --output-dir {params.output_dir} \
            --pop-count-max {params.pop_count_max} &> {log}
        """


# --- Per-variant GO dotplot (one page per variant; opt-in) ---
rule plot_enrichment:
    input:
        enrichment=f"{_ERAW}/go_enrichment_full_filtered.tsv",
    output:
        dotplot=f"{_ERAW}/go_enrichment_dotplot.pdf",
        preview=f"{_ERAW}/go_enrichment_dotplot.review.png",
    params:
        top_terms_per_cluster=_ENRICH["dotplot_top_terms"],
    log:
        "logs/3c_enrichment/plot_enrichment_{dataset}_{variant}_{pombase_version}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [enrichment] Plotting the GO dotplot for {wildcards.dataset}/{wildcards.variant}..."
    shell:
        """
        python workflow/scripts/enrichment/plot_enrichment.py \
            --enrichment {input.enrichment} \
            --output {output.dotplot} \
            --variant-label {wildcards.variant} \
            --top-terms-per-cluster {params.top_terms_per_cluster} &> {log}
        """


# --- Cross-variant comparison (heatmap + counts; opt-in) ---
# Depends on EVERY configured variant, so building it also builds the enrichment
# for variants that were never run on their own.
rule plot_variant_enrichment_comparison:
    input:
        filtered=expand(
            f"{_ERAW_PARTIAL}/go_enrichment_full_filtered.tsv",
            variant=_ENRICH_VARIANTS,
        ),
    output:
        heatmap="results/3c_enrichment/{dataset}/variant_enrichment_heatmap.pdf",
        heatmap_preview="results/3c_enrichment/{dataset}/variant_enrichment_heatmap.review.png",
        counts="results/3c_enrichment/{dataset}/variant_enrichment_counts.pdf",
        counts_preview="results/3c_enrichment/{dataset}/variant_enrichment_counts.review.png",
    params:
        output_dir="results/3c_enrichment/{dataset}",
        labels=" ".join(_ENRICH_VARIANTS),
        top_terms=_ENRICH["variant_heatmap_top_terms"],
        fdr_threshold=_ENRICH["fdr_threshold"],
    log:
        "logs/3c_enrichment/plot_variant_enrichment_comparison_{dataset}.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [enrichment] Plotting the cross-variant enrichment comparison for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/enrichment/plot_variant_enrichment_comparison.py \
            --input {input.filtered} \
            --label {params.labels} \
            --output-dir {params.output_dir} \
            --top-terms {params.top_terms} \
            --fdr-threshold {params.fdr_threshold} &> {log}
        """


# --- Network half: STRING-db (NON-deterministic, opt-in) ---
rule run_string_enrichment:
    input:
        all_genes=f"{_ERAW}/DIT_HAP_all_genes.txt",
    output:
        string_workbook=f"{_ENETWORK}/string_enrichment_results.xlsx",
    params:
        enrichment_dir=_ERAW,
        output_dir=_ENETWORK,
        cache_dir="resources/external/enrichment_cache/{dataset}",
    log:
        "logs/3c_enrichment/run_string_enrichment_{dataset}_{variant}_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** Running STRING network enrichment for {wildcards.dataset} (needs network on first run)..."
    shell:
        """
        python workflow/scripts/enrichment/run_string_enrichment.py \
            --enrichment-dir {params.enrichment_dir} \
            --output-dir {params.output_dir} \
            --cache-dir {params.cache_dir} &> {log}
        """


# --- Network half: REVIGO (NON-deterministic, opt-in) ---
rule run_revigo_annotation:
    input:
        go_full=f"{_ERAW}/go_enrichment_full.tsv",
    output:
        go_revigo=f"{_ENETWORK}/go_enrichment_full_revigo.tsv",
    params:
        enrichment_dir=_ERAW,
        output_dir=_ENETWORK,
        cache_dir="resources/external/enrichment_cache/{dataset}",
        revigo_cutoffs=" ".join(str(cutoff) for cutoff in _ENRICH["network"]["revigo_cutoffs"]),
    log:
        "logs/3c_enrichment/run_revigo_annotation_{dataset}_{variant}_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** Running REVIGO GO annotation for {wildcards.dataset} (needs network on first run)..."
    shell:
        """
        python workflow/scripts/enrichment/run_revigo_annotation.py \
            --enrichment-dir {params.enrichment_dir} \
            --output-dir {params.output_dir} \
            --cache-dir {params.cache_dir} \
            --revigo-cutoffs {params.revigo_cutoffs} &> {log}
        """
