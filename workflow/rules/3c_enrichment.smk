# =============================================================================
# enrichment.smk — Deterministic GO/FYPO/MONDO cluster enrichment
# =============================================================================
#
# Per-dataset: over-representation analysis for each merged cluster in the
# curated final_clusters.tsv. Pure-local (goatools + PomBase OBO/GAF); the
# STRING/REVIGO network steps are in enrichment_network.smk (design doc §5).
#
# Split by ontology, mirroring clustering.smk's method fan-out:
#   prepare  -> gene sets (background/nonWT/per-cluster) + id->name + gene lists
#   enrich   -> one job per ontology (GO / FYPO / MONDO), each writes its workbooks
#   finalize -> concat TSVs + the filtered GO target the network rule consumes
# Per-ontology frames are pickles under _work/; only the final TSVs/workbooks
# are user-facing.
#
# Fans out over clustering finalize VARIANTS: enrichment output nests under a
# {variant} dir so each variant's clusters get their own enrichment (compare them
# to pick one). final_clusters.tsv is sourced via final_clusters_path(dataset,
# variant) (clustering.smk) and BUILT by the matching finalize rule for every
# variant type: direct/auto_merge/grid via deterministic scripts, manual_merge via
# finalize_manual_merge (executes finalize_gene_clusters.ipynb headlessly).

_ENRICH_ONTOLOGIES = ["GO", "FYPO", "MONDO"]
_ERAW = "results/enrichment/raw/{dataset}/{variant}/{pombase_version}"
_EWORK = "results/enrichment/raw/{dataset}/{variant}/{pombase_version}/_work"

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
        wt_cluster=config.get("enrichment", {}).get("wt_cluster", 9),
    log:
        "logs/enrichment/prepare_genesets_{dataset}_{variant}_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [enrichment] Preparing gene sets for {wildcards.dataset}..."
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
        wt_cluster=config.get("enrichment", {}).get("wt_cluster", 9),
        fdr_threshold=config.get("enrichment", {}).get("fdr_threshold", 0.05),
    log:
        "logs/enrichment/enrich_{ontology}_{dataset}_{variant}_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [enrichment] Running {wildcards.ontology} enrichment for {wildcards.dataset}..."
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
        pop_count_max=config.get("enrichment", {}).get("pop_count_max", 400),
    log:
        "logs/enrichment/finalize_enrichment_{dataset}_{variant}_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [enrichment] Finalizing enrichment TSVs for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/enrichment/finalize_enrichment.py \
            --work-dir {params.work_dir} \
            --output-dir {params.output_dir} \
            --pop-count-max {params.pop_count_max} &> {log}
        """
