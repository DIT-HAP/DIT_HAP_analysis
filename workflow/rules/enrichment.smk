# =============================================================================
# enrichment.smk — Deterministic GO/FYPO/MONDO cluster enrichment
# =============================================================================

# Per-dataset: over-representation analysis for each merged cluster in the
# curated final_clusters.tsv. Pure-local (goatools + PomBase OBO/GAF); the
# STRING/REVIGO network steps are in enrichment_network.smk (design doc §5).
#
# final_clusters.tsv is an UN-BUILDABLE human-curated input (design doc §8): if
# it is missing, Snakemake reports "missing input" — run
# notebooks/clustering/finalize_gene_clusters.ipynb first.
rule run_ontology_enrichment:
    input:
        final_clusters="resources/curated/final_clusters.tsv",
        pombase_dir="resources/external/pombase/{pombase_version}",
        deletion_library_xlsx="resources/curated/deletion_library_categories.xlsx",
    output:
        go_filtered="results/enrichment/raw/{dataset}/{pombase_version}/go_enrichment_full_filtered.tsv",
        go_workbook="results/enrichment/raw/{dataset}/{pombase_version}/gene_ontology_enrichment_results.xlsx",
        fypo_workbook="results/enrichment/raw/{dataset}/{pombase_version}/fission_yeast_phenotype_ontology_enrichment_results.xlsx",
        mondo_workbook="results/enrichment/raw/{dataset}/{pombase_version}/mondo_disease_ontology_enrichment_results.xlsx",
    params:
        output_dir="results/enrichment/raw/{dataset}/{pombase_version}",
        intermediate_dir="results/enrichment/raw/{dataset}/{pombase_version}/_gaf",
        wt_cluster=config.get("enrichment", {}).get("wt_cluster", 9),
        pop_count_max=config.get("enrichment", {}).get("pop_count_max", 400),
        fdr_threshold=config.get("enrichment", {}).get("fdr_threshold", 0.05),
    log:
        "logs/enrichment/run_ontology_enrichment_{dataset}_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** Running GO/FYPO/MONDO enrichment for {wildcards.dataset}..."
    shell:
        """
        python workflow/scripts/enrichment/run_ontology_enrichment.py \
            --final-clusters {input.final_clusters} \
            --pombase-dir {input.pombase_dir} \
            --deletion-library-xlsx {input.deletion_library_xlsx} \
            --output-dir {params.output_dir} \
            --intermediate-dir {params.intermediate_dir} \
            --wt-cluster {params.wt_cluster} \
            --pop-count-max {params.pop_count_max} \
            --fdr-threshold {params.fdr_threshold} &> {log}
        """
