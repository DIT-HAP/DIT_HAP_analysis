# =============================================================================
# enrichment_network.smk — STRING + REVIGO network enrichment (optional)
# =============================================================================

# NON-deterministic: hits STRING-db + REVIGO web APIs. Not in `rule all` — invoke
# explicitly (design doc §5). Responses are cached under resources/external/
# enrichment_cache/{dataset}, so a second run is deterministic and offline.
#
# Consumes the deterministic enrichment outputs (gene lists + go_enrichment_full.tsv)
# from run_ontology_enrichment (enrichment.smk).
rule run_network_enrichment:
    input:
        go_full="results/enrichment/raw/{dataset}/{pombase_version}/go_enrichment_full.tsv",
        all_genes="results/enrichment/raw/{dataset}/{pombase_version}/DIT_HAP_all_genes.txt",
    output:
        string_workbook="results/enrichment/network/{dataset}/{pombase_version}/string_enrichment_results.xlsx",
        go_revigo="results/enrichment/network/{dataset}/{pombase_version}/go_enrichment_full_revigo.tsv",
    params:
        enrichment_dir="results/enrichment/raw/{dataset}/{pombase_version}",
        output_dir="results/enrichment/network/{dataset}/{pombase_version}",
        cache_dir="resources/external/enrichment_cache/{dataset}",
        wt_cluster=config.get("enrichment", {}).get("wt_cluster", 9),
        revigo_cutoffs=" ".join(
            str(c) for c in config.get("enrichment_network", {}).get("revigo_cutoffs", [0.7, 0.5])
        ),
    log:
        "logs/enrichment/run_network_enrichment_{dataset}_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** Running STRING/REVIGO network enrichment for {wildcards.dataset} (needs network on first run)..."
    shell:
        """
        python workflow/scripts/enrichment/run_network_enrichment.py \
            --enrichment-dir {params.enrichment_dir} \
            --output-dir {params.output_dir} \
            --cache-dir {params.cache_dir} \
            --wt-cluster {params.wt_cluster} \
            --revigo-cutoffs {params.revigo_cutoffs} &> {log}
        """
