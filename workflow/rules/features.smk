# =============================================================================
# features.smk — Pombe coding gene feature collection (dataset-independent)
# =============================================================================

# Depends only on the reference PomBase version, not on any DIT-HAP sequencing
# project — no `dataset` wildcard (design doc §8).
rule collect_pombe_features:
    input:
        pombase_dir="resources/external/pombase/{pombase_version}",
        alphafold_dir=DATASETS["reference"]["alphafold_dir"],
        literature_dir="resources/literature",
        deletion_library_xlsx="resources/curated/deletion_library_categories.xlsx",
        essentiality_verification_csv="resources/curated/essentiality_verification.csv",
        biogrid_tsv="resources/external/biogrid/BIOGRID-ORGANISM-Schizosaccharomyces_pombe_972h-5.0.251.tab3.txt",
        ensembl_paralogs_tsv="resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv",
    output:
        features="results/features/{pombase_version}/pombe_coding_gene_protein_features.tsv",
        codon_usage="results/features/{pombase_version}/codon_usage_matrix.tsv",
    log:
        "logs/features/collect_pombe_features_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** Collecting pombe coding gene features for PomBase {wildcards.pombase_version}..."
    shell:
        """
        python workflow/scripts/features/collect_pombe_features.py \
            --pombase-dir {input.pombase_dir} \
            --alphafold-dir {input.alphafold_dir} \
            --literature-dir {input.literature_dir} \
            --deletion-library-xlsx {input.deletion_library_xlsx} \
            --essentiality-verification-csv {input.essentiality_verification_csv} \
            --biogrid-tsv {input.biogrid_tsv} \
            --ensembl-paralogs-tsv {input.ensembl_paralogs_tsv} \
            --output {output.features} \
            --codon-usage-output {output.codon_usage} &> {log}
        """
