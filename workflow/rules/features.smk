# =============================================================================
# features.smk — Pombe coding gene feature collection (dataset-independent)
# =============================================================================
#
# Split by biological level: the DNA level is the "spine" (it builds the
# gffutils DB and enumerates the coding-gene set), and the other five levels
# read the DNA pickle to recover that set before assembling independently. A
# final merge rule joins the six per-level pickles into the feature matrix.
#
# Depends only on the reference PomBase version, not on any DIT-HAP sequencing
# project — no `dataset` wildcard (design doc §8). Per-level intermediates are
# pickles (not TSV) so bool/int dtypes and the intentional duplicate
# DeletionLibrary_essentiality column survive round-trip; only the final matrix
# is a TSV.

_LEVELS = "results/features/{pombase_version}/_levels"


# --- DNA level (the spine: builds the DB, enumerates coding genes, codon usage) ---
rule collect_dna_features:
    input:
        pombase_dir="resources/external/pombase/{pombase_version}",
        genome_landmarks="config/genome_landmarks.yaml",
    output:
        dna=f"{_LEVELS}/dna_features.parquet",
        codon_usage="results/features/{pombase_version}/codon_usage_matrix.tsv",
    log:
        "logs/features/collect_dna_features_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [DNA] Collecting DNA-level features for PomBase {wildcards.pombase_version}..."
    shell:
        """
        python workflow/scripts/features/collect_dna_features.py \
            --pombase-dir {input.pombase_dir} \
            --genome-landmarks {input.genome_landmarks} \
            --output {output.dna} \
            --codon-usage-output {output.codon_usage} &> {log}
        """


# --- RNA level ---
rule collect_rna_features:
    input:
        pombase_dir="resources/external/pombase/{pombase_version}",
        literature_dir="resources/literature",
        dna=f"{_LEVELS}/dna_features.parquet",
    output:
        rna=f"{_LEVELS}/rna_features.parquet",
    log:
        "logs/features/collect_rna_features_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [RNA] Collecting RNA-level features for PomBase {wildcards.pombase_version}..."
    shell:
        """
        python workflow/scripts/features/collect_rna_features.py \
            --pombase-dir {input.pombase_dir} \
            --literature-dir {input.literature_dir} \
            --dna-features {input.dna} \
            --output {output.rna} &> {log}
        """


# --- Protein level (slowest: walks the AlphaFold structure directory) ---
rule collect_protein_features:
    input:
        pombase_dir="resources/external/pombase/{pombase_version}",
        alphafold_dir=DATASETS["reference"]["alphafold_dir"],
        literature_dir="resources/literature",
        dna=f"{_LEVELS}/dna_features.parquet",
    output:
        protein=f"{_LEVELS}/protein_features.parquet",
    log:
        "logs/features/collect_protein_features_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [protein] Collecting protein-level features for PomBase {wildcards.pombase_version}..."
    shell:
        """
        python workflow/scripts/features/collect_protein_features.py \
            --pombase-dir {input.pombase_dir} \
            --alphafold-dir {input.alphafold_dir} \
            --literature-dir {input.literature_dir} \
            --dna-features {input.dna} \
            --output {output.protein} &> {log}
        """


# --- Evolutionary level ---
rule collect_evolutionary_features:
    input:
        pombase_dir="resources/external/pombase/{pombase_version}",
        literature_dir="resources/literature",
        ensembl_paralogs_tsv="resources/external/ensembl/pombe_paralog_from_ensemble_biomart_export.tsv",
        dna=f"{_LEVELS}/dna_features.parquet",
    output:
        evolutionary=f"{_LEVELS}/evolutionary_features.parquet",
    log:
        "logs/features/collect_evolutionary_features_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [evolutionary] Collecting evolutionary-level features for PomBase {wildcards.pombase_version}..."
    shell:
        """
        python workflow/scripts/features/collect_evolutionary_features.py \
            --pombase-dir {input.pombase_dir} \
            --literature-dir {input.literature_dir} \
            --ensembl-paralogs-tsv {input.ensembl_paralogs_tsv} \
            --dna-features {input.dna} \
            --output {output.evolutionary} &> {log}
        """


# --- Network level (loads GO DAG/GAF for term richness) ---
rule collect_network_features:
    input:
        pombase_dir="resources/external/pombase/{pombase_version}",
        biogrid_tsv="resources/external/biogrid/BIOGRID-ORGANISM-Schizosaccharomyces_pombe_972h-5.0.251.tab3.txt",
        dna=f"{_LEVELS}/dna_features.parquet",
    output:
        network=f"{_LEVELS}/network_features.parquet",
    log:
        "logs/features/collect_network_features_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [network] Collecting network-level features for PomBase {wildcards.pombase_version}..."
    shell:
        """
        python workflow/scripts/features/collect_network_features.py \
            --pombase-dir {input.pombase_dir} \
            --biogrid-tsv {input.biogrid_tsv} \
            --dna-features {input.dna} \
            --output {output.network} &> {log}
        """


# --- Phenotype level ---
rule collect_phenotype_features:
    input:
        pombase_dir="resources/external/pombase/{pombase_version}",
        literature_dir="resources/literature",
        deletion_library_xlsx="resources/curated/deletion_library_categories.xlsx",
        essentiality_verification_csv="resources/curated/essentiality_verification.csv",
        dna=f"{_LEVELS}/dna_features.parquet",
    output:
        phenotype=f"{_LEVELS}/phenotype_features.parquet",
    log:
        "logs/features/collect_phenotype_features_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [phenotype] Collecting phenotype-level features for PomBase {wildcards.pombase_version}..."
    shell:
        """
        python workflow/scripts/features/collect_phenotype_features.py \
            --pombase-dir {input.pombase_dir} \
            --literature-dir {input.literature_dir} \
            --deletion-library-xlsx {input.deletion_library_xlsx} \
            --essentiality-verification-csv {input.essentiality_verification_csv} \
            --dna-features {input.dna} \
            --output {output.phenotype} &> {log}
        """


# --- Merge (outer-join the six per-level pickles into the final matrix) ---
rule merge_pombe_features:
    input:
        pombase_dir="resources/external/pombase/{pombase_version}",
        dna=f"{_LEVELS}/dna_features.parquet",
        rna=f"{_LEVELS}/rna_features.parquet",
        protein=f"{_LEVELS}/protein_features.parquet",
        evolutionary=f"{_LEVELS}/evolutionary_features.parquet",
        network=f"{_LEVELS}/network_features.parquet",
        phenotype=f"{_LEVELS}/phenotype_features.parquet",
    output:
        features="results/features/{pombase_version}/pombe_coding_gene_protein_features.tsv",
    log:
        "logs/features/merge_pombe_features_{pombase_version}.log",
    conda:
        "../envs/biopython.yml"
    message:
        "*** [merge] Merging per-level features for PomBase {wildcards.pombase_version}..."
    shell:
        """
        python workflow/scripts/features/merge_features.py \
            --pombase-dir {input.pombase_dir} \
            --dna-features {input.dna} \
            --rna-features {input.rna} \
            --protein-features {input.protein} \
            --evolutionary-features {input.evolutionary} \
            --network-features {input.network} \
            --phenotype-features {input.phenotype} \
            --output {output.features} &> {log}
        """
