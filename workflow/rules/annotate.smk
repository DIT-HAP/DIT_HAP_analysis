# =============================================================================
# annotate.smk — Per-gene annotation reference table
# =============================================================================
#
# Parses PomBase + SGD sources once into a single wide annotation table, keyed by
# pombe systematic ID. Dataset-independent (like features.smk): it depends only
# on a PomBase version and an SGD snapshot, not on any DIT-HAP sequencing project.
#
# Annotating an actual table is NOT a rule — it is a general-purpose CLI you run
# by hand on whatever table you have:
#
#   python workflow/scripts/annotate/annotate_pombe_genes.py \
#       --input my_genes.tsv --gene-column gene_systematic_id \
#       --annotation-reference results/annotation/<pombase>/<sgd>/gene_annotation_reference.parquet \
#       --output my_genes.annotated.tsv
#
# The split exists because parsing the GO OBO/GAF and the 200k-row SGD phenotype
# table costs far more than the join; this rule pays that cost once and caches it
# as a parquet (per repo convention: intermediates are parquet, not TSV).
#
# resources/external/sgd/ is populated MANUALLY by
# workflow/scripts/annotate/fetch_sgd_data.sh — consistent with pombase/ and
# biogrid/, which are also hand-filled, git-ignored local caches. No rule in this
# repo makes network calls.


rule build_annotation_reference:
    input:
        pombase_dir="resources/external/pombase/{pombase_version}",
        sgd_dir="resources/external/sgd/{sgd_version}",
        deletion_library_xlsx="resources/curated/deletion_library_categories.xlsx",
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
            --output {output.reference} &> {log}
        """
