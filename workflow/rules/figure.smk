# =============================================================================
# figure.smk — Centralized plotting rules (computation-plotting decoupling)
# =============================================================================
#
# This file centralizes all figure-generation rules (rules that output PDF/PNG).
# Separation of concerns: computation rules (data generation) remain in their
# domain .smk files (clustering.smk, enrichment.smk, coherence.smk, etc.);
# plotting rules (consume data, output figures) are unified here.
#
# Naming conventions (aligned with upstream DIT_HAP_snakemake, ADR-0001):
#   Rule names:  plot_<feature> (e.g. plot_coherence, plot_variant_clusters)
#   Environment: conda: "../envs/cnsplots.yml"
#   Output paths: direct content description, no type suffix (coherence.pdf, not coherence_figure.pdf)
#
# Usage pattern:
#   rule plot_example:
#       conda: "../envs/cnsplots.yml"
#       input: "results/example/{dataset}/data.parquet"
#       output: "results/example/{dataset}/example.pdf"
#       script: "../scripts/example/plot_example.py"

# -----------------------------------------------------------------------------
# PCR Quality Control
# -----------------------------------------------------------------------------
_PCRWORK = "results/pcr_qc/_work"

rule plot_pcr_qc:
    input:
        pbl_pbr=f"{_PCRWORK}/pbl_pbr.parquet",
        tech=f"{_PCRWORK}/tech.parquet",
        bio=f"{_PCRWORK}/bio.parquet",
        spikein=f"{_PCRWORK}/spikein.parquet",
    output:
        "results/pcr_qc/PCR_quality_control.pdf",
    log:
        "logs/pcr_qc/plot_pcr_qc.log",
    conda:
        "../envs/cnsplots.yml"
    message:
        "*** [pcr_qc] Building 2x2 library-prep QC figure..."
    shell:
        """
        python workflow/scripts/pcr_qc/plot_pcr_qc.py \
            --pbl-pbr {input.pbl_pbr} \
            --tech {input.tech} \
            --bio {input.bio} \
            --spikein {input.spikein} \
            --output {output} &> {log}
        """
