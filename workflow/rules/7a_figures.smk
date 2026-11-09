# =============================================================================
# 7a_figures.smk — Final publication figures
# =============================================================================
#
# One rule per figure. Each rule reads finalized upstream tables and renders a
# complete multi-panel figure; every drawing decision (data columns, filters,
# labels, colours, panel sizes, canvas) lives in that figure's own script under
# workflow/scripts/figures/ and nowhere else. Nothing here imports
# workflow/src, so changing a figure cannot change another stage's output, and
# nothing here recomputes an upstream number: the inputs are the tables the
# analysis stages already produced.
#
# docs/figures.md describes each panel of each figure — inputs, transforms,
# drawing calls and layout constants — so the figure can be reproduced or
# retuned from that document alone.
#
# Figure 1: a (schematic, external vector PDF) + b (LFC coverage by FYPO class,
# Traditional Tn-Seq vs DIT-HAP) + c (spike-in dilution linearity).
# Figure S1: a PBL vs PBR + b technical replicate + c biological replicate, the
# first three panels of the 1a_pcr_qc QC figure re-rendered for print.

_FIGURES = "results/7a_figures"


rule plot_figure1:
    input:
        # Cross-dataset LFC coverage table, built by 2a_coverage.smk's
        # compare_lfc_coverage from the three released libraries.
        lfc_comparison="results/2a_coverage/cross_dataset/lfc_coverage_comparison.tsv",
        # Spike-in dilution series, built by 1a_pcr_qc.smk.
        spikein="results/1a_pcr_qc/spike_in_stats.tsv",
        # Panel A is drawn outside the pipeline; the page box is the content box.
        schematic="resources/figures/figure1a.pdf",
    output:
        pdf=f"{_FIGURES}/figure1.pdf",
        review=f"{_FIGURES}/figure1.review.png",
        # The canvas before the vector overlay, kept for layout debugging: it is
        # the file to open when a panel looks misplaced, since it shows the
        # reserved panel-A box instead of the schematic.
        canvas=f"{_FIGURES}/_work/figure1_canvas.pdf",
    log:
        "logs/7a_figures/plot_figure1.log",
    conda:
        "../envs/final_figures.yml"
    message:
        "*** [figures] Assembling Figure 1..."
    shell:
        """
        python workflow/scripts/figures/figure1.py \
            --lfc-comparison {input.lfc_comparison} \
            --spikein {input.spikein} \
            --schematic {input.schematic} \
            --output-pdf {output.pdf} \
            --output-review {output.review} \
            --canvas {output.canvas} &> {log}
        """


rule plot_figure_s1:
    input:
        # The 1a_pcr_qc stage publishes no table for panels a-c, only its 2x2
        # figure; these parquet intermediates are its prepared tables, so they
        # are the stable inputs here too.
        pbl_pbr="results/1a_pcr_qc/_work/pbl_pbr.parquet",
        tech="results/1a_pcr_qc/_work/tech.parquet",
        bio="results/1a_pcr_qc/_work/bio.parquet",
    output:
        pdf=f"{_FIGURES}/figure_s1.pdf",
        review=f"{_FIGURES}/figure_s1.review.png",
    log:
        "logs/7a_figures/plot_figure_s1.log",
    conda:
        "../envs/final_figures.yml"
    message:
        "*** [figures] Assembling Figure S1..."
    shell:
        """
        python workflow/scripts/figures/figure_s1.py \
            --pbl-pbr {input.pbl_pbr} \
            --tech {input.tech} \
            --bio {input.bio} \
            --output-pdf {output.pdf} \
            --output-review {output.review} &> {log}
        """
