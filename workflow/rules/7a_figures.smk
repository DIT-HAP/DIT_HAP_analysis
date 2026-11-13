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
# Traditional Tn-Seq vs DIT-HAP) + c (spike-in dilution linearity) + d (gene-level
# LFC trajectories by viability class) + e (FACS-validated genes against DIT-HAP)
# + f (the weighted aggregation from insertion level to the gene-level curve,
# worked on cdc7).
# Figure 2: a the nrs1 example: its measured gene-level LFCs and the released
# Gompertz fit (A/DR/DL) with the formula written in the panel.
# Figure S1: a PBL vs PBR + b technical replicate + c biological replicate, the
# first three panels of the 1a_pcr_qc QC figure re-rendered for print, plus
# d the per-gene insertion count distribution of HD_DIT_HAP.

_FIGURES = "results/7a_figures"


rule plot_figure1:
    input:
        # Cross-dataset LFC coverage table, built by 2a_coverage.smk's
        # compare_lfc_coverage from the three released libraries.
        lfc_comparison="results/2a_coverage/cross_dataset/lfc_coverage_comparison.tsv",
        # Spike-in dilution series, built by 1a_pcr_qc.smk.
        spikein="results/1a_pcr_qc/spike_in_stats.tsv",
        # Panel D: the release's gene-level LFC table, one row per gene and one
        # column per generation (YES0..YES4).
        gene_lfc=(
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets']['HD_DIT_HAP']['release_dir']}/gene_level/LFC.tsv"
        ),
        # Panel E: FACS depletion curves measured outside the pipeline, one row
        # per gene x strain x timepoint.
        facs="resources/curated/FACS_depletion_curves_from_XH.xlsx",
        # Panel F: the release's insertion-level tables provide the per-site LFCs,
        # adjusted p-values (lollipop weights) and gene-timepoint weights; the
        # PomBase primary-transcript BED provides cdc7's CDS/intron model.
        insertion_lfc=(
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets']['HD_DIT_HAP']['release_dir']}/insertion_level/LFC.tsv"
        ),
        insertion_padj=(
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets']['HD_DIT_HAP']['release_dir']}/insertion_level/padj.tsv"
        ),
        insertion_weights=(
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets']['HD_DIT_HAP']['release_dir']}/insertion_level/"
            "insertion_weights.tsv"
        ),
        gene_intervals=(
            f"resources/external/pombase/{DATASETS['reference']['pombase_version']}/"
            "genome_region/coding_gene_primary_transcripts.bed"
        ),
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
            --gene-lfc {input.gene_lfc} \
            --facs {input.facs} \
            --insertion-lfc {input.insertion_lfc} \
            --insertion-padj {input.insertion_padj} \
            --insertion-weights {input.insertion_weights} \
            --gene-intervals {input.gene_intervals} \
            --schematic {input.schematic} \
            --output-pdf {output.pdf} \
            --output-review {output.review} \
            --canvas {output.canvas} &> {log}
        """


rule plot_figure2:
    input:
        # Panel A: the release's gene-level curve fits, one row per gene. The
        # fitted curve is drawn from the row's A/DR/DL; its *_fitted columns
        # are the cross-check that this script uses the same model.
        gene_fitting=(
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets']['HD_DIT_HAP']['release_dir']}/gene_level/"
            "gene_level_fitting_statistics.tsv"
        ),
    output:
        pdf=f"{_FIGURES}/figure2.pdf",
        review=f"{_FIGURES}/figure2.review.png",
        # The canvas before the A4 mount, kept for layout debugging.
        canvas=f"{_FIGURES}/_work/figure2_canvas.pdf",
    log:
        "logs/7a_figures/plot_figure2.log",
    conda:
        "../envs/final_figures.yml"
    message:
        "*** [figures] Assembling Figure 2..."
    shell:
        """
        python workflow/scripts/figures/figure2.py \
            --gene-fitting {input.gene_fitting} \
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
        # Panel D: the HD_DIT_HAP release's per-insertion annotations (only the
        # Type and Systematic ID columns are read).
        annotations=(
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets']['HD_DIT_HAP']['release_dir']}/insertion_level/annotations.tsv.gz"
        ),
    output:
        pdf=f"{_FIGURES}/figure_s1.pdf",
        review=f"{_FIGURES}/figure_s1.review.png",
        # The canvas before the A4 mount, kept for layout debugging.
        canvas=f"{_FIGURES}/_work/figure_s1_canvas.pdf",
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
            --annotations {input.annotations} \
            --output-pdf {output.pdf} \
            --output-review {output.review} \
            --canvas {output.canvas} &> {log}
        """
