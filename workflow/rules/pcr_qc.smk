# =============================================================================
# pcr_qc.smk — Library-prep QC: spike-in linearity + PCR QC figure
# =============================================================================
#
# One stage, two halves (spikein.smk was folded in here 2026-09-16). The spike-in
# half produces the dilution-linearity table that the figure's panel (d) draws,
# so the two are wired as one module with a real DAG edge between them:
#
#   prepare_spikein_data   -> spike_in_stats parquet (single fan-out point for
#                             both downstream consumers)
#   compute_spikein_stats  -> long-form stats TSV (results/spikein/)
#   prepare_pcr_qc_data    -> pbl_pbr / tech / bio / spikein parquet intermediates
#   plot_pcr_qc            -> the 2x2 QC figure PDF
#
# The figure renders the spike-in stats itself (figure_render/spikein.py, panel
# d), so the former standalone plot_spikein_correlation rule and its
# results/spikein/spike_in_correlation.pdf are retired: one figure, not two.
#
# NO dataset wildcard: this QC compares a few specifically-named libraries
# against each other (LD1328-7 processed twice, LD1328-4 vs LD1328-8) plus the
# standalone Spikein calibration project — not a per-dataset generalization.
# Split per half so each analysis step is independently re-runnable; the
# figure's "load+merge" vs "render" boundary replaces verification.smk's
# "stats vs figure" split, since this module's only rendered artifact is a PDF.
# Ported from DIT_HAP_pipeline thesis_figures.ipynb ("2. PCR quality control")
# and spike_in.ipynb; see docs/plans/2026-07-19-pcr-qc-design.md.
#
# EXCEPTION to the release/ contract: spike-in reads the Spikein project's
# pre-release results/13_filtered/ table (release/ never packages it — see
# DIT_HAP_snakemake's packaging.smk RELEASE_MAP); panels (a)-(c) read upstream
# pre-release results/8_merged/ intermediates, reachable ONLY via
# merged_reads_path() for datasets that declare `results_dir` in datasets.yaml.

import json
import sys
sys.path.insert(0, workflow.basedir + "/..")  # repo root, so `workflow.src` imports resolve
from workflow.src.data_config import merged_reads_path

_PCR_QC = config["pcr_qc"]
_A = _PCR_QC["pbl_pbr"]
_B = _PCR_QC["technical_replicate"]
_C = _PCR_QC["biological_replicate"]

# Parquet intermediates.
_SPIKEWORK = "results/spikein/_work"
_PCRWORK = "results/pcr_qc/_work"


# ---------------------------------------------------------------------------
# Spike-in dilution linearity (upstream half)
# ---------------------------------------------------------------------------
rule prepare_spikein_data:
    input:
        raw_reads=(
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets']['Spikein']['results_dir']}/13_filtered/raw_reads.filtered.tsv"
        ),
    output:
        spike_in_stats=f"{_SPIKEWORK}/spike_in_stats.parquet",
    params:
        spike_in_sites_json=json.dumps(config.get("spikein", {}).get("coordinates", {})),
    log:
        "logs/spikein/prepare_spikein_data.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [pcr_qc] Preparing spike-in stats..."
    shell:
        """
        python workflow/scripts/spikein/prepare_spikein_data.py \
            --raw-reads {input.raw_reads} \
            --output-spike-in-stats {output.spike_in_stats} \
            --spike-in-sites-json '{params.spike_in_sites_json}' &> {log}
        """


rule compute_spikein_stats:
    input:
        spike_in_stats=f"{_SPIKEWORK}/spike_in_stats.parquet",
    output:
        stats="results/spikein/spike_in_stats.tsv",
    log:
        "logs/spikein/compute_spikein_stats.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [pcr_qc] Computing spike-in stats table..."
    shell:
        """
        python workflow/scripts/spikein/compute_spikein_stats.py \
            --spike-in-stats {input.spike_in_stats} \
            --output-stats {output.stats} &> {log}
        """


# ---------------------------------------------------------------------------
# Library-prep QC figure (downstream half)
# ---------------------------------------------------------------------------
rule prepare_pcr_qc_data:
    input:
        # Panel (a): PBL vs PBR of one library.
        pbl_pbr=merged_reads_path(_A["dataset"], _A["sample"], _A["timepoint"], _A["condition"]),
        # Panel (b): technical replicate — same sample in two upstream projects.
        tech_rep_1=merged_reads_path(_B["dataset_1"], _B["sample"], _B["timepoint"], _B["condition"]),
        tech_rep_2=merged_reads_path(_B["dataset_2"], _B["sample"], _B["timepoint"], _B["condition"]),
        # Panel (c): biological replicate — two samples in one project.
        bio_rep_1=merged_reads_path(_C["dataset"], _C["sample_1"], _C["timepoint"], _C["condition"]),
        bio_rep_2=merged_reads_path(_C["dataset"], _C["sample_2"], _C["timepoint"], _C["condition"]),
        # Panel (d): spike-in linearity — live output of compute_spikein_stats above.
        spikein="results/spikein/spike_in_stats.tsv",
    output:
        pbl_pbr=f"{_PCRWORK}/pbl_pbr.parquet",
        tech=f"{_PCRWORK}/tech.parquet",
        bio=f"{_PCRWORK}/bio.parquet",
        spikein=f"{_PCRWORK}/spikein.parquet",
    log:
        "logs/pcr_qc/prepare_pcr_qc_data.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [pcr_qc] Preparing merged tables..."
    shell:
        """
        python workflow/scripts/pcr_qc/prepare_pcr_qc_data.py \
            --pbl-pbr {input.pbl_pbr} \
            --tech-rep-1 {input.tech_rep_1} \
            --tech-rep-2 {input.tech_rep_2} \
            --bio-rep-1 {input.bio_rep_1} \
            --bio-rep-2 {input.bio_rep_2} \
            --spikein {input.spikein} \
            --output-pbl-pbr {output.pbl_pbr} \
            --output-tech {output.tech} \
            --output-bio {output.bio} \
            --output-spikein {output.spikein} &> {log}
        """


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
