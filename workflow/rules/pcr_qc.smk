# =============================================================================
# pcr_qc.smk — PCR / library-prep quality control figure (2x2 panels)
# =============================================================================
#
# Split into 2 rules so the loading/merging step and the figure-rendering step
# are independently re-runnable:
#   prepare_pcr_qc_data -> pbl_pbr / tech / bio / spikein parquet intermediates
#   plot_pcr_qc         -> the 2x2 QC figure PDF
# This module produces only a figure (no stats table), so the split follows
# the "load+merge" vs "render" boundary rather than verification.smk's
# "stats vs figure" boundary. NO dataset wildcard: this QC compares a few
# specifically-named libraries against each other (LD1328-7 processed twice,
# LD1328-4 vs LD1328-8), not a per-dataset generalization — same shape as the
# deferred spikein.smk. Ported from DIT_HAP_pipeline thesis_figures.ipynb
# ("2. PCR quality control"); see docs/plans/2026-07-19-pcr-qc-design.md.
#
# EXCEPTION to the release/ contract: panels (a)-(c) read upstream *pre-release*
# intermediates (results/8_merged/...), reachable ONLY via merged_reads_path()
# for datasets that declare `results_dir` in datasets.yaml. Panel (d) now reads
# the live spike-in stats table produced by spikein.smk (compute_spikein_stats):
# results/spikein/spike_in_stats.tsv — the earlier placeholder is retired. This
# makes spikein an upstream dependency of pcr_qc (a real DAG edge), fulfilling
# the "Phase 3+" interface promised in docs/plans/2026-07-19-pcr-qc-design.md §4.

import sys
sys.path.insert(0, workflow.basedir + "/..")  # repo root, so `workflow.src` imports resolve
from workflow.src.data_config import merged_reads_path

_PCR_QC = config["pcr_qc"]
_A = _PCR_QC["pbl_pbr"]
_B = _PCR_QC["technical_replicate"]
_C = _PCR_QC["biological_replicate"]

# Parquet intermediates shared between the two rules.
_PCRWORK = "results/pcr_qc/_work"


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
        # Panel (d): spike-in linearity — live output of spikein.smk (compute_spikein_stats).
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
        "../envs/statistics_and_figure_plotting.yml"
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
