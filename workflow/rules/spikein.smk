# =============================================================================
# spikein.smk — Spike-in dilution linearity QC
# =============================================================================
#
# Standalone (no {dataset} wildcard): reads the Spikein project's filtered
# insertion table, extracts 5 known spike-in coordinates, fits a log-log linear
# regression of reads vs known dilution ratios, and emits a stats TSV + PDF.
#
# EXCEPTION to the release/ contract (same shape as pcr_qc.smk): the filtered
# raw-reads table is a pre-release intermediate (results/13_filtered/...) that
# release/ never packages (see DIT_HAP_snakemake's packaging.smk RELEASE_MAP) —
# Spikein declares `results_dir` in datasets.yaml precisely so this rule can
# reach it.
#
# Split into 3 rules so each analysis step is independently re-runnable:
#   prepare_spikein_data      -> spike_in_stats parquet intermediate (the
#                                 single fan-out point)
#   compute_spikein_stats     -> long-form stats TSV
#   plot_spikein_correlation  -> log-log correlation PDF
# The two figure/table rules depend only on the prepared parquet, so editing
# e.g. the plot never forces the stats TSV to rebuild.

import json

# Parquet intermediate shared by the two downstream rules.
_SWORK = "results/spikein/_work"


rule prepare_spikein_data:
    input:
        raw_reads=(
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets']['Spikein']['results_dir']}/13_filtered/raw_reads.filtered.tsv"
        ),
    output:
        spike_in_stats=f"{_SWORK}/spike_in_stats.parquet",
    params:
        spike_in_sites_json=json.dumps(config.get("spikein", {}).get("coordinates", {})),
    log:
        "logs/spikein/prepare_spikein_data.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [spikein] Preparing spike-in stats..."
    shell:
        """
        python workflow/scripts/spikein/prepare_spikein_data.py \
            --raw-reads {input.raw_reads} \
            --output-spike-in-stats {output.spike_in_stats} \
            --spike-in-sites-json '{params.spike_in_sites_json}' &> {log}
        """


rule compute_spikein_stats:
    input:
        spike_in_stats=f"{_SWORK}/spike_in_stats.parquet",
    output:
        stats="results/spikein/spike_in_stats.tsv",
    log:
        "logs/spikein/compute_spikein_stats.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [spikein] Computing spike-in stats table..."
    shell:
        """
        python workflow/scripts/spikein/compute_spikein_stats.py \
            --spike-in-stats {input.spike_in_stats} \
            --output-stats {output.stats} &> {log}
        """


rule plot_spikein_correlation:
    input:
        spike_in_stats=f"{_SWORK}/spike_in_stats.parquet",
    output:
        figure="results/spikein/spike_in_correlation.pdf",
    log:
        "logs/spikein/plot_spikein_correlation.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [spikein] Plotting spike-in correlation..."
    shell:
        """
        python workflow/scripts/spikein/plot_spikein_correlation.py \
            --spike-in-stats {input.spike_in_stats} \
            --output-figure {output.figure} &> {log}
        """
