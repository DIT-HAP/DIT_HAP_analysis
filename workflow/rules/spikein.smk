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
# Single rule (no prepare/compute split — data is tiny and self-contained).

import json

rule run_spikein_analysis:
    input:
        raw_reads=(
            f"{DATASETS['snakemake_repo']}/"
            f"{DATASETS['datasets']['Spikein']['results_dir']}/13_filtered/raw_reads.filtered.tsv"
        ),
    output:
        stats="results/spikein/spike_in_stats.tsv",
        figure="results/spikein/spike_in_correlation.pdf",
    params:
        spike_in_sites_json=json.dumps(config.get("spikein", {}).get("coordinates", {})),
    log:
        "logs/spikein/run_spikein_analysis.log",
    conda:
        "../envs/statistics_and_figure_plotting.yml"
    message:
        "*** [spikein] Running spike-in linearity QC..."
    shell:
        """
        python workflow/scripts/spikein/run_spikein_analysis.py \
            --raw-reads {input.raw_reads} \
            --output-stats {output.stats} \
            --output-figure {output.figure} \
            --spike-in-sites-json '{params.spike_in_sites_json}' &> {log}
        """
