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

# Plotting rules will be migrated here in subsequent tickets (tickets 03-04: pure plotting script migration)
