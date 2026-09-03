#!/usr/bin/env bash
# =============================================================================
# fetch_sgd_data.sh — download the SGD tables used by gene annotation
# =============================================================================
#
# Downloads into resources/external/sgd/<version>/:
#   SGD_features.tab    per-ORF systematic name, standard name, qualifier, description
#   phenotype_data.tab  per-allele phenotype records (null-mutant viability)
#
# Driven by `rule fetch_sgd_data` (annotate.smk) — the only rule in this repo that
# makes a network call, and not reachable from `rule all`. Like pombase/ and
# biogrid/, the SGD copies are git-ignored regenerable caches; the rule exists so
# build_annotation_reference can declare the snapshot directory as an input.
#
# Usage
# -----
#     snakemake --cores 2 resources/external/sgd/2026-08-11    # preferred
#     bash workflow/scripts/annotate/fetch_sgd_data.sh 2026-08-11
#     bash workflow/scripts/annotate/fetch_sgd_data.sh         # defaults to today
#
# Note: phenotype_data.tab is ~39 MB and SGD can be slow; --retry covers dropouts
# (the rule adds `retries: 2` on top).

set -euo pipefail

readonly SGD_BASE_URL="https://downloads.yeastgenome.org/curation"
readonly VERSION="${1:-$(date +%Y-%m-%d)}"
readonly REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
readonly TARGET_DIR="${REPO_ROOT}/resources/external/sgd/${VERSION}"

mkdir -p "${TARGET_DIR}"

echo "Fetching SGD data into ${TARGET_DIR}"

# The progress meter is only useful interactively; under the fetch_sgd_data rule
# stderr is redirected to a log file, where it would add ~60 lines of carriage-
# returned noise per download. `--silent --show-error` rather than the tidier
# --no-progress-meter: the latter needs curl >= 7.67 and this host has 7.61.
# --show-error keeps real failures visible, which plain --silent would swallow.
PROGRESS_FLAG="--silent --show-error"
[[ -t 2 ]] && PROGRESS_FLAG=""
readonly PROGRESS_FLAG

# Download to a .part file then move into place, so an interrupted run cannot leave
# a truncated table looking complete. Deliberately NOT using --continue-at: SGD
# answers a range request past EOF with HTTP 416, which would make every re-run of
# an already-complete download fail.
fetch() {
    local url="$1" target="$2"
    curl --fail --location --retry 5 --retry-delay 5 ${PROGRESS_FLAG} \
        --output "${target}.part" "${url}"
    mv "${target}.part" "${target}"
    echo "  $(basename "${target}")  $(wc -l < "${target}") lines"
}

fetch "${SGD_BASE_URL}/chromosomal_feature/SGD_features.tab" "${TARGET_DIR}/SGD_features.tab"
fetch "${SGD_BASE_URL}/literature/phenotype_data.tab" "${TARGET_DIR}/phenotype_data.tab"

echo "Done. Build the annotation reference with:"
echo "  snakemake --use-conda --cores 4 results/annotation/2026-06-01/${VERSION}/gene_annotation_reference.parquet"
