#!/usr/bin/env bash
# =============================================================================
# fetch_sgd_data.sh — download the SGD tables used by gene annotation
# =============================================================================
#
# Run this ONCE by hand (not via Snakemake). Like resources/external/pombase and
# resources/external/biogrid, the SGD copies are regenerable local caches that
# are git-ignored and populated manually — this script exists to record how they
# were obtained, not to put a network call into the DAG.
#
# Downloads into resources/external/sgd/<date>/:
#   SGD_features.tab    per-ORF systematic name, standard name, qualifier, description
#   phenotype_data.tab  per-allele phenotype records (null-mutant viability)
#
# Usage
# -----
#     bash workflow/scripts/annotate/fetch_sgd_data.sh              # dated dir, e.g. 2026-08-11
#     bash workflow/scripts/annotate/fetch_sgd_data.sh 2026-08-11   # explicit version
#
# Note: phenotype_data.tab is ~39 MB and SGD can be slow; --retry covers dropouts.

set -euo pipefail

readonly SGD_BASE_URL="https://downloads.yeastgenome.org/curation"
readonly VERSION="${1:-$(date +%Y-%m-%d)}"
readonly REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
readonly TARGET_DIR="${REPO_ROOT}/resources/external/sgd/${VERSION}"

mkdir -p "${TARGET_DIR}"

echo "Fetching SGD data into ${TARGET_DIR}"

# Download to a .part file then move into place, so an interrupted run cannot leave
# a truncated table looking complete. Deliberately NOT using --continue-at: SGD
# answers a range request past EOF with HTTP 416, which would make every re-run of
# an already-complete download fail.
fetch() {
    local url="$1" target="$2"
    curl --fail --location --retry 5 --retry-delay 5 \
        --output "${target}.part" "${url}"
    mv "${target}.part" "${target}"
    echo "  $(basename "${target}")  $(wc -l < "${target}") lines"
}

fetch "${SGD_BASE_URL}/chromosomal_feature/SGD_features.tab" "${TARGET_DIR}/SGD_features.tab"
fetch "${SGD_BASE_URL}/literature/phenotype_data.tab" "${TARGET_DIR}/phenotype_data.tab"

echo "Done. Build the annotation reference with:"
echo "  snakemake --use-conda --cores 4 results/annotation/2026-06-01/${VERSION}/gene_annotation_reference.parquet"
