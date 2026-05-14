#!/usr/bin/env bash
# run_pipeline.sh
# Usage: bash run_pipeline.sh banzai_params.sh

MAIN_DIR="$(dirname "$0")"
echo "Main dir is ${MAIN_DIR}"
LIBS_DIR="${MAIN_DIR}/lib"
MODULE_DIR="${MAIN_DIR}/modules"

# Load all helper functions
for file in "${LIBS_DIR}"/*.sh; do
    source "$file"
done

param_file="$1"
echo "Reading parameters from ${param_file}"

source "${MODULE_DIR}/read_params.sh"    "$param_file"
source "${MODULE_DIR}/check_metadata.sh"

echo "Creating output folder in ${OUTPUT_DIRECTORY}"
source "${MODULE_DIR}/init_folders.sh"
source "${MODULE_DIR}/extract_metadata_columns.sh"

echo "Initial setup complete. Proceeding to demultiplexing..."

if [[ "${ALREADY_DEMULTIPLEXED}" != "YES" ]]; then
    source "${MODULE_DIR}/validate_files.sh"
    source "${MODULE_DIR}/demultiplex.sh"
    source "${MODULE_DIR}/remove_primers.sh"
else
    echo "Skipping demultiplexing — reusing output from: ${DEMULT_OUTPUT}"
    cp "${DEMULT_OUTPUT}/sample_trans.tmp" "${OUTPUT_DIR}/"
    cp "${DEMULT_OUTPUT}/summary.csv"      "${OUTPUT_DIR}/"
    NOPRIMERS_DIR="${DEMULT_OUTPUT}/noprimers"
    export NOPRIMERS_DIR
fi

[[ "${SEARCH_ASVs}"       == "YES" ]] && source "${MODULE_DIR}/cluster_dada2.sh"
[[ "${SEARCH_Unoise}"     == "YES" ]] && source "${MODULE_DIR}/cluster_vsearch.sh"
[[ "${SECONDARY_SWARM}"   == "YES" ]] && source "${MODULE_DIR}/cluster_swarm.sh"

source "${MODULE_DIR}/hoard_control.sh"

echo "Pipeline complete. Output in: ${OUTPUT_DIR}"
