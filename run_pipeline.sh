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

# LIBRARY_PREP selects how samples were multiplexed:
#   CS_TAGS  (default) - CS1/CS2 tags inside the reads; demultiplex with cutadapt
#   ILLUMINA           - only i5/i7 indices; reads arrive already split per sample,
#                        Fwd primer on R1 and Rev primer on R2
# Older parameter files don't set it, so they keep running the CS_TAGS route.
LIBRARY_PREP="${LIBRARY_PREP:-CS_TAGS}"
if [[ "${LIBRARY_PREP}" != "CS_TAGS" && "${LIBRARY_PREP}" != "ILLUMINA" ]]; then
    echo "ERROR: LIBRARY_PREP must be CS_TAGS or ILLUMINA (got '${LIBRARY_PREP}'). Aborting."
    exit 1
fi
export LIBRARY_PREP
echo "Library prep: ${LIBRARY_PREP}"

source "${MODULE_DIR}/check_metadata.sh"

echo "Creating output folder in ${OUTPUT_DIRECTORY}"
source "${MODULE_DIR}/init_folders.sh"
if [[ "${LIBRARY_PREP}" == "ILLUMINA" ]]; then
    source "${MODULE_DIR}/extract_metadata_columns_illumina.sh"
else
    source "${MODULE_DIR}/extract_metadata_columns.sh"
fi

echo "Initial setup complete. Proceeding to demultiplexing..."

if [[ "${ALREADY_DEMULTIPLEXED}" != "YES" && "${LIBRARY_PREP}" == "ILLUMINA" ]]; then
    echo "Reads are already demultiplexed by the sequencer: going straight to primer removal."
    source "${MODULE_DIR}/validate_files_illumina.sh"
    source "${MODULE_DIR}/remove_primers_illumina.sh"
elif [[ "${ALREADY_DEMULTIPLEXED}" != "YES" ]]; then
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

if [[ "${SEARCH_ASVs}" == "YES" ]]; then
    if [[ "${LIBRARY_PREP}" == "ILLUMINA" ]]; then
        source "${MODULE_DIR}/cluster_dada2_illumina.sh"
    else
        source "${MODULE_DIR}/cluster_dada2.sh"
    fi
fi
[[ "${SEARCH_Unoise}"     == "YES" ]] && source "${MODULE_DIR}/cluster_vsearch.sh"
[[ "${SECONDARY_SWARM}"   == "YES" ]] && source "${MODULE_DIR}/cluster_swarm.sh"

source "${MODULE_DIR}/hoard_control.sh"

echo "Pipeline complete. Output in: ${OUTPUT_DIR}"
