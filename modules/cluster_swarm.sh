# modules/cluster_swarm.sh
# Runs the three-step Swarm secondary clustering:
#   1. prepare_for_swarm.R  — builds per-sample FASTA files from ASV table
#   2. swarm.sh             — runs Swarm on each FASTA
#   3. parsing_swarm.R      — parses centroids into final OTU table

echo "Running Swarm secondary clustering..."

echo "  [1/3] Preparing input FASTAs..."
Rscript "${MODULE_DIR}/clustering/prepare_for_swarm.R" "${OUTPUT_DIR}"
if [[ $? -ne 0 ]]; then
    echo "ERROR: prepare_for_swarm.R failed. Aborting." >&2
    exit 1
fi

echo "  [2/3] Running Swarm..."
bash "${MODULE_DIR}/clustering/swarm.sh" "${OUTPUT_DIR}/swarm_input"
if [[ $? -ne 0 ]]; then
    echo "ERROR: swarm.sh failed. Aborting." >&2
    exit 1
fi

echo "  [3/3] Parsing Swarm output..."
Rscript "${MODULE_DIR}/clustering/parsing_swarm.R" "${OUTPUT_DIR}"
if [[ $? -ne 0 ]]; then
    echo "ERROR: parsing_swarm.R failed. Aborting." >&2
    exit 1
fi

echo "Swarm clustering complete."
