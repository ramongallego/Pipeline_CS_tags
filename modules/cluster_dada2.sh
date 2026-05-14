# modules/cluster_dada2.sh
# Calls the DADA2 R script for ASV inference.

echo "Running DADA2 clustering..."

Rscript "${MODULE_DIR}/clustering/code_dada2_cluster.r" \
    "${OUTPUT_DIR}" \
    "${NOPRIMERS_DIR}" \
    "${USE_HASH}" \
    "${LENR1}" \
    "${LENR2}"

if [[ $? -ne 0 ]]; then
    echo "ERROR: DADA2 clustering failed. Aborting." >&2
    exit 1
fi

echo "DADA2 clustering complete."
