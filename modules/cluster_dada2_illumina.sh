# modules/cluster_dada2_illumina.sh
# ILLUMINA mode counterpart of cluster_dada2.sh.
# All reads are in one orientation (Fwd primer on R1), so the DADA2 script
# learns two error models (R1, R2) instead of four and skips the reorienting step.

echo "Running DADA2 clustering (ILLUMINA mode)..."

Rscript "${MODULE_DIR}/clustering/code_dada2_cluster_illumina.r" \
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
