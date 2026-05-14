# modules/cluster_vsearch.sh
# Calls the VSEARCH/UNOISE clustering pipeline.
echo "############################################################################"
echo "Running VSEARCH/UNOISE clustering..."
echo ""

bash "${MODULE_DIR}/clustering/vsearch_clustering.sh" \
    "${OUTPUT_DIR}" \
    "${NOPRIMERS_DIR}" \
    "${USE_HASH}" \
    "${LENR1}" \
    "${LENR2}"

if [[ $? -ne 0 ]]; then
    echo "ERROR: VSEARCH clustering failed. Aborting." >&2
    exit 1
fi

echo "Tidying vsearch summary files..."
 
Rscript "${MODULE_DIR}/clustering/tidy_vsearch_summary.R" "${OUTPUT_DIR}"
 
if [[ $? -ne 0 ]]; then
    echo "ERROR: tidy_vsearch_summary.R failed. Aborting." >&2
    exit 1
fi
echo ""
echo "VSEARCH clustering complete."
echo "############################################################################"
echo ""