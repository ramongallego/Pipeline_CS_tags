# modules/extract_metadata_columns.sh
# Resolves column names (defined in params) to column numbers in the metadata CSV.
# All COLNUM_* variables are exported for use by downstream modules.

echo "Resolving metadata column numbers..."

COLNUM_FILE1=$(      get_colnum "${COLNAME_FILE1}"       "${SEQUENCING_METADATA}")
COLNUM_FILE2=$(      get_colnum "${COLNAME_FILE2}"       "${SEQUENCING_METADATA}")
COLNUM_ID1=$(        get_colnum "${COLNAME_ID1_NAME}"    "${SEQUENCING_METADATA}")
COLNUM_ID2=$(        get_colnum "${COLNAME_ID2_SEQ}"     "${SEQUENCING_METADATA}")
COLNUM_ID2_R=$(      get_colnum "${COLNAME_ID2_SEQ_R}"   "${SEQUENCING_METADATA}")
COLNUM_ID2_NAME=$(   get_colnum "${COLNAME_ID2_NAME}"    "${SEQUENCING_METADATA}")
COLNUM_ID2_NAME_R=$( get_colnum "${COLNAME_ID2_NAME_R}"  "${SEQUENCING_METADATA}")
COLNUM_SAMPLE=$(     get_colnum "${COLNAME_SAMPLE_ID}"   "${SEQUENCING_METADATA}")
COLNUM_PRIMER1=$(    get_colnum "${COLNAME_PRIMER1}"     "${SEQUENCING_METADATA}")
COLNUM_PRIMER2=$(    get_colnum "${COLNAME_PRIMER2}"     "${SEQUENCING_METADATA}")
COLNUM_LOCUS=$(      get_colnum "${COLNAME_LOCUS}"       "${SEQUENCING_METADATA}")

all_columns=(
    COLNUM_FILE1 COLNUM_FILE2
    COLNUM_ID1   COLNUM_ID2    COLNUM_ID2_R
    COLNUM_ID2_NAME COLNUM_ID2_NAME_R
    COLNUM_SAMPLE
    COLNUM_PRIMER1 COLNUM_PRIMER2 COLNUM_LOCUS
)

all_ok=true
for col in "${all_columns[@]}"; do
    if [[ "${!col}" -gt 0 ]]; then
        echo -ne "  OK: ${col} = ${!col}\r"
    else
        echo "  ERROR: could not find column '${col}' (looked for name '${!col//COLNUM_/}')."
        all_ok=false
    fi
done

if [[ "$all_ok" == false ]]; then
    echo "One or more required metadata columns are missing. Aborting."
    exit 1
fi

export COLNUM_FILE1 COLNUM_FILE2 \
       COLNUM_ID1 COLNUM_ID1_SEQ \
       COLNUM_ID2 COLNUM_ID2_R \
       COLNUM_ID2_NAME COLNUM_ID2_NAME_R \
       COLNUM_SAMPLE \
       COLNUM_PRIMER1 COLNUM_PRIMER2 COLNUM_LOCUS

echo "All metadata columns resolved successfully."
