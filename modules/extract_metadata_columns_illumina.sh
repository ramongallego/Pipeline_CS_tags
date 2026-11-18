# modules/extract_metadata_columns_illumina.sh
# ILLUMINA mode counterpart of extract_metadata_columns.sh.
# Only the columns this mode needs are required: no CS tag or library columns.
# Exports the same COLNUM_* names as the CS module so later steps can share them.

echo "Resolving metadata column numbers (ILLUMINA mode)..."

COLNUM_FILE1=$(   get_colnum "${COLNAME_FILE1}"     "${SEQUENCING_METADATA}")
COLNUM_FILE2=$(   get_colnum "${COLNAME_FILE2}"     "${SEQUENCING_METADATA}")
COLNUM_SAMPLE=$(  get_colnum "${COLNAME_SAMPLE_ID}" "${SEQUENCING_METADATA}")
COLNUM_PRIMER1=$( get_colnum "${COLNAME_PRIMER1}"   "${SEQUENCING_METADATA}")
COLNUM_PRIMER2=$( get_colnum "${COLNAME_PRIMER2}"   "${SEQUENCING_METADATA}")
COLNUM_LOCUS=$(   get_colnum "${COLNAME_LOCUS}"     "${SEQUENCING_METADATA}")

# pairs of "COLNUM variable : COLNAME variable" so errors can name the missing column
required=(
    COLNUM_FILE1:COLNAME_FILE1
    COLNUM_FILE2:COLNAME_FILE2
    COLNUM_SAMPLE:COLNAME_SAMPLE_ID
    COLNUM_PRIMER1:COLNAME_PRIMER1
    COLNUM_PRIMER2:COLNAME_PRIMER2
    COLNUM_LOCUS:COLNAME_LOCUS
)

all_ok=true
for pair in "${required[@]}"; do
    num_var="${pair%%:*}"
    name_var="${pair##*:}"
    if [[ "${!num_var}" =~ ^[0-9]+$ && "${!num_var}" -gt 0 ]]; then
        echo "  OK: ${!name_var} is column ${!num_var}"
    else
        echo "  ERROR: column '${!name_var}' (${name_var}) not found in ${SEQUENCING_METADATA}"
        all_ok=false
    fi
done

if [[ "$all_ok" == false ]]; then
    echo "One or more required metadata columns are missing. Aborting."
    exit 1
fi

export COLNUM_FILE1 COLNUM_FILE2 COLNUM_SAMPLE \
       COLNUM_PRIMER1 COLNUM_PRIMER2 COLNUM_LOCUS

echo "All metadata columns resolved successfully."
