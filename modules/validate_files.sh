# modules/validate_files.sh
# Checks fastq pairs, secondary indexes, primers, and sample name uniqueness.
# Populates FILE1, FILE2, and SAMPLE_NAMES arrays; writes sample_trans.tmp.

echo "Validating input files and metadata..."

# ── Fastq pairs ───────────────────────────────────────────────────────────────
FILE1=($(awk -F',' -v COL="$COLNUM_FILE1" 'NR>1{print $COL}' "${SEQUENCING_METADATA}" | sort -u))
FILE2=($(awk -F',' -v COL="$COLNUM_FILE2" 'NR>1{print $COL}' "${SEQUENCING_METADATA}" | sort -u))
NFILE1="${#FILE1[@]}"
NFILE2="${#FILE2[@]}"

if [[ "$NFILE1" != "$NFILE2" ]]; then
    echo "ERROR: Mismatched R1/R2 file counts (${NFILE1} vs ${NFILE2}). Aborting."
    exit 1
fi

echo "  Found ${NFILE1} library pair(s):"
for (( i=0; i < NFILE1; i++ )); do
    if [[ ! -f "${PARENT_DIR}/${FILE1[i]}" ]]; then
        echo "  ERROR: R1 file not found: ${PARENT_DIR}/${FILE1[i]}"
        exit 1
    fi
    if [[ ! -f "${PARENT_DIR}/${FILE2[i]}" ]]; then
        echo "  ERROR: R2 file not found: ${PARENT_DIR}/${FILE2[i]}"
        exit 1
    fi
    printf '    R1: %s\n    R2: %s\n' "${FILE1[i]}" "${FILE2[i]}"
done

# ── Secondary indexes ─────────────────────────────────────────────────────────
if [[ "${SECONDARY_INDEX}" == "YES" ]]; then
    ID2S=($(  awk -F',' -v COL="$COLNUM_ID2"   'NR>1{print $COL}' "${SEQUENCING_METADATA}" | sort -u))
    ID2S_R=($(awk -F',' -v COL="$COLNUM_ID2_R" 'NR>1{print $COL}' "${SEQUENCING_METADATA}" | sort -u))

    if [[ "${#ID2S[@]}" -lt 2 ]]; then
        echo "ERROR: Only ${#ID2S[@]} forward secondary index found. Aborting."
        exit 1
    fi
    if [[ "${#ID2S_R[@]}" -lt 2 ]]; then
        echo "ERROR: Only ${#ID2S_R[@]} reverse secondary index found. Aborting."
        exit 1
    fi
    echo "  Fwd secondary indexes: ${#ID2S[@]} unique"
    echo "  Rev secondary indexes: ${#ID2S_R[@]} unique"
fi

# ── Primers ───────────────────────────────────────────────────────────────────
PRIMER1=($(awk -F',' -v COL="$COLNUM_PRIMER1" 'NR>1{print $COL}' "${SEQUENCING_METADATA}" | sort -u))
PRIMER2=($(awk -F',' -v COL="$COLNUM_PRIMER2" 'NR>1{print $COL}' "${SEQUENCING_METADATA}" | sort -u))
LOCUS=($(  awk -F',' -v COL="$COLNUM_LOCUS"   'NR>1{print $COL}' "${SEQUENCING_METADATA}" | sort -u))

if [[ -z "${PRIMER1[*]}" || -z "${PRIMER2[*]}" ]]; then
    echo "ERROR: Could not read primers from metadata. Aborting."
    exit 1
fi
echo "  Loci:        ${LOCUS[*]}"
echo "  Fwd primers: ${PRIMER1[*]}"
echo "  Rev primers: ${PRIMER2[*]}"

# ── Sample name uniqueness ────────────────────────────────────────────────────
SAMPLE_NAMES=($(awk -F',' -v COL="$COLNUM_SAMPLE" 'NR>1{print $COL}' "${SEQUENCING_METADATA}"))
UNIQ_SAMPLES=($( printf '%s\n' "${SAMPLE_NAMES[@]}" | sort -u))
if [[ "${#SAMPLE_NAMES[@]}" != "${#UNIQ_SAMPLES[@]}" ]]; then
    echo "ERROR: Duplicate sample names detected. Aborting."
    exit 1
fi
echo "  Sample names: ${#SAMPLE_NAMES[@]} unique ✓"

# ── Write sample translation file ────────────────────────────────────────────
# Columns (tab-separated): Sample | prefix_Fwd | prefix_Rev | ID_COMBO
# Used by R scripts to map demux filenames back to sample names.
SAMPLE_TRANS_FILE="${OUTPUT_DIR}/sample_trans.tmp"
printf "Sample\tfastq_header\tlocus\n" > "${SAMPLE_TRANS_FILE}"
awk -F',' \
    -v ID1="$COLNUM_ID1"      \
    -v ID2N="$COLNUM_ID2_NAME"\
    -v SN="$COLNUM_SAMPLE"    \
    -v LOCUS="$COLNUM_LOCUS"  \
    'NR>1 {
        id2 = $ID2N
        gsub(/_CS[12]/, "", id2)
        printf "%s\t%s_%s\t%s\n", $SN, $ID1, id2, $LOCUS
    }' \
    "${SEQUENCING_METADATA}" >> "${SAMPLE_TRANS_FILE}"

export FILE1 FILE2 NFILE1 SAMPLE_NAMES SAMPLE_TRANS_FILE

echo "Validation complete."
