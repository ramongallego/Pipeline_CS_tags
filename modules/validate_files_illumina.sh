# modules/validate_files_illumina.sh
# ILLUMINA mode counterpart of validate_files.sh.
#
# Here every R1/R2 pair already belongs to one sample (the sequencer split the
# reads by their i5/i7 indices), so there is nothing to demultiplex. A sample can
# still carry several loci: one metadata row per sample x locus, rows of the same
# sample pointing to the same pair of files.
#
# Checks:
#   - every file pair exists, and R1/R2 are paired row by row (not sorted apart)
#   - one sample <-> one file pair
#   - sample x locus is unique
#   - sample and locus names are safe to use in file names
#   - primers are filled in
# Populates FILE1, FILE2, NFILE1 and writes sample_trans.tmp
# (fastq_header = sample name, which is how the noprimers files are named).

echo "Validating input files and metadata (ILLUMINA mode)..."

# One line per metadata row: sample <tab> locus <tab> file1 <tab> file2 <tab> primerF <tab> primerR
ROWS=$(awk -F',' -v S="$COLNUM_SAMPLE" -v L="$COLNUM_LOCUS" \
                 -v F1="$COLNUM_FILE1" -v F2="$COLNUM_FILE2" \
                 -v P1="$COLNUM_PRIMER1" -v P2="$COLNUM_PRIMER2" \
           'NR>1 && NF>1 {print $S"\t"$L"\t"$F1"\t"$F2"\t"$P1"\t"$P2}' "${SEQUENCING_METADATA}")

if [[ -z "${ROWS}" ]]; then
    echo "ERROR: no sample rows found in ${SEQUENCING_METADATA}. Aborting."
    exit 1
fi

errors=$(printf '%s\n' "${ROWS}" | awk -F'\t' '
    {
        s=$1; l=$2; f1=$3; f2=$4
        if (s !~ /^[A-Za-z0-9.-][A-Za-z0-9._-]*$/) print "  sample name \"" s "\" has characters other than letters, digits, . _ -"
        if (s ~ /_Locus_/)                    print "  sample name \"" s "\" contains _Locus_ (used internally as a separator)"
        if (l !~ /^[A-Za-z0-9.-]+$/)            print "  locus \"" l "\" (sample " s ") must be letters, digits, . or - (no underscores)"
        if ($5 == "" || $6 == "")             print "  sample " s ": missing Fwd or Rev primer"
        if (f1 == "" || f2 == "")             print "  sample " s ": missing file1 or file2"
        if (f1 == f2)                         print "  sample " s ": file1 and file2 are the same file"
        if (seen[s"\t"l]++)                   print "  sample " s " / locus " l " appears more than once"
        pair = f1"\t"f2
        if ((s in s2pair) && s2pair[s] != pair) print "  sample " s " points to more than one pair of files"
        if ((pair in pair2s) && pair2s[pair] != s) print "  files " f1 " / " f2 " are used by samples " pair2s[pair] " and " s
        if ((f1 in r1pair) && r1pair[f1] != f2) print "  R1 file " f1 " is paired with more than one R2 file"
        s2pair[s] = pair; pair2s[pair] = s; r1pair[f1] = f2
    }')

if [[ -n "${errors}" ]]; then
    echo "ERROR: problems in the metadata file:"
    echo "${errors}"
    echo "Aborting."
    exit 1
fi

# ── Fastq pairs (kept together, one entry per sample) ─────────────────────────
FILE1=()
FILE2=()
while IFS=$'\t' read -r f1 f2; do
    FILE1+=("$f1")
    FILE2+=("$f2")
done < <(printf '%s\n' "${ROWS}" | awk -F'\t' '{print $3"\t"$4}' | sort -u)
NFILE1="${#FILE1[@]}"

echo "  Found ${NFILE1} sample file pair(s)"
missing=0
for (( i=0; i < NFILE1; i++ )); do
    for f in "${FILE1[i]}" "${FILE2[i]}"; do
        if [[ ! -f "${PARENT_DIR}/${f}" ]]; then
            echo "  ERROR: file not found: ${PARENT_DIR}/${f}"
            missing=$((missing+1))
        fi
    done
done
if [[ "${missing}" -gt 0 ]]; then
    echo "ERROR: ${missing} fastq file(s) missing. Aborting."
    exit 1
fi

LOCUS=($(printf '%s\n' "${ROWS}" | cut -f2 | sort -u))
echo "  Loci: ${LOCUS[*]}"
echo "  Samples: $(printf '%s\n' "${ROWS}" | cut -f1 | sort -u | wc -l) unique"

# ── Sample translation file ──────────────────────────────────────────────────
# Same columns as in CS mode: Sample | fastq_header | locus
SAMPLE_TRANS_FILE="${OUTPUT_DIR}/sample_trans.tmp"
printf "Sample\tfastq_header\tlocus\n" > "${SAMPLE_TRANS_FILE}"
printf '%s\n' "${ROWS}" | awk -F'\t' '{printf "%s\t%s\t%s\n", $1, $1, $2}' >> "${SAMPLE_TRANS_FILE}"

export NFILE1 SAMPLE_TRANS_FILE

echo "Validation complete."
