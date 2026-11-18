# modules/remove_primers_illumina.sh
# ILLUMINA mode counterpart of remove_primers.sh.
#
# Reads come one pair of files per sample, always in the same orientation:
# Fwd primer at the start of R1, Rev primer at the start of R2. After cutting
# sequencing adapters (see pass 1 below), one cutadapt pass per sample,
# with one linked adapter per locus:
#
#   R1:  ^PrimerF ... revcomp(PrimerR)
#   R2:  ^PrimerR ... revcomp(PrimerF)
#
# The 5' primer is anchored (it must sit at the start of the read) and required.
# The 3' part is optional and trims read-through into the opposite primer and
# the sequencing adapter, which happens with short inserts and primer dimers.
# It needs at least 10 matching bases, so random matches at the end of normal
# reads are not trimmed.
#
# Pairs where either read ends up shorter than MIN_LENGTH (default 50) are
# dropped, which removes primer dimers.
#
# Output (same naming as the CS route, so the clustering scripts can read it):
#   ${NOPRIMERS_DIR}/<sample>_Locus_<locus>_Fwd.R1.fastq / .R2.fastq
# Read counts go to summary.csv as steps 1.Init and 3.Primers.

echo "############################################################################"
echo "Starting primer removal (ILLUMINA mode)..."
echo ""

MIN_LENGTH="${MIN_LENGTH:-50}"

OUTPUT_SUMMARY="${OUTPUT_DIR}/summary.csv"
printf "fastq_header,locus,step,direction,nReads\n" > "${OUTPUT_SUMMARY}"

SCRATCH_DIR="${OUTPUT_DIR}/_scratch"
mkdir -p "${SCRATCH_DIR}"

# cutadapt does not know inosine (I); it pairs with anything, so treat it as N
clean_primer() { echo "$1" | tr 'a-z' 'A-Z' | tr 'I' 'N'; }

count_reads() { [[ -s "$1" ]] && echo $(( $(wc -l < "$1") / 4 )) || echo 0; }

for (( i=0; i < NFILE1; i++ )); do

    READ1="${PARENT_DIR}/${FILE1[i]}"
    READ2="${PARENT_DIR}/${FILE2[i]}"

    SAMPLE=$(awk -F',' -v CF="$COLNUM_FILE1" -v VAL="${FILE1[i]}" -v S="$COLNUM_SAMPLE" \
                 'NR>1 && $CF==VAL {print $S; exit}' "${SEQUENCING_METADATA}")

    echo -ne "  [$((i+1))/${NFILE1}] ${SAMPLE}                    \r"

    # ── Linked primer FASTA files for this sample (one entry per locus) ──────
    # Both files list the loci in the same order: --pair-adapters pairs them by
    # position, and the R1 name becomes {name} in the output file names.
    primers_R1="${BARCODES_DIR}/primers_${SAMPLE}_R1.fasta"
    primers_R2="${BARCODES_DIR}/primers_${SAMPLE}_R2.fasta"
    : > "${primers_R1}"
    : > "${primers_R2}"

    while IFS=$'\t' read -r locus pF pR; do
        pF=$(clean_primer "$pF")
        pR=$(clean_primer "$pR")
        rcF=$(revcom "$pF")
        rcR=$(revcom "$pR")
        printf ">Locus_%s_Fwd\n^%s...%s;min_overlap=10\n" "$locus" "$pF" "$rcR" >> "${primers_R1}"
        printf ">Locus_%s_Fwd\n^%s...%s;min_overlap=10\n" "$locus" "$pR" "$rcF" >> "${primers_R2}"
    done < <(awk -F',' -v CF="$COLNUM_FILE1" -v VAL="${FILE1[i]}" \
                 -v L="$COLNUM_LOCUS" -v P1="$COLNUM_PRIMER1" -v P2="$COLNUM_PRIMER2" \
                 'NR>1 && $CF==VAL {print $L"\t"$P1"\t"$P2}' "${SEQUENCING_METADATA}" | sort -u)

    LOG_ADAPT="${SCRATCH_DIR}/adapters_${SAMPLE}.log"
    LOG="${SCRATCH_DIR}/primers_${SAMPLE}.log"

    # Pass 1: cut sequencing adapters (Nextera and TruSeq) from the 3' end.
    # Primer dimers often carry only part of the opposite primer before the
    # adapter, which the linked primer in pass 2 cannot recognise; once the
    # adapter is gone they are short enough for MIN_LENGTH to drop them.
    # Pass 2: primers, as described above. Reads go from one to the other
    # through a pipe (interleaved), so nothing extra is written to disk.
    # (with -o - cutadapt writes its report to stderr)
    cutadapt \
        -a CTGTCTCTTATACACATCT -a AGATCGGAAGAGC \
        -A CTGTCTCTTATACACATCT -A AGATCGGAAGAGC \
        -j 2 --interleaved -o - \
        "${READ1}" "${READ2}" 2> "${LOG_ADAPT}" |
    cutadapt \
        --interleaved \
        -a "file:${primers_R1}" \
        -A "file:${primers_R2}" \
        --pair-adapters \
        --discard-untrimmed \
        --minimum-length "${MIN_LENGTH}" \
        -j 2 \
        -o "${NOPRIMERS_DIR}/${SAMPLE}_{name}.R1.fastq" \
        -p "${NOPRIMERS_DIR}/${SAMPLE}_{name}.R2.fastq" \
        - > "${LOG}"
    status=("${PIPESTATUS[@]}")
    status_adapt="${status[0]}"
    status_primers="${status[1]}"

    if [[ "${status_adapt}" -ne 0 || "${status_primers}" -ne 0 ]]; then
        echo ""
        echo "ERROR: cutadapt failed for sample ${SAMPLE}. See ${LOG_ADAPT} and ${LOG}. Aborting." >&2
        exit 1
    fi

    # ── Read counts ───────────────────────────────────────────────────────────
    awk -v S="${SAMPLE}" '
        /^Total read pairs processed:/ { gsub(/,/, ""); print S",all_loci,1.Init,both,"$NF }
    ' "${LOG_ADAPT}" >> "${OUTPUT_SUMMARY}"

    # Counted from the output files: cutadapt reports linked adapters differently
    # from plain ones, and this way the number is what really goes downstream.
    awk -F',' -v CF="$COLNUM_FILE1" -v VAL="${FILE1[i]}" -v L="$COLNUM_LOCUS" \
        'NR>1 && $CF==VAL {print $L}' "${SEQUENCING_METADATA}" | sort -u |
    while read -r locus; do
        n=$(count_reads "${NOPRIMERS_DIR}/${SAMPLE}_Locus_${locus}_Fwd.R1.fastq")
        echo "${SAMPLE},${locus},3.Primers,Fwd,${n}" >> "${OUTPUT_SUMMARY}"
    done
done
echo ""

# cutadapt can leave empty files for loci with no reads; they would break clustering
find "${NOPRIMERS_DIR}" -name '*.fastq' -empty -delete

export NOPRIMERS_DIR SCRATCH_DIR OUTPUT_SUMMARY
echo ""
echo "Primer removal complete."
echo "############################################################################"
echo ""
