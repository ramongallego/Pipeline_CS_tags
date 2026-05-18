# modules/remove_primers.sh
#
# Builds per-library primer FASTA files, then runs cutadapt primer removal
# in two loops based on orientation encoded in the demux filename:
#
#   *_CS1.R1.fastq → CS1 hit R1 → Fwd primer on R1
#                    -g primers_Fwd.fasta -G primers_Rev.fasta
#
#   *_CS2.R1.fastq → CS2 hit R1 → Rev primer on R1
#                    -g primers_Rev.fasta -G primers_Fwd.fasta
#
# Output filenames have CS1/CS2 stripped via outbase, so the locus tag
# from cutadapt {name} lands cleanly:
#   e.g. LIB1_Adap_001_Locus_COI_Fwd.R1.fastq
#
# Originals are moved to DEMULT_DIR after trimming.
# Read counts are appended to OUTPUT_SUMMARY as step 3.Primers.
echo "############################################################################"
echo "Starting primer removal..."
echo ""
for (( i=0; i < NFILE1; i++ )); do

    ID1S=$(awk -F',' -v CF="$COLNUM_FILE1" -v VAL="${FILE1[i]}" \
                     -v ID="$COLNUM_ID1" \
               'NR>1 && $CF==VAL {print $ID; exit}' "${SEQUENCING_METADATA}")

    echo "  [$((i+1))/${NFILE1}] Library: ${ID1S}"

    DIR_LIB="${DEMULT_DIR}/${ID1S}"

    # ── Build primer FASTA files for this library ─────────────────────────────
    # primers_Fwd.fasta: one entry per locus, Fwd primer sequence
    # primers_Rev.fasta: one entry per locus, Rev primer sequence

    primers_Fwd="${BARCODES_DIR}/primers_${ID1S}_Fwd.fasta"
    primers_Rev="${BARCODES_DIR}/primers_${ID1S}_Rev.fasta"

   awk -F',' -v CF="$COLNUM_FILE1" -v VAL="${FILE1[i]}" \
          -v LOCUS="$COLNUM_LOCUS" -v FWD="$COLNUM_PRIMER1" \
    'NR>1 && $CF==VAL { print $LOCUS, $FWD }' \
    "${SEQUENCING_METADATA}" | sort -u | \
awk '{ printf ">Locus_%s_Fwd\n%s\n", $1, $2 }' > "${primers_Fwd}"

awk -F',' -v CF="$COLNUM_FILE1" -v VAL="${FILE1[i]}" \
          -v LOCUS="$COLNUM_LOCUS" -v REV="$COLNUM_PRIMER2" \
    'NR>1 && $CF==VAL { print $LOCUS, $REV }' \
    "${SEQUENCING_METADATA}" | sort -u | \
awk '{ printf ">Locus_%s_Rev\n%s\n", $1, $2 }' > "${primers_Rev}"

    # ── Loop 1: CS1 on R1 → Fwd primer on R1 ────────────────────────────────
    n=0
    for r1 in "${DIR_LIB}"/*CS1.R1.fastq; do
        [[ -e "$r1" ]] || continue
        n=$((n+1))
        echo -ne "    Fwd: sample ${n}\r"

        r2="${r1/R1.fastq/R2.fastq}"
        outbase=$(basename "${r1}" | sed 's/CS1\.R1\.fastq$//')

        cutadapt \
            -g "file:${primers_Fwd}" \
            -G "file:${primers_Rev}" \
            --discard-untrimmed \
            -o "${NOPRIMERS_DIR}/${outbase}{name}.R1.fastq" \
            -p "${NOPRIMERS_DIR}/${outbase}{name}.R2.fastq" \
            -j 4 --pair-adapters \
            "${r1}" "${r2}" \
            > "${SCRATCH_DIR}/primers_Fwd_${outbase}.log"

        # Per-locus trimmed counts
        awk -v S="${outbase}" '
            /^=== First read: Adapter/ { split($0,a," "); name=a[5] }
            /^Sequence:/ && name {
                for (i=1; i<=NF; i++) {
                    if ($i == "Trimmed:") { n=$(i+1); break }
                }
                print S","name",3.Primers,Fwd,"n
                name=""
            }
        ' "${SCRATCH_DIR}/primers_Fwd_${outbase}.log" >> "${OUTPUT_SUMMARY}"

        mv "${r1}" "${r2}" "${DEMULT_DIR}/"
    done
    echo

    # ── Loop 2: CS2 on R1 → Rev primer on R1 ────────────────────────────────
    n=0
    for r1 in "${DIR_LIB}"/*CS2.R1.fastq; do
        [[ -e "$r1" ]] || continue
        n=$((n+1))
        echo -ne "    Rev: sample ${n}\r"

        r2="${r1/R1.fastq/R2.fastq}"
        outbase=$(basename "${r1}" | sed 's/CS2\.R1\.fastq$//')

        cutadapt \
            -g "file:${primers_Rev}" \
            -G "file:${primers_Fwd}" \
            --discard-untrimmed \
            -o "${NOPRIMERS_DIR}/${outbase}{name}.R1.fastq" \
            -p "${NOPRIMERS_DIR}/${outbase}{name}.R2.fastq" \
            -j 4 --pair-adapters \
            "${r1}" "${r2}" \
            > "${SCRATCH_DIR}/primers_Rev_${outbase}.log"


        # Per-locus trimmed counts
        awk -v S="${outbase}" '
            /^=== First read: Adapter/ { split($0,a," "); name=a[5] }
            /^Sequence:/ && name {
                for (i=1; i<=NF; i++) {
                    if ($i == "Trimmed:") { n=$(i+1); break }
                }
                print S","name",3.Primers,Rev,"n
                name=""
            }
        ' "${SCRATCH_DIR}/primers_Rev_${outbase}.log" >> "${OUTPUT_SUMMARY}"

        mv "${r1}" "${r2}" "${DEMULT_DIR}/"
    done
    echo

    echo "    Primer removal done for ${ID1S}."
done

export NOPRIMERS_DIR
echo ""
echo "Primer removal complete."
echo "############################################################################"
echo ""