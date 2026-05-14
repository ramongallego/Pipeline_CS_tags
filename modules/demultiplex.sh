# modules/demultiplex.sh
#
# Runs ONE cutadapt demultiplexing pass per library using combined barcode files:
#
#   barcodes_F.fasta — all CS1 entries (NAME=CS1, SEQ=CS1_seq)
#                    + all CS2 entries (NAME=CS2, SEQ=CS2_seq)
#
#   barcodes_R.fasta — same sequences, names swapped:
#                      NAME=CS1 carries CS2_seq
#                      NAME=CS2 carries CS1_seq
#
# With --pair-adapters, cutadapt matches R1 against barcodes_F and R2 against
# barcodes_R BY POSITION (same adapter name must hit both reads). This resolves
# both orientations in one pass:
#   Standard: CS1 on R1 → "CS1" hit in barcodes_F; CS2 on R2 → "CS1" entry
#             in barcodes_R carries CS2_seq → match ✓
#   Flipped:  CS2 on R1 → "CS2" hit in barcodes_F; CS1 on R2 → "CS2" entry
#             in barcodes_R carries CS1_seq → match ✓
#
# Output files land in DEMULT_DIR/${ID1S}/ and are named by the R1 barcode hit,
# which encodes the orientation (_CS1 = standard, _CS2 = flipped).
# Demux read counts are appended to OUTPUT_SUMMARY.
echo "############################################################################"
echo "Starting demultiplexing..."
echo ""

OUTPUT_SUMMARY="${OUTPUT_DIR}/summary.csv"
# Update header
printf "fastq_header,locus,step,direction,nReads\n" > "${OUTPUT_SUMMARY}"

SCRATCH_DIR="${OUTPUT_DIR}/_scratch"
mkdir -p "${SCRATCH_DIR}"

for (( i=0; i < NFILE1; i++ )); do

    READ1="${PARENT_DIR}/${FILE1[i]}"
    READ2="${PARENT_DIR}/${FILE2[i]}"

    ID1S=$(awk -F',' -v CF="$COLNUM_FILE1" -v VAL="${FILE1[i]}" \
                     -v ID="$COLNUM_ID1" \
               'NR>1 && $CF==VAL {print $ID; exit}' "${SEQUENCING_METADATA}")

    echo "  [$((i+1))/${NFILE1}] Library: ${ID1S}"

    # ── Build combined barcode FASTA files ────────────────────────────────────
    #
    # barcodes_F: each sample contributes two entries, its own name and sequence:
    #   >NAME_CS1  SEQ_CS1
    #   >NAME_CS2  SEQ_CS2
    #
    # barcodes_R: same sequences, names swapped — cutadapt's --pair-adapters
    # requires the SAME name to match on both R1 and R2, so swapping the names
    # here means a CS1 hit on R1 will look for a CS1-named entry on R2, which
    # carries SEQ_CS2 → correctly matches the opposite end:
    #   >NAME_CS1  SEQ_CS2   (swapped)
    #   >NAME_CS2  SEQ_CS1   (swapped)

    Barcodes_F="${BARCODES_DIR}/barcodes_${ID1S}_F.fasta"
    Barcodes_R="${BARCODES_DIR}/barcodes_${ID1S}_R.fasta"

    awk -F',' -v CF="$COLNUM_FILE1" -v VAL="${FILE1[i]}" \
              -v NAME1="$COLNUM_ID2_NAME"   -v SEQ1="$COLNUM_ID2" \
              -v NAME2="$COLNUM_ID2_NAME_R" -v SEQ2="$COLNUM_ID2_R" \
        'NR>1 && $CF==VAL {
            printf ">%s\n%s\n", $NAME1, $SEQ1
            printf ">%s\n%s\n", $NAME2, $SEQ2
        }' "${SEQUENCING_METADATA}" > "${Barcodes_F}"

    awk -F',' -v CF="$COLNUM_FILE1" -v VAL="${FILE1[i]}" \
              -v NAME1="$COLNUM_ID2_NAME"   -v SEQ1="$COLNUM_ID2" \
              -v NAME2="$COLNUM_ID2_NAME_R" -v SEQ2="$COLNUM_ID2_R" \
        'NR>1 && $CF==VAL {
            printf ">%s\n%s\n", $NAME1, $SEQ2
            printf ">%s\n%s\n", $NAME2, $SEQ1
        }' "${SEQUENCING_METADATA}" > "${Barcodes_R}"

    # ── Output dir for this library ───────────────────────────────────────────
    DIR_LIB="${DEMULT_DIR}/${ID1S}"
    mkdir -p "${DIR_LIB}"

    # ── Single cutadapt demultiplexing pass ───────────────────────────────────
    cutadapt \
        -g "file:${Barcodes_F};min_overlap=20" \
        -G "file:${Barcodes_R};min_overlap=20" \
        -o "${DIR_LIB}/${ID1S}_{name}.R1.fastq" \
        -p "${DIR_LIB}/${ID1S}_{name}.R2.fastq" \
        "${READ1}" "${READ2}" \
        --discard-untrimmed -j 6 -e 1 --pair-adapters \
        > "${SCRATCH_DIR}/demux_${ID1S}.log"

    # ── Append demux counts to summary ───────────────────────────────────────
 # Initial read count
awk -v LIB="${ID1S}" '
    /^Total read pairs processed:/ {
        gsub(/,/, "")
        print LIB",all_loci,1.Init,both,"$NF
    }
' "${SCRATCH_DIR}/demux_${ID1S}.log" >> "${OUTPUT_SUMMARY}"

# Per-adapter demux counts
awk -v LIB="${ID1S}" '
    /^=== First read: Adapter/ {
        split($0, a, " "); name=a[5]
        direction = (name ~ /CS1/) ? "Fwd" : (name ~ /CS2/) ? "Rev" : "unknown"
    }
    /^Sequence:/ && name {
        for (i=1; i<=NF; i++) {
            if ($i == "Trimmed:") { n=$(i+1); break }
        }
        print LIB"_"name",all_loci,2.demultiplexing,"direction","n
        name=""
    }
' "${SCRATCH_DIR}/demux_${ID1S}.log" >> "${OUTPUT_SUMMARY}"

    echo "    Demultiplexing done for ${ID1S}."
done

export SCRATCH_DIR OUTPUT_SUMMARY
echo ""
echo "Demultiplexing complete."
echo "###########################################################################"
echo ""