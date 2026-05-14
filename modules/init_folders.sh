# scripts/init_folders.sh
START_TIME=$(date +%Y%m%d_%H%M)
OUTPUT_DIR="${OUTPUT_DIRECTORY}/demultiplexed_${START_TIME}"
BARCODES_DIR="${OUTPUT_DIR}"/barcodes_and_primers
DEMULT_DIR="${OUTPUT_DIR}"/demultiplexed
NOPRIMERS_DIR="${OUTPUT_DIR}"/noprimers

mkdir -p "$OUTPUT_DIR"/{demultiplexed,noprimers,barcodes_and_primers}
cp "$SEQUENCING_METADATA" "$OUTPUT_DIR/metadata.csv"
cp "$param_file" "$OUTPUT_DIR/banzai_params.sh"

LOGFILE="$OUTPUT_DIR/logfile.txt"
exec > >(tee "$LOGFILE") 2>&1
