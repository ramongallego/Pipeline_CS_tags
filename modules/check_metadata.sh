# scripts/check_metadata.sh
fix_crlf() {
    local file="$1"
    if [[ $(file "$file") == *"CRLF"* ]]; then
        echo "Fixing CRLF endings in $file..."
        tr -d '\r' < "$file" > "${file%.csv}_fixed.csv"
        echo "${file%.csv}_fixed.csv"
    else
        echo "$file"
    fi
}

if [[ ! -s "$SEQUENCING_METADATA" ]]; then
    echo "Metadata file not found or empty: $SEQUENCING_METADATA"
    exit 1
fi

SEQUENCING_METADATA=$(fix_crlf "$SEQUENCING_METADATA")
export SEQUENCING_METADATA
