#!/usr/local/bin

# USAGE   bash swarm.sh <folder with all demultiplexed fastas>

folder=$1

CENTROIDS="${folder}"/centroids
OUTPUTS="${folder}"/outputs
mkdir "${CENTROIDS}"
mkdir "${OUTPUTS}"
FILE1=("${folder}"/*fasta)

 n=0
for fasta_file in "${folder}"/*fasta; do
    n=$((n+1))
    echo -ne "Processing sample: ${n} of ${#FILE1[@]}\r"
    sample="${fasta_file%.*}"
    output=$(basename "${sample}")

    swarm  -t 4 -z -w "${CENTROIDS}"/"${output}".centroids.fasta "${fasta_file}" -d 2  -o "${OUTPUTS}"/"${output}".output.txt 2>> "${folder}/swarm.log"

done