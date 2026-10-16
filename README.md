# demulting_CS
We have now a different set of barcodes that we can add to our amplicons before adding library indices and sequence samples using Illumina Platforms. Because more and more Illumina sequencers are moving towards binned Qscores, I am adding a selection of clustering algorithms - either DADA2 or UNOISE, and the option to collapse those ASVs/OTUs with swarm *within* each sample. At present, all clustering parameters are hard coded into the script, besides the trimming lengths of the amplicons after primer removal. 
The system allows for multiple primer sets at the same time, which is not super usual but Ok. In this version, the trimming parameters are global, but in future versions we should allow for a primer-specific setup.

## Dependencies

This script has been tested in Windows using WSL2. It should work in Linux as well, and probably in a Mac - I have to check if all the dependencies are available for Mac

The script relies also on R packages. 

#### Stand alone dependencies

* cutadapt (https://cutadapt.readthedocs.io/en/stable/index.html).
* vsearch (https://github.com/torognes/vsearch)
* flash2 (https://github.com/dstreett/FLASH2)
* seqkit (https://github.com/shenwei356/seqkit)
* swarm (https://github.com/torognes/swarm)

### R packages

Checking whether you have the right R packages installed before running R is complicated - and you don't want to realise you don't have them by the time the pipeline gets there. So please check that you have these packages installed:

* tidyverse
* dada2
* Biostrings
* digest
* eDNAfuns
* tictoc

### Usage

I would start by creating a metadata file in which, for each sample you fill in the information the pipeline needs: fastq files, barcode indices (include the link tags in those), Fwd and Rev primers... Then, create a copy of the parameter file and fill in the study-specific parameters

Then run the job with bash <path_to_pipeline>/run_pipeline.sh <path_to_parameters_file>

### Libraries without CS tags (ILLUMINA mode)

If samples were multiplexed only with Illumina i5/i7 indices, the sequencer already delivers one R1/R2 pair per sample, with the Fwd primer on R1 and the Rev primer on R2. Set `LIBRARY_PREP="ILLUMINA"` in the parameters file (see `data/params_illumina_MiFish.sh`) and the pipeline skips demultiplexing and goes straight to primer removal, then to the clustering you choose. Parameter files without `LIBRARY_PREP` keep running the CS route exactly as before.

Metadata: one row per sample and locus, with at least these columns (names are set in the parameters file):

| column | example | notes |
|---|---|---|
| sample | `CWH_2603_01_R1` | letters, digits, `.` `_` `-`; sample + locus must be unique |
| PrimerF / PrimerR | `GYCGGTAAAACTCGTGCCAGC` | IUPAC codes allowed; inosine (I) is read as N |
| locus | `MiFish` | no underscores |
| file1 / file2 | `..._1.fq` / `..._2.fq` | relative to `PARENT_DIR`; one pair of files per sample |

Primer removal in this mode:
1. cuts Nextera and TruSeq adapter read-through from the 3' end;
2. removes the primers with one linked adapter per locus (`^PrimerF...revcomp(PrimerR)` on R1, `^PrimerR...revcomp(PrimerF)` on R2). Pairs without both primers at the start of the reads are discarded;
3. drops pairs where either read is shorter than `MIN_LENGTH` (default 50), which removes primer dimers.

DADA2 uses its own script in this mode (`code_dada2_cluster_illumina.r`), with two error models instead of four. UNOISE and swarm use the same scripts as the CS route.

### Output

Depending on the options chosen, the output would be either an ASV table or an OTU table. Either way, the basic output files are:
  * ASV_table.csv: Abundances tables in a long format, with columns `Hash`, `Sample` and `nReads`
  * Hash_key.[csv|fasta]: Representative sequences in a csv or fasta format
  * metadata.csv: needed to repeat the analysis
  * params.sh : ditto
  * summary.csv: number of sequences kept at the end of each step (demultiplexing and primer finding)
  * summary_[dada|unoise].csv: number of sequences at the end of each step of the clustering processes
  * pipeline_summary.jpg: a barplot with the number of sequences after each step (useful for troubleshooting)

If you choose secondary clustering with `swarm` it generates a new set of ASV table and Hash key files - not included in the summary file as swarm does not drop any reads
### Authors
Made by Robert Cassidy and Ramón Gallego
