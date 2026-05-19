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

### Authors
Made by Robert Cassidy and Ramón Gallego
