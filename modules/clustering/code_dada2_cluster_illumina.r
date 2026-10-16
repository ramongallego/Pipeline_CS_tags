# code_dada2_cluster_illumina.r
# ILLUMINA mode counterpart of code_dada2_cluster.r
#
# usage: Rscript code_dada2_cluster_illumina.r <output_dir> <noprimers_dir> <USE_HASH> <LENR1> <LENR2>
#
# Differences from the CS version:
#   - one orientation only (<sample>_Locus_<locus>_Fwd.R1/.R2.fastq), so two
#     error models (R1, R2) instead of four, and no reorienting / joining step
#   - samples are matched to the sample map by fastq_header AND locus
#   - insect is not needed

arguments <- commandArgs(TRUE)

params <- list(folder = arguments[1],
               fastqs = arguments[2],
               hash   = arguments[3],
               len1   = as.numeric(arguments[4]),
               len2   = as.numeric(arguments[5]))
params

suppressPackageStartupMessages({
  library (tidyverse)
  library (dada2)
  library (digest)
  library (eDNAfuns)
  library (tictoc)
})

# dada() returns a single object instead of a list when given one sample
as_list <- function(x) if (is.list(x) && !inherits(x, c("dada", "derep", "data.frame"))) x else list(x)

sample.map <- read_delim(file.path(params$folder, "sample_trans.tmp"),
                         delim = "\t", show_col_types = FALSE)
print("The sample map looks like this")
head(sample.map)

## ---- listing files ---------------------------------------------------------
files.noprimers <- tibble(files = list.files(params$fastqs, pattern = "_Fwd\\.R[12]\\.fastq$", full.names = TRUE)) |>
  mutate(locus        = str_extract(basename(files), "(?<=_Locus_)[^_]+"),
         read         = str_extract(basename(files), "R[12](?=\\.fastq$)"),
         fastq_header = str_extract(basename(files), "^.+(?=_Locus_)")) |>
  pivot_wider(names_from = read, values_from = files) |>
  drop_na() |>
  inner_join(sample.map, by = c("fastq_header", "locus"))

print("Files matched to samples (one row per sample x locus)")
files.noprimers

if (nrow(files.noprimers) == 0) stop("No primer-trimmed files matched the sample map.")

## ---- filter and trim -------------------------------------------------------
filt_path <- file.path(params$folder, "filtered")
dir.create(filt_path, showWarnings = FALSE)

filt_function <- function(file1, file2) {
  filterAndTrim(file1, file.path(filt_path, basename(file1)),
                file2, file.path(filt_path, basename(file2)),
                truncLen = c(params$len1, params$len2),
                maxN = 0, maxEE = c(2, 2),
                truncQ = 2, rm.phix = TRUE,
                compress = TRUE, multithread = TRUE) |>
    as_tibble()
}

tic("filtering")
files.noprimers |>
  mutate(filt1 = file.path(filt_path, basename(R1)),
         filt2 = file.path(filt_path, basename(R2)),
         out   = map2_dfr(R1, R2, filt_function)) -> files.noprimers
toc()

# discard those with 100 reads or fewer after filtering
goodqs    <- files.noprimers |> filter(out$reads.out > 100)
discarded <- files.noprimers |> anti_join(goodqs |> select(fastq_header, locus),
                                          by = c("fastq_header", "locus"))
rm(files.noprimers)

if (nrow(goodqs) == 0) stop("No sample kept more than 100 reads after filtering - check LENR1/LENR2.")

## ---- learning errors (one model per read, shared by all samples) -----------
tic("Learning errors")
err1 <- learnErrors(goodqs$filt1, multithread = TRUE, verbose = 0, nbases = 100e6)
err2 <- learnErrors(goodqs$filt2, multithread = TRUE, verbose = 0, nbases = 100e6)
toc()
saveRDS(list(err1, err2), file = file.path(params$folder, "all.errors.rds"))

## ---- dereplication and denoising -------------------------------------------
tic("Dereplicating")
goodqs <- goodqs |>
  mutate(derep1 = map(filt1, derepFastq),
         derep2 = map(filt2, derepFastq))
toc()

tic("dada")
goodqs$dada_R1 <- as_list(dada(goodqs$derep1, err = err1, multithread = TRUE, verbose = FALSE, pool = "pseudo"))
goodqs$dada_R2 <- as_list(dada(goodqs$derep2, err = err2, multithread = TRUE, verbose = FALSE, pool = "pseudo"))
toc()

## ---- merging pairs and chimeras --------------------------------------------
tic("merging")
goodqs <- goodqs |>
  mutate(mergers = pmap(list(dada_R1, derep1, dada_R2, derep2), mergePairs, minOverlap = 10),
         nochim  = map(mergers, ~ if (nrow(.x) == 0) .x else removeBimeraDenovo(.x, multithread = TRUE)))
toc()

## ---- tidying and writing ---------------------------------------------------
goodqs |>
  select(Sample, locus, nochim) |>
  unnest(nochim) |>
  select(Sample, locus, sequence, nReads = abundance) -> Abundance_table

write_rds(Abundance_table, file.path(params$folder, "Abundance_table.rds"))

if (grepl("yes", params$hash, ignore.case = TRUE)) {

  Abundance_table |>
    distinct(sequence) |>
    mutate(Hash = map_chr(sequence, digest::sha1)) -> Hash_key

  write_csv(Hash_key, file.path(params$folder, "Hash_key.csv"))
  eDNAfuns::fasta_writer(df = Hash_key,
                         sequence = sequence,
                         header = Hash,
                         file.out = file.path(params$folder, "Hash_key.fasta"))

  Abundance_table |>
    inner_join(Hash_key, by = "sequence") |>
    select(Sample, locus, Hash, nReads) |>
    write_csv(file.path(params$folder, "ASV_table.csv"))

} else {
  write_csv(Abundance_table, file.path(params$folder, "ASV_table.csv"))
}

## ---- summary ---------------------------------------------------------------
getN <- function(x) sum(getUniques(x))

summary.good <- goodqs |>
  transmute(Sample, fastq_header, locus,
            input      = out$reads.in,
            filtered   = out$reads.out,
            dada_R1    = map_dbl(dada_R1, getN),
            dada_R2    = map_dbl(dada_R2, getN),
            merged     = map_dbl(mergers, ~ if (nrow(.x) == 0) 0 else getN(.x)),
            nochim     = map_dbl(nochim, ~ if (nrow(.x) == 0) 0 else getN(.x)))

summary.all <- bind_rows(summary.good,
                         discarded |>
                           transmute(Sample, fastq_header, locus,
                                     input    = out$reads.in,
                                     filtered = out$reads.out))

write_csv(summary.all, file.path(params$folder, "summary_dada2.csv"))

# Long format, same steps as the CS route (no 2.demultiplexing, no 7.reorienting)
summary.all |>
  pivot_longer(c(filtered, dada_R1, dada_R2, merged, nochim), names_to = "name", values_to = "nReads") |>
  mutate(step = case_when(name == "filtered"           ~ "4.filtering",
                          name %in% c("dada_R1", "dada_R2") ~ "5.dadaing",
                          name == "merged"             ~ "6.merging",
                          name == "nochim"             ~ "8.chimera"),
         direction = if_else(step == "8.chimera", "both", "Fwd")) |>
  replace_na(list(nReads = 0)) |>
  # for dada_R1 / dada_R2 keep the larger count, as in the CS route
  group_by(Sample, fastq_header, locus, direction, step) |>
  summarise(nReads = max(nReads), .groups = "drop") |>
  bind_rows(read_csv(file.path(params$folder, "summary.csv"), show_col_types = FALSE)) -> summary

write_csv(summary, file.path(params$folder, "pipeline_summary.csv"))

summary |>
  ggplot(aes(x = step, y = nReads, fill = direction)) +
  geom_col() -> p
ggsave(file.path(params$folder, "pipeline_summary.png"), p, width = 8, height = 6)
