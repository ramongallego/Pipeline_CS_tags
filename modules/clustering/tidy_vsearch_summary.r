# tidy_vsearch_summary.R
# Tidies vsearch and merge summary files into the same long format as summary.csv:
#   fastq_header, locus, step, direction, nReads
# Run after Parse_Abundances.R has finished writing ASV_table.csv

suppressPackageStartupMessages({library(tidyverse)})

arguments <- commandArgs(TRUE)
folder <- arguments[1]

summary_file <- file.path(folder, "summary.csv")
vsearch_file <- file.path(folder, "vsearch_summary.csv")
merge_file   <- file.path(folder, "merge_summary.csv")
new_summary_file <- file.path(folder, "pipeline_summary.csv")

# ── Helper: split sample column into fastq_header, locus, direction ──────────
# Handles two formats:
#   MEDITS_Lib1_A10_Locus_COI_Fwd → fastq_header=MEDITS_Lib1_A10, locus=COI, direction=Fwd
#   MEDITS_Lib1_E1_Locus_COI      → fastq_header=MEDITS_Lib1_E1,  locus=COI, direction=both

split_sample_col <- function(df, col = "Sample") {
  df |>
    separate(.data[[col]], 
             into = c("fastq_header", "rest"), 
             sep  = "_Locus_") |>
    separate(rest, 
             into = c("locus", "direction"), 
             sep  = "_", 
             extra = "drop",
             fill  = "right") |>
    mutate(direction = case_when(
      direction == "Fwd" ~ "Fwd",
      direction == "Rev" ~ "Rev",
      TRUE               ~ "both"
    ))
}

# ── vsearch_summary.csv ───────────────────────────────────────────────────────
# Long format: Sample, step, nReads
# step values (as written by vsearch_clustering.sh):
#   filtering, redirecting, denoising, chimeras
vsearch_summary <- read_csv(vsearch_file,
                             trim_ws   = TRUE, show_col_types = FALSE) |>
    split_sample_col()

# ── Append to existing summary.csv ───────────────────────────────────────────
read_csv(summary_file, show_col_types = FALSE) |>
  bind_rows(vsearch_summary) |>
  arrange(fastq_header, locus, direction, step) |>
  write_csv(new_summary_file)

 read_csv(new_summary_file, show_col_types = FALSE) |>
    ggplot(aes(x=step, y = nReads, fill= direction))+
    geom_col() -> p
ggsave(file.path(folder, "pipeline_summary.png"), p, width = 8, height = 6)
message("vsearch summaries tidied and appended to: ", new_summary_file)