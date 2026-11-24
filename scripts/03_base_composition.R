# ============================================================
# 03_base_composition.R
# Function: Calculate base composition (A,U,G,C,GC,AUGC,Deviation)
# Input: data/SINV/{cds,nsp,sp}.fasta
# Output: results/tables/BaseComposition_{CDS,NSP,SP}.csv
#         results/tables/BaseComposition_Summary.csv
# Dependencies: here, seqinr, tidyverse
# ============================================================

library(seqinr)
library(tidyverse)

# Assuming you run this from the project root directory
read_sinv <- function(fname) {
  s <- read.fasta(file.path("data", "SINV", fname), seqtype = "DNA", as.string = TRUE)
  # Remove version numbers from sequence names
  names(s) <- sub("\\.\\d+$", "", names(s))
  s
}

cds <- read_sinv("cds.fasta")
nsp <- read_sinv("nsp.fasta")
sp  <- read_sinv("sp.fasta")

base_composition <- function(seq_list) {
  tibble(
    Accession = names(seq_list),
    Length    = sapply(seq_list, function(x) nchar(x[1])),
    A = sapply(seq_list, function(x) str_count(toupper(x[1]), "A")),
    U = sapply(seq_list, function(x) str_count(toupper(x[1]), "T")),
    G = sapply(seq_list, function(x) str_count(toupper(x[1]), "G")),
    C = sapply(seq_list, function(x) str_count(toupper(x[1]), "C"))
  ) %>%
    mutate(
      # Calculate valid length (ignoring N, R, Y, etc.)
      ValidLength = A + U + G + C,
      `A%`        = A / ValidLength * 100,
      `U%`        = U / ValidLength * 100,
      `G%`        = G / ValidLength * 100,
      `C%`        = C / ValidLength * 100,
      `GC%`       = (G + C) / ValidLength * 100,
      `AUGC%`     = (A + U + G + C) / Length * 100,
      `Deviation` = 100 - `AUGC%`
    ) %>%
    select(Accession, Length, `A%`, `U%`, `G%`, `C%`, `GC%`, `AUGC%`, `Deviation`)
}

append_summary <- function(df, file_name, dataset_name) {
  mean_row <- df %>% summarise(across(where(is.numeric), ~ mean(.x))) %>%
    mutate(Accession = paste0(dataset_name, "_Mean"), .before = 1)
  sd_row <- df %>% summarise(across(where(is.numeric), ~ sd(.x))) %>%
    mutate(Accession = paste0(dataset_name, "_SD"), .before = 1)

  out <- bind_rows(df, mean_row, sd_row) %>%
    mutate(
      Length = ifelse(grepl("_Mean|_SD", Accession), sprintf("%.1f", Length), sprintf("%.0f", Length)),
      across(c(`A%`, `U%`, `G%`, `C%`, `GC%`, `AUGC%`, `Deviation`), ~ sprintf("%.4f", .x))
    )
  
  out_dir <- file.path("results", "tables")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(out, file.path(out_dir, file_name), row.names = FALSE, quote = FALSE)
}

cat("Calculating base composition...\n")
append_summary(base_composition(cds), "BaseComposition_CDS.csv", "CDS")
append_summary(base_composition(nsp), "BaseComposition_NSP.csv", "NSP")
append_summary(base_composition(sp),  "BaseComposition_SP.csv",  "SP")

cat("Base composition saved to results/tables/\n")

# Print the CDS summary to console
cat("\n--- CDS Summary (from CSV) ---\n")
cds_out <- read.csv(file.path("results", "tables", "BaseComposition_CDS.csv"))
print(tail(cds_out, 2))
