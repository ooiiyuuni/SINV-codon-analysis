# ============================================================
# 02_degenerate_bases.R
# Function: Extract detailed table of degenerate bases for CDS only.
# Input:  data/SINV/cds.fasta
# Output: results/tables/DegenerateBases_CDS.csv
# ============================================================

library(seqinr)
library(tidyverse)

in_dir  <- "data/SINV"
out_dir <- "results/tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fname <- "cds.fasta"
filepath <- file.path(in_dir, fname)
seqs <- read.fasta(filepath, seqtype = "DNA", as.string = TRUE)

res <- lapply(names(seqs), function(id) {
  seq_str <- tolower(as.character(seqs[[id]]))
  seq_vec <- strsplit(seq_str, "")[[1]]
  seq_len <- length(seq_vec)
  
  # Find positions of degenerate bases
  deg_idx <- which(!(seq_vec %in% c("a", "c", "g", "t")))
  
  if (length(deg_idx) == 0) return(NULL)
  
  symbols <- toupper(seq_vec[deg_idx])
  positions <- deg_idx
  
  codons <- unique(ceiling(deg_idx / 3))
  n_codons <- length(codons)
  
  pct_codons <- (n_codons / (seq_len / 3)) * 100
  
  # Remove version number from id
  clean_id <- sub("\\.\\d+$", "", id)
  
  data.frame(
    Accession = clean_id,
    `Length (bp)` = seq_len,
    `No. of degenerate bases` = length(deg_idx),
    `Degenerate bases` = paste(unique(symbols), collapse = "; "),
    Positions = paste(positions, collapse = "; "),
    `No. of affected codons` = n_codons,
    `% of codons` = sprintf("%.3f", pct_codons),
    check.names = FALSE
  )
})

df <- bind_rows(res)

if (!is.null(df) && nrow(df) > 0) {
  # Calculate Total row
  total_len <- sum(df$`Length (bp)`)
  total_deg <- sum(df$`No. of degenerate bases`)
  total_codons <- sum(df$`No. of affected codons`)
  total_pct <- (total_codons / (total_len / 3)) * 100
  
  total_row <- data.frame(
    Accession = "Total",
    `Length (bp)` = total_len,
    `No. of degenerate bases` = total_deg,
    `Degenerate bases` = "—",
    Positions = "—",
    `No. of affected codons` = total_codons,
    `% of codons` = sprintf("%.3f", total_pct),
    check.names = FALSE
  )
  
  df <- bind_rows(df, total_row)
  
  write.csv(df, file.path(out_dir, "DegenerateBases_CDS.csv"), row.names = FALSE, quote = FALSE)
  cat("Saved to results/tables/DegenerateBases_CDS.csv\n\n")
  
  # Print the table to console so you can see it
  print(as.data.frame(df))
} else {
  cat("No degenerate bases found in CDS.\n")
}
