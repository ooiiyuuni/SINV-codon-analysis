# ============================================================
# 12_prepare_web_input.R
# Purpose: Prepare CAIcal/RCDI web input files.
#          1. Clean virus CDS, NSP, SP by removing codons with 
#             ambiguous bases or stop codons.
#          2. Calculate host codon frequency and counts, and
#             format for CAIcal reference input.
# Input:   data/SINV/cds.fasta
#          data/SINV/nsp.fasta
#          data/SINV/sp.fasta
#          ~/Desktop/Culex pipiens cds.fna
#          ~/Desktop/Turdidae.fna
#          ~/Desktop/Passeriformes.fna
#          ~/Desktop/equine cds.fna
#          ~/Desktop/homo cds.fna
# Output:  results/new_data/virus/SINV_cds_cleaned.fasta
#          results/new_data/virus/SINV_nsp_cleaned.fasta
#          results/new_data/virus/SINV_sp_cleaned.fasta
#          results/new_data/host_freq/{species}_CAIcal.txt
# Dependencies: here, Biostrings, dplyr, tibble
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Biostrings)
  library(dplyr)
  library(tibble)
})

# ---------- 1. Directories ----------
out_virus_dir <- here("results", "new_data", "virus")
out_host_dir  <- here("results", "new_data", "host_freq")
dir.create(out_virus_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_host_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 2. Process Virus Sequences ----------
clean_sequences <- function(input_path, output_name) {
  message("Processing ", basename(input_path), "...")
  if (!file.exists(input_path)) {
    message("File not found: ", input_path, " skipping.")
    return(invisible(NULL))
  }
  
  dna <- readDNAStringSet(input_path)
  clean_list <- list()
  
  for (i in seq_along(dna)) {
    s <- toupper(as.character(dna[[i]]))
    name <- names(dna)[i]
    
    n_codons <- nchar(s) %/% 3
    if (n_codons == 0) next
    codons <- substring(s, seq(1, n_codons*3, 3), seq(3, n_codons*3, 3))
    
    bad_idx <- grepl("[^ACGT]", codons) | codons %in% c("TAA", "TAG", "TGA")
    
    if (all(bad_idx)) next
    clean_list[[name]] <- paste0(codons[!bad_idx], collapse = "")
  }
  
  out_dna <- DNAStringSet(unlist(clean_list))
  out_path <- file.path(out_virus_dir, output_name)
  writeXStringSet(out_dna, out_path)
  message(sprintf("  %s: %d retained out of %d", output_name, length(out_dna), length(dna)))
}

clean_sequences(here("data", "SINV", "cds.fasta"), "SINV_cds_cleaned.fasta")
clean_sequences(here("data", "SINV", "nsp.fasta"), "SINV_nsp_cleaned.fasta")
clean_sequences(here("data", "SINV", "sp.fasta"),  "SINV_sp_cleaned.fasta")

# ---------- 3. Process Host CDS ----------
message("Processing host frequencies...")
sense_codons <- c(
  "UUU","UUC","UUA","UUG","UCU","UCC","UCA","UCG",
  "UAU","UAC","UGU","UGC","UGG",
  "CUU","CUC","CUA","CUG","CCU","CCC","CCA","CCG",
  "CAU","CAC","CAA","CAG","CGU","CGC","CGA","CGG",
  "AUU","AUC","AUA","AUG","ACU","ACC","ACA","ACG",
  "AAU","AAC","AAA","AAG","AGU","AGC","AGA","AGG",
  "GUU","GUC","GUA","GUG","GCU","GCC","GCA","GCG",
  "GAU","GAC","GAA","GAG","GGU","GGC","GGA","GGG"
)
stop_codons <- c("UAA","UAG","UGA")
codon_order <- c(
  "UUU","UUC","UUA","UUG","UCU","UCC","UCA","UCG",
  "UAU","UAC","UAA","UAG","UGU","UGC","UGA","UGG",
  "CUU","CUC","CUA","CUG","CCU","CCC","CCA","CCG",
  "CAU","CAC","CAA","CAG","CGU","CGC","CGA","CGG",
  "AUU","AUC","AUA","AUG","ACU","ACC","ACA","ACG",
  "AAU","AAC","AAA","AAG","AGU","AGC","AGA","AGG",
  "GUU","GUC","GUA","GUG","GCU","GCC","GCA","GCG",
  "GAU","GAC","GAA","GAG","GGU","GGC","GGA","GGG"
)

calc_counts <- function(fasta_path) {
  dna <- readDNAStringSet(fasta_path)
  total <- 0L
  counts <- setNames(integer(length(sense_codons)), sense_codons)
  for (i in seq_along(dna)) {
    s <- gsub("T", "U", as.character(dna[[i]]))
    n <- nchar(s)
    if (n < 3) next
    starts <- seq(1, n - 2, by = 3)
    codons <- substring(s, starts, starts + 2)
    codons <- codons[!codons %in% stop_codons]
    if (length(codons) == 0) next
    tab <- table(codons)
    hit <- intersect(names(tab), sense_codons)
    counts[hit] <- counts[hit] + as.integer(tab[hit])
    total <- total + length(codons)
  }
  tibble(Codon = sense_codons,
         Count = as.integer(counts),
         Freq  = as.numeric(counts / total))
}

hosts <- tribble(
  ~species,         ~path,
  "Culex_pipiens",  "~/Desktop/Culex pipiens cds.fna",
  "Turdidae",       "~/Desktop/Turdidae.fna",
  "Passeriformes",  "~/Desktop/Passeriformes.fna",
  "Equus_caballus", "~/Desktop/equine cds.fna",
  "Homo_sapiens",   "~/Desktop/homo cds.fna"
)

for (i in seq_len(nrow(hosts))) {
  sp <- hosts$species[i]
  p <- path.expand(hosts$path[i])
  if (!file.exists(p)) stop("File not found: ", p)
  
  df <- calc_counts(p)
  
  df_full <- bind_rows(df, tibble(Codon = stop_codons, Count = 0L, Freq = 0)) %>%
    slice(match(codon_order, Codon)) %>%
    mutate(token = sprintf("%s %4.1f(%d)", Codon, Freq * 1000, Count))
  
  m <- matrix(df_full$token, nrow = 16, ncol = 4, byrow = TRUE)
  lines_out <- apply(m, 1, function(x) paste(x, collapse = "  "))
  
  out_file <- file.path(out_host_dir, paste0(sp, "_CAIcal.txt"))
  writeLines(lines_out, out_file)
  message(sprintf("[OK] Host: %s -> %s", sp, basename(out_file)))
}

message("Done. Outputs in results/new_data/")
