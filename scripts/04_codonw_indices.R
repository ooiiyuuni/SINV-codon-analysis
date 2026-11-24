# ============================================================
# 04_codonw_indices.R
# Function: Parse CodonW .out files and generate summary table
# Input: results/tables/codonw/{cds,nsp,sp}.out
# Output: results/tables/codonw/CodonUsageIndices_{CDS,NSP,SP}.csv
#         results/tables/codonw/CodonUsageIndices_Summary.csv
# Dependencies: tidyverse
# ============================================================

library(tidyverse)

in_dir  <- file.path("results", "tables", "codonw")
out_dir <- file.path("results", "tables", "codonw")

read_codonw_out <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  
  df <- df %>%
    select(title, A3s, T3s, G3s, C3s, GC3s, Nc, Gravy, Aromo) %>%
    rename(Accession = title, U3s = T3s)
  
  df$Accession <- trimws(df$Accession)
  df <- df[!is.na(df$Accession) & df$Accession != "", ]
  # Remove version numbers
  df$Accession <- sub("\\.\\d+$", "", df$Accession)
  df
}

append_mean_sd <- function(df, dataset_name) {
  num_cols <- c("A3s", "U3s", "G3s", "C3s", "GC3s", "Nc", "Gravy", "Aromo")
  
  mean_row <- df %>% summarise(across(all_of(num_cols), ~ mean(.x, na.rm = TRUE))) %>%
    mutate(Accession = paste0(dataset_name, "_Mean"), .before = 1)
  sd_row <- df %>% summarise(across(all_of(num_cols), ~ sd(.x, na.rm = TRUE))) %>%
    mutate(Accession = paste0(dataset_name, "_SD"), .before = 1)
  
  bind_rows(df, mean_row, sd_row)
}

format_4 <- function(df) {
  df %>% mutate(across(c(A3s, U3s, G3s, C3s, GC3s, Nc, Gravy, Aromo), ~ sprintf("%.4f", .x)))
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

datasets <- list(CDS = "cds.out", NSP = "nsp.out", SP = "sp.out")
summary_rows <- list()

for (name in names(datasets)) {
  file_path <- file.path(in_dir, datasets[[name]])
  if (!file.exists(file_path)) {
    cat(sprintf("Warning: %s not found, skipping.\n", file_path))
    next
  }
  
  df <- read_codonw_out(file_path)
  df_out <- append_mean_sd(df, name)
  
  write.csv(format_4(df_out), file.path(out_dir, paste0("CodonUsageIndices_", name, ".csv")), row.names = FALSE, quote = FALSE)
  
  mean_row <- df_out %>% filter(Accession == paste0(name, "_Mean"))
  sd_row   <- df_out %>% filter(Accession == paste0(name, "_SD"))
  
  summary_rows[[name]] <- tibble(
    Dataset = name,
    A3s   = sprintf("%.4f ± %.4f", mean_row$A3s,   sd_row$A3s),
    U3s   = sprintf("%.4f ± %.4f", mean_row$U3s,   sd_row$U3s),
    G3s   = sprintf("%.4f ± %.4f", mean_row$G3s,   sd_row$G3s),
    C3s   = sprintf("%.4f ± %.4f", mean_row$C3s,   sd_row$C3s),
    GC3s  = sprintf("%.4f ± %.4f", mean_row$GC3s,  sd_row$GC3s),
    Nc    = sprintf("%.4f ± %.4f", mean_row$Nc,    sd_row$Nc),
    Gravy = sprintf("%.4f ± %.4f", mean_row$Gravy, sd_row$Gravy),
    Aromo = sprintf("%.4f ± %.4f", mean_row$Aromo, sd_row$Aromo)
  )
}

if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  write.csv(summary_df, file.path(out_dir, "CodonUsageIndices_Summary.csv"), row.names = FALSE, quote = FALSE)
  cat("\n--- CodonW Indices Summary ---\n")
  print(as.data.frame(summary_df))
} else {
  cat("\nNo .out files found. Please run CodonW first and place files in results/tables/codonw/\n")
}
