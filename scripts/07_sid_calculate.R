# =============================================================================
# 07_sid_calculate.R
# Purpose: Calculate per-sequence SiD between SINV and hosts.
# Input: results/tables/RSCU/SINV/RSCU_long_*.csv
#        results/tables/RSCU/hosts/*_codon_freq.txt
# Output: results/tables/SiD/SiD_long.csv
# =============================================================================
suppressPackageStartupMessages({ library(here); library(tidyverse) })
set.seed(12345)

read_sinv_long <- function(path, region) {
  read_csv(path, show_col_types = FALSE) %>% select(Accession, Codon, RSCU) %>%
    filter(!Codon %in% c("AUG", "UGG")) %>% mutate(Region = region)
}
sinv <- bind_rows(
  read_sinv_long(here("results/tables/RSCU/SINV/RSCU_long_CDS.csv"), "CDS"),
  read_sinv_long(here("results/tables/RSCU/SINV/RSCU_long_NSP.csv"), "NSP"),
  read_sinv_long(here("results/tables/RSCU/SINV/RSCU_long_SP.csv"),  "SP")
)

host_files <- c("Culex pipiens" = "Culex_pipiens_codon_freq.txt", "Turdidae" = "Turdidae_codon_freq.txt",
                "Passeriformes" = "Passeriformes_codon_freq.txt", "Equus caballus" = "Equus_caballus_codon_freq.txt",
                "Homo sapiens" = "Homo_sapiens_codon_freq.txt")
codon_aa <- c(UUU="Phe", UUC="Phe", UUA="Leu", UUG="Leu", CUU="Leu", CUC="Leu", CUA="Leu", CUG="Leu", AUU="Ile", AUC="Ile", AUA="Ile", AUG="Met", GUU="Val", GUC="Val", GUA="Val", GUG="Val", UCU="Ser", UCC="Ser", UCA="Ser", UCG="Ser", CCU="Pro", CCC="Pro", CCA="Pro", CCG="Pro", ACU="Thr", ACC="Thr", ACA="Thr", ACG="Thr", GCU="Ala", GCC="Ala", GCA="Ala", GCG="Ala", UAU="Tyr", UAC="Tyr", UAA="Stop", UAG="Stop", CAU="His", CAC="His", CAA="Gln", CAG="Gln", AAU="Asn", AAC="Asn", AAA="Lys", AAG="Lys", GAU="Asp", GAC="Asp", GAA="Glu", GAG="Glu", UGU="Cys", UGC="Cys", UGA="Stop", UGG="Trp", CGU="Arg", CGC="Arg", CGA="Arg", CGG="Arg", AGU="Ser", AGC="Ser", AGA="Arg", AGG="Arg", GGU="Gly", GGC="Gly", GGA="Gly", GGG="Gly")

read_host_rscu <- function(fname, host_name) {
  read_tsv(here("results/tables/RSCU/hosts", fname), col_names = c("Codon", "Freq"), show_col_types = FALSE) %>%
    filter(!Codon %in% c("AUG", "UGG")) %>% mutate(AA = codon_aa[Codon]) %>%
    group_by(AA) %>% mutate(n_syn = n(), RSCU_host = Freq * n_syn / sum(Freq)) %>%
    ungroup() %>% select(Codon, RSCU_host) %>% mutate(Host = host_name)
}
host_rscu <- purrr::map2_dfr(host_files, names(host_files), read_host_rscu)

cosine_sim <- function(a, b) { sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2))) }

sinv_wide <- sinv %>% select(Accession, Region, Codon, RSCU) %>% pivot_wider(names_from = Codon, values_from = RSCU)
host_wide <- host_rscu %>% select(Host, Codon, RSCU_host) %>% pivot_wider(names_from = Codon, values_from = RSCU_host)
codon_order <- sort(unique(sinv$Codon))

sinv_mat <- as.matrix(sinv_wide[, codon_order])
host_mat <- as.matrix(host_wide[, codon_order])
rownames(sinv_mat) <- paste(sinv_wide$Accession, sinv_wide$Region, sep = "__")
rownames(host_mat) <- host_wide$Host

R_mat <- outer(seq_len(nrow(sinv_mat)), seq_len(nrow(host_mat)), Vectorize(function(i, j) cosine_sim(sinv_mat[i, ], host_mat[j, ])))
rownames(R_mat) <- rownames(sinv_mat); colnames(R_mat) <- rownames(host_mat)

R_long <- as.data.frame(as.table(R_mat)) %>% rename(Vector = Var1, Host = Var2, R = Freq) %>%
  separate(Vector, into = c("Accession", "Region"), sep = "__") %>% mutate(SiD = (1 - R) / 2) %>% select(Accession, Region, Host, R, SiD)

dir.create(here("results/tables/SiD"), recursive = TRUE, showWarnings = FALSE)
write_csv(R_long, here("results/tables/SiD/SiD_long.csv"))
message("SiD_long.csv written: ", nrow(R_long), " rows")
