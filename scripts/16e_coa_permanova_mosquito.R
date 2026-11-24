# ============================================================
# 16e_coa_permanova_mosquito.R
# Purpose : Supplementary PERMANOVA on SINV COA coordinates using
#           mosquito-only sequences (n = 136), to test whether the
#           geographic effect is independent of the 143 rodent
#           sequences that are collinear with Americas. GeoRegion
#           Unknown is excluded (mosquito-only already excludes
#           rodent; Unknown has no country). No other filtering.
# Input   : results/tables/COA/SINV/COA_coords_{CDS,NSP,SP}.csv
#           results/tables/metadata.csv
# Output  : results/tables/COA/SINV/COA_permanova_mosquito_only.csv
# Dependencies: here, tidyverse, vegan
# Random seed : 12345 (permutations)
# Run     : Rscript scripts/16e_coa_permanova_mosquito.R
# ============================================================

suppressPackageStartupMessages({
  library(here); library(tidyverse); library(vegan)
})
set.seed(12345)

n_expected <- 315L
coa_dir <- here("results", "tables", "COA", "SINV")

# ---------- metadata ----------
meta <- read_csv(here("results", "tables", "metadata.csv"), show_col_types = FALSE) %>%
  select(Accession, Host, Country)

host_group <- function(x) {
  x_low <- tolower(trimws(x))
  ifelse(is.na(x) | x_low == "", "Unknown",
    case_when(
      str_detect(x_low, "^(culex|aedes|ochlerotatus|anopheles|culiseta|mansonia|coquillettidia|culicidae|culicoides|mimomyia)") |
        x_low %in% c("mosquito", "midge") |
        str_detect(x_low, "^phlebotomine") ~ "Mosquito",
      x_low %in% c("cricetulus griseus", "mesocricetus auratus", "cricetus cricetus") ~ "Rodent",
      x_low == "homo sapiens" ~ "Human",
      str_detect(x_low, "^(streptopelia|somateria|nycticorax|corvus|columba|bubulcus)") ~ "Bird",
      TRUE ~ "Other"
    ))
}

region_group <- function(x) {
  prefix <- str_extract(x, "^[^:]+") %>% str_trim()
  ifelse(is.na(x) | trimws(x) == "", "Unknown",
    case_when(
      prefix %in% c("USA", "Brazil") ~ "Americas",
      prefix %in% c("Sweden", "Spain", "Finland", "Germany", "Italy",
                    "Slovenia", "Russia", "Slovakia", "Norway",
                    "Armenia", "Azerbaijan") ~ "Europe",
      prefix %in% c("Australia", "Papua New Guinea") ~ "Oceania",
      prefix %in% c("South Africa", "Central African Republic", "Kenya",
                    "Uganda", "Algeria", "Cameroon", "Egypt",
                    "Senegal", "Zambia") ~ "Africa",
      prefix %in% c("Saudi Arabia", "Israel", "China", "Malaysia") ~ "Asia",
      TRUE ~ "Unknown"
    ))
}

meta <- meta %>%
  mutate(
    HostGroup = factor(host_group(Host),
                       levels = c("Mosquito", "Rodent", "Human", "Bird", "Other", "Unknown")),
    GeoRegion = factor(region_group(Country),
                       levels = c("Americas", "Europe", "Oceania", "Africa", "Asia", "Unknown"))
  )
stopifnot(nrow(meta) == n_expected)

# ---------- mosquito-only PERMANOVA ----------
results <- list()

for (r in c("CDS", "NSP", "SP")) {
  d <- read_csv(file.path(coa_dir, paste0("COA_coords_", r, ".csv")),
                show_col_types = FALSE) %>%
    inner_join(meta, by = "Accession") %>%
    filter(HostGroup == "Mosquito", GeoRegion != "Unknown")

  cat("\n", r, " mosquito-only n =", nrow(d), "\n")
  print(table(d$GeoRegion))

  if (nrow(d) >= 10 && length(unique(d$GeoRegion)) >= 2) {
    xy <- d %>% select(Axis1, Axis2) %>% as.data.frame()
    ad <- adonis2(xy ~ GeoRegion, data = d, method = "euclidean", permutations = 999)
    results[[length(results) + 1]] <- tibble(
      Region = r, Factor = "GeoRegion (mosquito only)", n_used = nrow(d),
      df = ad$Df[1], SumSq = ad$SumOfSqs[1], R2 = ad$R2[1],
      F = ad$F[1], P = ad$`Pr(>F)`[1],
      Note = "Mosquito-only subset; Unknown excluded; independent geographic test"
    )
  } else {
    cat("  -> skipped: too few samples or groups\n")
  }
}

res <- bind_rows(results) %>%
  mutate(R2 = round(R2, 3), F = round(F, 3), P = signif(P, 3))

write_csv(res, file.path(coa_dir, "COA_permanova_mosquito_only.csv"))

cat("\nDone.\n")
cat("Mosquito-only GeoRegion PERMANOVA:\n")
print(res, n = 10)
cat("\nFile written to:", file.path(coa_dir, "COA_permanova_mosquito_only.csv"), "\n")
