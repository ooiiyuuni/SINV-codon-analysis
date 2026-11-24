# ============================================================
# 16d_coa_permanova.R
# Purpose : PERMANOVA on SINV COA coordinates for CDS, NSP and SP,
#           testing whether host group or geographic region
#           explains codon usage variation. No sequences are
#           filtered: missing Host or Country is coded as
#           "Unknown" and kept as a group, so n_used = 315 for
#           both factors. A Note column flags host/geo collinearity.
# Input   : results/tables/COA/SINV/COA_coords_{CDS,NSP,SP}.csv
#           results/tables/metadata.csv
# Output  : results/tables/COA/SINV/COA_permanova_results.csv
# Dependencies: here, tidyverse, vegan
# Random seed : 12345 (permutations)
# Run     : Rscript scripts/16d_coa_permanova.R
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
stopifnot(!any(is.na(meta$HostGroup)))
stopifnot(!any(is.na(meta$GeoRegion)))

# sanity: each group must have >= 2 samples for adonis2
print(table(meta$HostGroup))
print(table(meta$GeoRegion))

# ---------- run PERMANOVA (no filtering) ----------
results <- list()

for (r in c("CDS", "NSP", "SP")) {
  d <- read_csv(file.path(coa_dir, paste0("COA_coords_", r, ".csv")),
                show_col_types = FALSE) %>%
    inner_join(meta, by = "Accession")
  stopifnot(nrow(d) == n_expected)
  stopifnot(!any(is.na(d$Axis1)))
  stopifnot(!any(is.na(d$Axis2)))

  xy <- d %>% select(Axis1, Axis2) %>% as.data.frame()

  ad_host <- adonis2(xy ~ HostGroup, data = d, method = "euclidean", permutations = 999)
  results[[length(results) + 1]] <- tibble(
    Region = r, Factor = "HostGroup", n_used = nrow(d),
    df = ad_host$Df[1], SumSq = ad_host$SumOfSqs[1], R2 = ad_host$R2[1],
    F = ad_host$F[1], P = ad_host$`Pr(>F)`[1],
    Note = "All 315 sequences; Unknown kept as a group; host/geo collinear on USA/Rodent"
  )

  ad_geo <- adonis2(xy ~ GeoRegion, data = d, method = "euclidean", permutations = 999)
  results[[length(results) + 1]] <- tibble(
    Region = r, Factor = "GeoRegion", n_used = nrow(d),
    df = ad_geo$Df[1], SumSq = ad_geo$SumOfSqs[1], R2 = ad_geo$R2[1],
    F = ad_geo$F[1], P = ad_geo$`Pr(>F)`[1],
    Note = "All 315 sequences; Unknown (n=2) kept as a group; host/geo collinear on USA/Rodent"
  )
}

res <- bind_rows(results) %>%
  mutate(R2 = round(R2, 3), F = round(F, 3), P = signif(P, 3))

write_csv(res, file.path(coa_dir, "COA_permanova_results.csv"))

cat("Done.\n")
cat("PERMANOVA results:\n")
print(res, n = 10)
cat("\nFile written to:", file.path(coa_dir, "COA_permanova_results.csv"), "\n")
