# ============================================================
# 16c_coa_clustering.R
# Purpose : k-means (k = 3-6) + silhouette widths on SINV COA
#           space for CDS, NSP and SP. k = 4 is chosen for
#           interpretation (k = 5 only splits an n = 6 noise
#           cluster). Host-group and geographic composition are
#           tabulated at k = 4, plus cross-region membership.
# Input   : results/tables/COA/SINV/COA_coords_{CDS,NSP,SP}.csv
#           results/tables/metadata.csv
# Output  : results/tables/COA/SINV/COA_cluster_silhouette.csv
#           results/tables/COA/SINV/COA_cluster_membership.csv
#           results/tables/COA/SINV/COA_cluster_host_geo_counts.csv
#           results/tables/COA/SINV/COA_cluster_cross_region.csv
# Dependencies: here, tidyverse, cluster
# Random seed : 12345 (k-means starts)
# Run     : Rscript scripts/16c_coa_clustering.R
# ============================================================

suppressPackageStartupMessages({
  library(here); library(tidyverse); library(cluster)
})
set.seed(12345)

n_expected <- 315L
coa_dir <- here("results", "tables", "COA", "SINV")

region_group <- function(x) {
  prefix <- str_extract(x, "^[^:]+") %>% str_trim()
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
  )
}

# ---------- metadata with host / geo groups ----------
meta <- read_csv(here("results", "tables", "metadata.csv"), show_col_types = FALSE) %>%
  select(Accession, Host, Country)

# host classification from 09_enc_gc3s.R (host_class_summary.csv)
summary_df <- read.csv(here("results", "tables", "host_class_summary.csv"),
                       stringsAsFactors = FALSE)
map_list <- list()
for (i in 1:nrow(summary_df)) {
  accs <- unlist(strsplit(summary_df$All_accessions[i], "; "))
  map_list[[i]] <- data.frame(Accession = accs, HostGroup = summary_df$Host_class[i])
}
host_map <- bind_rows(map_list)

meta <- meta %>%
  left_join(host_map, by = "Accession") %>%
  mutate(
    HostGroup = factor(HostGroup,
                       levels = c("Human", "Experimental_passage", "Other_vertebrate",
                                  "Bird", "Mosquito", "Other_arthropod", "Unknown")),
    GeoRegion = factor(region_group(Country),
                       levels = c("Americas", "Europe", "Oceania", "Africa", "Asia", "Unknown"))
  )
stopifnot(nrow(meta) == n_expected)
stopifnot(!any(is.na(meta$HostGroup)))
stopifnot(!any(is.na(meta$GeoRegion)))

# ---------- read coordinates ----------
coords <- lapply(c(CDS = "CDS", NSP = "NSP", SP = "SP"), function(r) {
  read_csv(file.path(coa_dir, paste0("COA_coords_", r, ".csv")), show_col_types = FALSE)
})
for (nm in names(coords)) {
  stopifnot(nrow(coords[[nm]]) == n_expected)
  stopifnot(!any(is.na(coords[[nm]]$Axis1)))
  stopifnot(!any(is.na(coords[[nm]]$Axis2)))
}

# ---------- k-means and silhouette ----------
ks <- 3:6
sil_tab <- list(); membership <- list()

for (r in names(coords)) {
  d <- coords[[r]]
  xy <- as.matrix(d[, c("Axis1", "Axis2")])
  rownames(xy) <- d$Accession

  for (k in ks) {
    km <- kmeans(xy, centers = k, nstart = 50)
    sil <- silhouette(km$cluster, dist(xy))
    mean_sil <- mean(sil[, "sil_width"])

    sil_tab[[length(sil_tab) + 1]] <- tibble(
      Region = r, k = k, Mean_silhouette = mean_sil
    )
    membership[[length(membership) + 1]] <- tibble(
      Region = r, k = k, Accession = rownames(xy),
      Cluster = km$cluster, Sil_width = sil[, "sil_width"]
    )
  }
}

sil_tab <- bind_rows(sil_tab)
membership <- bind_rows(membership)

# ---------- silhouette table ----------
sil_out <- sil_tab %>%
  mutate(Mean_silhouette = round(Mean_silhouette, 3)) %>%
  pivot_wider(names_from = Region, values_from = Mean_silhouette) %>%
  select(k, CDS, NSP, SP)
write_csv(sil_out, file.path(coa_dir, "COA_cluster_silhouette.csv"))

# ---------- attach metadata ----------
membership <- membership %>% left_join(meta, by = "Accession")
stopifnot(nrow(membership) == n_expected * 3 * length(ks))
stopifnot(!any(is.na(membership$HostGroup)))
write_csv(membership, file.path(coa_dir, "COA_cluster_membership.csv"))

# ---------- host / geo composition per cluster (k = 4) ----------
counts <- membership %>%
  filter(k == 4) %>%
  group_by(Region, k, Cluster) %>%
  summarise(
    n = n(),
    Human                = sum(HostGroup == "Human"),
    Experimental_passage = sum(HostGroup == "Experimental_passage"),
    Other_vertebrate     = sum(HostGroup == "Other_vertebrate"),
    Bird                 = sum(HostGroup == "Bird"),
    Mosquito             = sum(HostGroup == "Mosquito"),
    Other_arthropod      = sum(HostGroup == "Other_arthropod"),
    Unknown_host         = sum(HostGroup == "Unknown"),
    Americas = sum(GeoRegion == "Americas"),
    Europe   = sum(GeoRegion == "Europe"),
    Oceania  = sum(GeoRegion == "Oceania"),
    Africa   = sum(GeoRegion == "Africa"),
    Asia     = sum(GeoRegion == "Asia"),
    Unknown_geo = sum(GeoRegion == "Unknown"),
    .groups = "drop"
  )
write_csv(counts, file.path(coa_dir, "COA_cluster_host_geo_counts.csv"))

# ---------- cross-region cluster membership (k = 4) ----------
cross <- membership %>%
  filter(k == 4) %>%
  select(Region, Accession, Cluster) %>%
  pivot_wider(names_from = Region, values_from = Cluster, names_prefix = "Cluster_")
write_csv(cross, file.path(coa_dir, "COA_cluster_cross_region.csv"))

# ---------- console report ----------
cat("Done.\n")
cat("Silhouette widths (mean, k = 3-6):\n")
print(sil_out)
cat("\nCluster composition at k = 4:\n")
print(counts)
cat("\nCross-region membership (head):\n")
print(head(cross, 10))
cat("\nFiles written to:", coa_dir, "\n")
