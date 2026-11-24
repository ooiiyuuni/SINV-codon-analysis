# ============================================================
# 16b_coa_biplot.R
# Purpose : COA biplot for SINV CDS, NSP and SP, colored by host
#           group and by geographic region. Host groups and colors
#           are strictly aligned with 09_enc_gc3s.R.
# Input   : results/tables/COA/SINV/COA_coords_{CDS,NSP,SP}.csv
#           results/tables/COA/SINV/COA_axis_contribution.csv
#           results/tables/codonw/CodonUsageIndices_{CDS,NSP,SP}.csv
#           results/tables/BaseComposition_{CDS,NSP,SP}.csv
#           results/tables/host_class_summary.csv
# Output  : results/figures/COA_biplot_by_host_SINV.pdf
#           results/figures/COA_biplot_by_region_SINV.pdf
#           results/tables/COA/SINV/COA_arrow_correlations.csv
#           results/tables/COA/SINV/COA_axis_contribution_wide.csv
# Dependencies: here, tidyverse, patchwork, ggrepel
# Run     : Rscript scripts/16b_coa_biplot.R
# ============================================================

suppressPackageStartupMessages({
  library(here); library(tidyverse); library(patchwork); library(ggrepel)
})

n_expected <- 315L
coa_dir <- here("results", "tables", "COA", "SINV")
out_fig <- here("results", "figures")
tab_dir <- here("results", "tables")
if (!dir.exists(out_fig)) dir.create(out_fig, recursive = TRUE)

# ---------- 1. COA coordinates ----------
coords <- map_dfr(c("CDS", "NSP", "SP"), function(r) {
  read_csv(file.path(coa_dir, paste0("COA_coords_", r, ".csv")), show_col_types = FALSE) %>%
    mutate(Region = r)
}) %>% mutate(Region = factor(Region, levels = c("CDS", "NSP", "SP")))
stopifnot(nrow(coords) == n_expected * 3)

# ---------- 2. metadata and grouping (aligned with 09) ----------
# Read host mapping from 09's output
summary_df <- read.csv(file.path(tab_dir, "host_class_summary.csv"), stringsAsFactors = FALSE)
map_list <- list()
for (i in 1:nrow(summary_df)) {
  accs <- unlist(strsplit(summary_df$All_accessions[i], "; "))
  map_list[[i]] <- data.frame(Accession = accs, Group = summary_df$Host_class[i])
}
meta <- bind_rows(map_list) %>% mutate(meta_row = row_number())

# Read original metadata for Country
orig_meta <- read_csv(here("results", "tables", "metadata.csv"), show_col_types = FALSE) %>%
  select(Accession, Country)

meta <- meta %>% left_join(orig_meta, by = "Accession")

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

meta <- meta %>%
  mutate(
    HostGroup = factor(Group, levels = c("Human", "Experimental_passage", "Other_vertebrate", "Bird", "Mosquito", "Other_arthropod", "Unknown")),
    GeoRegion = factor(region_group(Country), levels = c("Americas", "Europe", "Oceania", "Africa", "Asia", "Unknown"))
  )

stopifnot(nrow(meta) == n_expected)
stopifnot(!any(is.na(meta$HostGroup)))
stopifnot(!any(is.na(meta$GeoRegion)))

# ---------- 3. codon usage indices + base composition ----------
read_vars <- function(region) {
  cod <- read_csv(here("results", "tables", "codonw",
                       paste0("CodonUsageIndices_", region, ".csv")),
                  show_col_types = FALSE) %>%
    filter(!str_detect(Accession, "_Mean$|_SD$")) %>%
    select(Accession, A3s, U3s, G3s, C3s, GC3s, Nc, Gravy, Aromo)

  bc <- read_csv(here("results", "tables",
                      paste0("BaseComposition_", region, ".csv")),
                 show_col_types = FALSE) %>%
    filter(!str_detect(Accession, "_Mean$|_SD$")) %>%
    select(Accession, `A%`, `U%`, `G%`, `C%`, `GC%`) %>%
    rename(A = `A%`, U = `U%`, G = `G%`, C = `C%`, GC = `GC%`)

  inner_join(cod, bc, by = "Accession")
}
vars_by_region <- lapply(c(CDS = "CDS", NSP = "NSP", SP = "SP"), read_vars)
for (nm in names(vars_by_region)) {
  stopifnot(nrow(vars_by_region[[nm]]) == n_expected)
  stopifnot(!any(is.na(vars_by_region[[nm]])))
}

dat <- map_dfr(c("CDS", "NSP", "SP"), function(r) {
  coords %>% filter(Region == r) %>%
    inner_join(vars_by_region[[r]], by = "Accession")
}) %>%
  left_join(meta, by = "Accession") %>%
  mutate(Region = factor(Region, levels = c("CDS", "NSP", "SP")))

stopifnot(nrow(dat) == n_expected * 3)
stopifnot(!any(is.na(dat$HostGroup)))
stopifnot(!any(is.na(dat$GeoRegion)))

# ---------- 4. Spearman arrow correlations ----------
var_names <- c("GC3s", "Nc", "Gravy", "Aromo",
               "A", "U", "G", "C", "GC",
               "A3s", "U3s", "G3s", "C3s")
arrow_tab <- map_dfr(c("CDS", "NSP", "SP"), function(r) {
  d <- dat %>% filter(Region == r)
  map_dfr(var_names, function(v) {
    x <- d[[v]]
    ct1 <- suppressWarnings(cor.test(x, d$Axis1, method = "spearman"))
    ct2 <- suppressWarnings(cor.test(x, d$Axis2, method = "spearman"))
    tibble(
      Region = r, Variable = v,
      r_Axis1 = unname(ct1$estimate), P_Axis1 = ct1$p.value,
      r_Axis2 = unname(ct2$estimate), P_Axis2 = ct2$p.value
    )
  })
})
write_csv(arrow_tab, file.path(coa_dir, "COA_arrow_correlations.csv"))

# ---------- 5. axis labels and wide contribution ----------
contrib <- read_csv(file.path(coa_dir, "COA_axis_contribution.csv"), show_col_types = FALSE)
axis_lab <- function(region) {
  v <- contrib %>% filter(Region == !!region, Axis %in% c(1, 2)) %>% arrange(Axis)
  paste0("Axis ", v$Axis, " (", sprintf("%.1f", v$Contribution_percent), "%)")
}
axis_wide <- contrib %>%
  filter(Axis <= 4) %>%
  select(Region, Axis, Contribution_percent, Cumulative_percent) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  pivot_wider(names_from = Region,
              values_from = c(Contribution_percent, Cumulative_percent),
              names_glue = "{Region}_{.value}") %>%
  select(Axis, starts_with("CDS_"), starts_with("NSP_"), starts_with("SP_"))
write_csv(axis_wide, file.path(coa_dir, "COA_axis_contribution_wide.csv"))

# ---------- 6. plot limits + isotropic arrow scale ----------
r_x <- max(abs(dat$Axis1)) * 1.08
r_y <- max(abs(dat$Axis2)) * 1.08
x_lim <- c(-r_x, r_x); y_lim <- c(-r_y, r_y)
br_x <- pretty(x_lim, n = 5); br_y <- pretty(y_lim, n = 5)
scale_factor <- min(r_x, r_y) * 0.8

# ---------- 7. disjoint palettes + n in legend labels (aligned with 09) ----------
host_cols <- c(Human = "#E41A1C", Experimental_passage = "#681F72", Other_vertebrate = "#A65628",
               Bird = "#DDCE43", Mosquito = "#64A9E2", Other_arthropod = "#33A02C", Unknown = "#808080")
geo_cols  <- c(Americas = "#1B9E77", Europe = "#D95F02", Oceania = "#7570B3",
               Africa = "#E7298A", Asia = "#66A61E", Unknown = "#BDBDBD")

host_n <- dat %>% distinct(Accession, HostGroup) %>% count(HostGroup, .drop = FALSE)
geo_n  <- dat %>% distinct(Accession, GeoRegion)  %>% count(GeoRegion,  .drop = FALSE)
host_labels <- setNames(paste0(host_n$HostGroup, " (n = ", host_n$n, ")"), host_n$HostGroup)
geo_labels  <- setNames(paste0(geo_n$GeoRegion,  " (n = ", geo_n$n,  ")"), geo_n$GeoRegion)

# ---------- 8. theme ----------
common_theme <- theme_classic(base_size = 10, base_family = "ArialMT") +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0, family = "ArialMT"),
    axis.text = element_text(size = 18, color = "black", family = "ArialMT"),
    axis.title = element_text(size = 18, family = "ArialMT"),
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.title.x = element_text(margin = margin(t = 12)),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, family = "ArialMT"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# ---------- 9. plot factory ----------
make_plot <- function(region, colour_by) {
  d <- dat %>% filter(Region == region)

  a <- arrow_tab %>%
    filter(Region == !!region) %>%
    mutate(r_norm = sqrt(r_Axis1^2 + r_Axis2^2)) %>%
    filter(r_norm >= 0.5) %>%
    mutate(xend = r_Axis1 * scale_factor,
           yend = r_Axis2 * scale_factor)

  xl <- axis_lab(region)[1]; yl <- axis_lab(region)[2]
  n_host <- d %>% distinct(Accession) %>% nrow()

  if (colour_by == "host") {
    p <- ggplot(d, aes(Axis1, Axis2, colour = HostGroup)) +
      geom_point(size = 2.2, alpha = 0.9) +
      scale_colour_manual(values = host_cols, labels = host_labels, drop = FALSE)
  } else {
    p <- ggplot(d, aes(Axis1, Axis2, colour = GeoRegion)) +
      geom_point(size = 2.2, alpha = 0.9) +
      scale_colour_manual(values = geo_cols, labels = geo_labels, drop = FALSE)
  }

  p +
    annotate("segment", x = -r_x, y = 0, xend = r_x, yend = 0,
             arrow = arrow(length = unit(0.02, "npc")),
             colour = "black", linewidth = 0.75 * 0.352778) +
    annotate("segment", x = 0, y = -r_y, xend = 0, yend = r_y,
             arrow = arrow(length = unit(0.02, "npc")),
             colour = "black", linewidth = 0.75 * 0.352778) +
    geom_segment(data = a,
                 aes(x = 0, y = 0, xend = xend, yend = yend),
                 inherit.aes = FALSE,
                 arrow = arrow(length = unit(0.02, "npc")),
                 colour = "grey20", linewidth = 0.4) +
    geom_text_repel(data = a,
                    aes(x = xend, y = yend, label = Variable),
                    inherit.aes = FALSE,
                    size = 10 * 0.352778, family = "ArialMT", colour = "grey20",
                    box.padding = 0.35, point.padding = 0.2,
                    min.segment.length = 0,
                    segment.colour = "grey60", segment.size = 0.3,
                    max.overlaps = Inf) +
    scale_x_continuous(limits = x_lim, breaks = br_x, expand = c(0, 0)) +
    scale_y_continuous(limits = y_lim, breaks = br_y, expand = c(0, 0)) +
    labs(title = paste0(region, " - ",
                        ifelse(colour_by == "host", "by host group", "by geographic region"),
                        " (n = ", n_host, ")"),
         x = xl, y = yl) +
    common_theme
}

# ---------- 10. combine and export ----------
p_host <- make_plot("CDS", "host") + make_plot("NSP", "host") + make_plot("SP", "host") +
  plot_layout(guides = "collect") & theme(legend.position = "right")
ggsave(file.path(out_fig, "COA_biplot_by_host_SINV.pdf"), p_host,
       width = 24, height = 7, device = "pdf")

p_geo <- make_plot("CDS", "region") + make_plot("NSP", "region") + make_plot("SP", "region") +
  plot_layout(guides = "collect") & theme(legend.position = "right")
ggsave(file.path(out_fig, "COA_biplot_by_region_SINV.pdf"), p_geo,
       width = 24, height = 7, device = "pdf")

cat("Done.\n")
cat("Figures written to:", out_fig, "\n")
cat("Arrow correlations:", file.path(coa_dir, "COA_arrow_correlations.csv"), "\n")
