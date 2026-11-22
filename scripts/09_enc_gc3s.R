# ============================================================
# 09_enc_gc3s.R
# Purpose : ENc-GC3s plots for SINV CDS, NSP and SP, colored by host group;
#           compute expected ENc (ENCe), ENC ratio (ENCr) and perform
#           one-sample Wilcoxon signed-rank test of ENCr vs 0.
# Input   : results/tables/codonw/CodonUsageIndices_CDS.csv
#           results/tables/codonw/CodonUsageIndices_NSP.csv
#           results/tables/codonw/CodonUsageIndices_SP.csv
#           results/tables/host_class_summary.csv
# Output  : results/figures/ENc_GC3s_all.pdf
#           results/tables/ENc_GC3s/ENc_GC3s_per_sequence_CDS.csv
#           results/tables/ENc_GC3s/ENc_GC3s_per_sequence_NSP.csv
#           results/tables/ENc_GC3s/ENc_GC3s_per_sequence_SP.csv
#           results/tables/ENc_GC3s/ENc_GC3s_wilcox.csv
# Dependencies: here, tidyverse, patchwork
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(patchwork)
})

fam <- "ArialMT"
tab_dir <- here("results", "tables")
fig_dir <- here("results", "figures")
out_tab <- here("results", "tables", "ENc_GC3s")
dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)

summary_df <- read.csv(file.path(tab_dir, "host_class_summary.csv"), stringsAsFactors = FALSE)
map_list <- list()
for (i in 1:nrow(summary_df)) {
  accs <- unlist(strsplit(summary_df$All_accessions[i], "; "))
  map_list[[i]] <- data.frame(Accession = accs, Group = summary_df$Host_class[i])
}
meta <- bind_rows(map_list) %>% mutate(meta_row = row_number())

read_region <- function(region) {
  ind <- read.csv(file.path(tab_dir, "codonw", paste0("CodonUsageIndices_", region, ".csv")),
                  check.names = FALSE, stringsAsFactors = FALSE) %>%
    filter(!grepl("_Mean$|_SD$", Accession)) %>%
    select(Accession, ENCo = Nc, GC3s)
  ind %>%
    inner_join(meta, by = "Accession") %>%
    mutate(Region = region,
           Group = factor(Group, levels = c("Human", "Experimental_passage", "Other_vertebrate", "Bird", "Mosquito", "Other_arthropod", "Unknown"))) %>%
    mutate(ENCe = 2 + GC3s + 29 / (GC3s^2 + (1 - GC3s)^2),
           ENCr = (ENCe - ENCo) / ENCe)
}

regions <- c("CDS", "NSP", "SP")
dat_all <- bind_rows(lapply(regions, read_region))

n_labels <- paste0(summary_df$Host_class, " (n=", summary_df$N_sequences, ")")
names(n_labels) <- summary_df$Host_class

# Full plot ranges
x_full <- c(0.35, 0.65)
y_full <- c(45, 62)
curve_full <- data.frame(GC3s = seq(x_full[1], x_full[2], length.out = 1000)) %>%
  mutate(ENCe = 2 + GC3s + 29 / (GC3s^2 + (1 - GC3s)^2))

# Zoom plot ranges
x_zoom <- c(0.53, 0.59)
y_zoom <- c(50, 60)

grp_cols <- c(Human = "#E41A1C", Experimental_passage = "#681F72", Other_vertebrate = "#A65628",
              Bird = "#DDCE43", Mosquito = "#64A9E2", Other_arthropod = "#33A02C", Unknown = "#808080")

base_theme <- theme_classic(base_size = 10, base_family = fam) +
  theme(
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5, family = fam, margin = margin(b = 8)),
    axis.title.x = element_text(size = 16, face = "bold", family = fam, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, face = "bold", family = fam, margin = margin(r = 20)),
    axis.text    = element_text(size = 10, color = "black", family = fam),
    axis.line    = element_line(color = "black", linewidth = 0.75 * 0.352778),
    axis.ticks   = element_line(color = "black", linewidth = 0.75 * 0.352778),
    legend.title = element_text(size = 10, face = "bold", family = fam),
    legend.text  = element_text(size = 10, family = fam),
    plot.margin  = margin(8, 8, 8, 8),
    panel.grid.minor = element_blank()
  )
grid_theme <- theme(panel.grid.major = element_line(color = "grey85", linewidth = 0.3))

make_plot_full <- function(df, curve_df, x_rng, y_rng, title, show_legend = FALSE) {
  p <- ggplot() +
    geom_line(data = curve_df, aes(x = GC3s, y = ENCe), color = "black", linewidth = 0.8) +
    geom_point(data = df, aes(x = GC3s, y = ENCo, color = Group), shape = 16, size = 2.4, alpha = 0.9) +
    scale_color_manual(values = grp_cols, name = "Host group", labels = n_labels) +
    coord_cartesian(xlim = x_rng, ylim = y_rng) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(title = title, x = "GC3s", y = "Effective number of codons (ENc)") +
    base_theme + theme(panel.grid = element_blank())
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

make_plot_zoom <- function(df, x_rng, y_rng, title, show_legend = FALSE) {
  p <- ggplot() +
    geom_point(data = df, aes(x = GC3s, y = ENCo, color = Group), shape = 16, size = 2.4, alpha = 0.9) +
    scale_color_manual(values = grp_cols, name = "Host group", labels = n_labels) +
    coord_cartesian(xlim = x_rng, ylim = y_rng) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(title = title, x = "GC3s", y = "Effective number of codons (ENc)") +
    base_theme + grid_theme
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

plots <- list()
for (r in regions) {
  df_r <- dat_all %>% filter(Region == r)
  plots[[paste0(r, "_full")]] <- make_plot_full(df_r, curve_full, x_full, y_full, paste0(r, " (full range)"), show_legend = (r == "CDS"))
  plots[[paste0(r, "_zoom")]] <- make_plot_zoom(df_r, x_zoom, y_zoom, paste0(r, " (zoom)"), show_legend = FALSE)
}

combined <- (plots[["CDS_full"]] | plots[["NSP_full"]] | plots[["SP_full"]]) /
            (plots[["CDS_zoom"]] | plots[["NSP_zoom"]] | plots[["SP_zoom"]]) +
  plot_layout(widths = c(1, 1, 1), heights = c(1, 1), guides = "collect") & theme(legend.position = "right")

ggsave(file.path(fig_dir, "ENc_GC3s_all.pdf"), combined, width = 18, height = 10, device = "pdf")

# 1. Per-sequence tables
for (r in regions) {
  df_r <- dat_all %>% filter(Region == r)
  seq_out <- df_r %>%
    arrange(meta_row) %>%
    select(Accession, Group, GC3s, ENCo, ENCe, ENCr)
  write.csv(seq_out, file.path(out_tab, paste0("ENc_GC3s_per_sequence_", r, ".csv")), row.names = FALSE)
}

# 2. Wilcoxon summary table
wilcox_all <- bind_rows(lapply(regions, function(r) {
  df_r <- dat_all %>% filter(Region == r)
  wt <- wilcox.test(df_r$ENCr, mu = 0, alternative = "two.sided", exact = FALSE)
  data.frame(region = r, n = nrow(df_r),
             median_ENCr = median(df_r$ENCr),
             mean_ENCr = mean(df_r$ENCr),
             n_ENCr_gt_0.05 = sum(df_r$ENCr > 0.05),
             prop_ENCr_gt_0.05 = mean(df_r$ENCr > 0.05),
             statistic = unname(wt$statistic),
             p_value = wt$p.value)
}))

wilcox_all$p_value_fmt <- formatC(wilcox_all$p_value, format = "e", digits = 2)
wilcox_out <- wilcox_all %>% select(region, n, median_ENCr, mean_ENCr, n_ENCr_gt_0.05, prop_ENCr_gt_0.05, statistic, p_value = p_value_fmt)

write.csv(wilcox_out, file.path(out_tab, "ENc_GC3s_wilcox.csv"), row.names = FALSE)

cat("Done.\n")
cat("Figure:", file.path(fig_dir, "ENc_GC3s_all.pdf"), "\n")
cat("Tables written to:", out_tab, "\n")
