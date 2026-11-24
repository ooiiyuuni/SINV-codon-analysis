# =============================================================================
# 07b_sid_figure.R
# Purpose : Summarise mean SiD by host and genomic region, write the summary
#           table, and generate a grouped bar chart of mean SiD between SINV
#           sequences and five host reference profiles.
# Input   : results/tables/SiD/SiD_long.csv
# Output  : results/tables/SiD/SiD_summary.csv
#           results/figures/SiD_barplot.pdf
# Dependencies: here, tidyverse
# Run     : Rscript scripts/07b_sid_figure.R
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

set.seed(12345)

# ---------------------------------------------------------------------------
# 1. Read per-sequence SiD
# ---------------------------------------------------------------------------
sid <- read_csv(here("results/tables/SiD/SiD_long.csv"), show_col_types = FALSE)

# ---------------------------------------------------------------------------
# 2. Fixed host order and region order
# ---------------------------------------------------------------------------
host_order <- c("Culex pipiens", "Turdidae", "Passeriformes", "Equus caballus", "Homo sapiens")
region_order <- c("CDS", "NSP", "SP")

# ---------------------------------------------------------------------------
# 3. Summarise mean SiD by Host and Region, write summary table
# ---------------------------------------------------------------------------
sid_summary <- sid %>%
  mutate(
    Host   = factor(Host,   levels = host_order),
    Region = factor(Region, levels = region_order)
  ) %>%
  group_by(Host, Region) %>%
  summarise(mean_SiD = mean(SiD), .groups = "drop")

dir.create(here("results/tables/SiD"), recursive = TRUE, showWarnings = FALSE)
write_csv(sid_summary, here("results/tables/SiD/SiD_summary.csv"))

# ---------------------------------------------------------------------------
# 4. Plot
# ---------------------------------------------------------------------------
PT <- 0.352778

df_long <- sid_summary %>%
  mutate(Label = sprintf("%.4f", mean_SiD))

fill_colors <- c("CDS" = "#53629B", "NSP" = "#FFDCE4", "SP" = "#C46C81")

host_labels <- c(
  "Culex pipiens"     = expression(italic("Culex pipiens")),
  "Turdidae"          = expression(italic("Turdidae")),
  "Passeriformes"     = expression(italic("Passeriformes")),
  "Equus caballus"    = expression(italic("Equus caballus")),
  "Homo sapiens"      = expression(italic("Homo sapiens"))
)

p <- ggplot(df_long, aes(x = Host, y = mean_SiD, fill = Region)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  geom_text(aes(label = Label, y = mean_SiD + 0.0015), 
            position = position_dodge(width = 0.8), 
            size = 14/2.845, family = "Arial", vjust = 0) +
  scale_fill_manual(values = fill_colors, name = "Region") +
  scale_x_discrete(labels = host_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = NULL, 
    y = "Similarity Index (SiD)",
    title = "Similarity Index (SiD) of SINV to Hosts"
  ) +
  theme_classic(base_size = 18, base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", size = 18, color = "black"),
    axis.text.x = element_text(family = "Arial", size = 18, color = "black"),
    axis.text.y = element_text(family = "Arial", size = 18, color = "black"),
    axis.title.y = element_text(family = "Arial", size = 18, color = "black"),
    axis.line = element_line(color = "black", linewidth = 1 * PT),
    axis.ticks = element_line(color = "black", linewidth = 1 * PT),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(family = "Arial", size = 24, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

# ---------------------------------------------------------------------------
# 5. Save
# ---------------------------------------------------------------------------
out_dir <- here("results","figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(here(out_dir, "SiD_barplot.pdf"), p, width = 22, height = 8, device = cairo_pdf)

message("SiD_summary.csv and SiD_barplot.pdf written")
