# Purpose: Generate RCDI table and barplot for SINV.
# Input: Data from RCDI/eRCDI web server (Puigbo et al., 2010).
# Output: results/tables/RCDI/RCDI_web_summary.csv
#         results/figures/RCDI_barplot.pdf
# Dependencies: here, tidyverse
# Random seed: no stochastic step

library(here)
library(tidyverse)

out_tab <- here("results", "tables", "RCDI")
out_fig <- here("results", "figures")
dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

host_levels <- c("Culex_pipiens", "Turdidae", "Passeriformes", 
                 "Equus_caballus", "Homo_sapiens")
region_levels <- c("NSP", "SP")

dat <- data.frame(
  Region = factor(rep(region_levels, each = 5), levels = region_levels),
  Host = factor(rep(host_levels, 2), levels = host_levels),
  Average_RCDI = c(
    1.404, 1.322, 1.302, 1.290, 1.323,
    1.352, 1.326, 1.310, 1.298, 1.347
  ),
  eRCDI_threshold = c(
    1.555, 1.454, 1.430, 1.410, 1.452,
    1.493, 1.467, 1.444, 1.422, 1.505
  )
)

readr::write_csv(dat, file.path(out_tab, "RCDI_web_summary.csv"))

# 修改了颜色
region_colors <- c(NSP = "#FFDCE4", SP = "#C46C81")

common_theme <- theme_classic(base_size = 18, base_family = "ArialMT") +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0, margin = margin(b = 10), family = "ArialMT"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 18, color = "black", family = "ArialMT"),
    axis.text.y = element_text(size = 18, color = "black", family = "ArialMT"),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 12), family = "ArialMT"),
    legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 18, family = "ArialMT"),
    axis.line = element_line(linewidth = 0.75 * 0.352778, colour = "black"),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3), panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.text = element_text(size = 18, face = "bold", family = "ArialMT")
  )

p_rcdi <- ggplot(dat, aes(x = Host, y = Average_RCDI)) +
  geom_bar(aes(fill = Region), stat = "identity", width = 0.8, colour = NA) +
  geom_crossbar(aes(y = eRCDI_threshold, ymin = eRCDI_threshold, ymax = eRCDI_threshold), 
                width = 0.8, colour = "#E41A1C", linetype = "dashed", linewidth = 0.5, fatten = 0) +
  geom_text(aes(label = sprintf("%.3f", Average_RCDI), y = Average_RCDI + 0.02), 
            size = 18 * 0.352778, family = "ArialMT", vjust = 0, hjust = 0.5) +
  geom_text(aes(label = sprintf("%.3f", eRCDI_threshold), y = eRCDI_threshold + 0.02), 
            size = 18 * 0.352778, family = "ArialMT", vjust = 0, hjust = 0.5, colour = "#E41A1C") +
  facet_wrap(~ Region, nrow = 1) +
  labs(title = "Relative Codon Deoptimization Index (RCDI) across hosts", x = NULL, y = "RCDI") +
  scale_fill_manual(values = region_colors) +
  scale_x_discrete(labels = function(x) parse(text = paste0("italic('", gsub("_", " ", x), "')"))) +
  scale_y_continuous(limits = c(0, 2.0), breaks = seq(0, 2.0, 0.5), expand = c(0, 0)) +
  common_theme

ggsave(file.path(out_fig, "RCDI_barplot.pdf"), p_rcdi, width = 18, height = 7, device = "pdf")
