# ============================================================
# 06d_rscu_heatmap.R
# Purpose: Generate a single RSCU heatmap for SINV and five hosts.
#          Groups codons by amino acid with bracket lines.
#          Enlarged fonts, centered text, added title.
# Input:
#   results/tables/RSCU/RSCU.csv
# Output:
#   results/figures/RSCU_heatmap.pdf
# ============================================================

set.seed(12345)

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr)
  library(readr); library(ggplot2)
})

PT <- 0.352778

# ---------- 1. Read RSCU table ----------
rscu <- read_csv(here("results","tables","RSCU","RSCU.csv"),
                 show_col_types = FALSE,
                 col_types = cols(.default = col_character()))
num_cols <- setdiff(names(rscu), c("AA","Codon"))
rscu[num_cols] <- lapply(rscu[num_cols], as.numeric)

rscu59 <- rscu %>% filter(!Codon %in% c("AUG","UGG"))

# ---------- 2. Datasets ----------
data_cols <- c("Complete","NSP","SP",
               "Culex pipiens","Turdidae","Passeriformes",
               "Equus caballus","Homo sapiens")
dataset_levels <- rev(data_cols)

# ---------- 3. Setup codon order & index ----------
rscu59 <- rscu59 %>%
  mutate(Codon = factor(Codon, levels = Codon)) %>%
  arrange(Codon) %>%
  mutate(CodonIdx = row_number())

# ---------- 4. Long format ----------
long <- rscu59 %>%
  select(Codon, CodonIdx, all_of(data_cols)) %>%
  pivot_longer(cols = all_of(data_cols), names_to = "Dataset", values_to = "RSCU")
long$Dataset <- factor(long$Dataset, levels = dataset_levels)

# ---------- 5. Discretized colour scale ----------
long <- long %>%
  mutate(RSCU_cut = cut(RSCU,
                        breaks = c(-Inf, 0.6, 0.8, 1.0, 1.2, 1.6, Inf),
                        labels = c("<0.6","0.6-0.8","0.8-1.0","1.0-1.2","1.2-1.6",">1.6"),
                        right = FALSE))
fill_colors <- c("<0.6"="#053061","0.6-0.8"="#92C5DE","0.8-1.0"="#F7F7F7",
                 "1.0-1.2"="#FDDBC7","1.2-1.6"="#D6604D",">1.6"="#67001F")

# ---------- 6. Bracket bars for amino acids ----------
gap <- 0.2
aa_bars <- rscu59 %>% 
  group_by(AA) %>%
  summarise(start = min(CodonIdx) - 0.5 + gap,
            end   = max(CodonIdx) + 0.5 - gap,
            mid   = mean(CodonIdx), .groups = "drop")

# ---------- 7. Y-axis labels ----------
dataset_labels <- c(
  "Complete"          = "CDS",
  "NSP"               = "NSP",
  "SP"                = "SP",
  "Culex pipiens"     = expression(italic("Culex pipiens")),
  "Turdidae"          = expression(italic("Turdidae")),
  "Passeriformes"     = expression(italic("Passeriformes")),
  "Equus caballus"    = expression(italic("Equus caballus")),
  "Homo sapiens"      = expression(italic("Homo sapiens"))
)

# ---------- 8. Plot ----------
n_ds <- length(data_cols)
y_codon   <- 0
y_aa_bar  <- -1.0
y_aa_text <- -1.8

p <- ggplot(long, aes(x = CodonIdx, y = Dataset, fill = RSCU_cut)) +
  geom_tile(color = NA) +
  scale_fill_manual(values = fill_colors, name = "RSCU", drop = FALSE) +
  
  geom_text(data = rscu59, aes(x = CodonIdx, y = y_codon, label = Codon),
            inherit.aes = FALSE, size = 18/2.845, angle = 90, family = "Arial",
            color = "black", vjust = 0.5, hjust = 0.5) +
            
  geom_segment(data = aa_bars,
               aes(x = start, xend = end, y = y_aa_bar, yend = y_aa_bar),
               inherit.aes = FALSE, linewidth = 1.5 * PT, color = "black") +
               
  geom_text(data = aa_bars,
            aes(x = mid, y = y_aa_text, label = AA),
            inherit.aes = FALSE, size = 18/2.845,
            fontface = "plain", family = "Arial", color = "black") +

  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  scale_y_discrete(expand = c(0, 0), labels = dataset_labels) +
  coord_cartesian(ylim = c(y_aa_text - 1.5, n_ds + 0.5), clip = "off") +
  
  labs(x = NULL, y = NULL, 
       title = "Relative Synonymous Codon Usage (RSCU) of SINV and Hosts",
       caption = "RSCU > 1 indicates positive bias; RSCU < 1 indicates negative bias") +
       
  theme_bw(base_size = 18, base_family = "Arial") +
  theme(
    text          = element_text(family = "Arial", size = 18, color = "black", hjust = 0.5),
    axis.text.y   = element_text(family = "Arial", size = 24, color = "black"),
    axis.text.x   = element_blank(),
    plot.title    = element_text(family = "Arial", size = 24, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    plot.caption  = element_text(family = "Arial", size = 14, hjust = 0.5, margin = margin(t = 15)),
    axis.ticks.x  = element_blank(),
    axis.ticks.y  = element_line(color = "black", linewidth = 1 * PT),
    axis.line     = element_line(color = "black", linewidth = 1 * PT),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1 * PT),
    panel.grid    = element_blank(),
    legend.title  = element_text(size = 18),
    legend.text   = element_text(size = 18),
    legend.key    = element_rect(color = "black", linewidth = 0.75 * PT),
    legend.position = "right"
  )

# ---------- 9. Export ----------
out_dir <- here("results","figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(here(out_dir, "RSCU_heatmap.pdf"),
       p, width = 28, height = 10, device = cairo_pdf)

cat("\n===== RSCU heatmap completed =====\n")
cat("Output: results/figures/RSCU_heatmap.pdf\n")
