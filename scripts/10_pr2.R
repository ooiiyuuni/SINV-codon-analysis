# ============================================================
# 10_pr2.R
# Purpose : PR2 plots (G3/(G3+C3) vs A3/(A3+U3)) for SINV CDS, NSP
#           and SP, colored by host group.
# Input   : data/SINV/cds.fasta, nsp.fasta, sp.fasta
#           results/tables/host_class_summary.csv
# Output  : results/figures/PR2_all.pdf
#           results/tables/PR2/PR2_coords_CDS.csv
#           results/tables/PR2/PR2_coords_NSP.csv
#           results/tables/PR2/PR2_coords_SP.csv
# Dependencies: here, Biostrings, tidyverse, patchwork
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Biostrings)
  library(tidyverse)
  library(patchwork)
})

fam <- "ArialMT"
tab_dir <- here("results", "tables")
fig_dir <- here("results", "figures")
out_tab <- here("results", "tables", "PR2")
dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)

summary_df <- read.csv(file.path(tab_dir, "host_class_summary.csv"), stringsAsFactors = FALSE)
map_list <- list()
for (i in 1:nrow(summary_df)) {
  accs <- unlist(strsplit(summary_df$All_accessions[i], "; "))
  map_list[[i]] <- data.frame(Accession = accs, Group = summary_df$Host_class[i])
}
meta <- bind_rows(map_list) %>% mutate(meta_row = row_number())

compute_pr2 <- function(fa_path) {
  dna <- readDNAStringSet(fa_path)
  res <- lapply(seq_along(dna), function(i) {
    s <- toupper(as.character(dna[[i]]))
    n <- nchar(s)
    cod <- substring(s, seq(1, n - 2, by = 3), seq(3, n, by = 3))
    t3 <- substr(cod, 3, 3)
    A3 <- mean(t3 == "A"); U3 <- mean(t3 == "T")
    G3 <- mean(t3 == "G"); C3 <- mean(t3 == "C")
    data.frame(Accession = names(dna)[i],
               x = G3 / (G3 + C3),
               y = A3 / (A3 + U3))
  })
  bind_rows(res) %>%
    mutate(Accession = sub("\\..*$", "", Accession))
}

regions <- c("CDS", "NSP", "SP")
files <- c(CDS = "cds.fasta", NSP = "nsp.fasta", SP = "sp.fasta")

dat_all <- bind_rows(lapply(regions, function(r) {
  d <- compute_pr2(here("data", "SINV", files[[r]]))
  d$Region <- r
  d
})) %>%
  inner_join(meta, by = "Accession") %>%
  mutate(Group = factor(Group, levels = c("Human", "Experimental_passage", "Other_vertebrate", "Bird", "Mosquito", "Other_arthropod", "Unknown")))

stopifnot(nrow(dat_all) == 315 * 3)
for (r in regions) {
  d <- dat_all %>% filter(Region == r)
  stopifnot(nrow(d) == 315)
  stopifnot(!any(is.na(d$Group)), !any(is.na(d$x)), !any(is.na(d$y)))
}

# ---- 统一的三图坐标轴范围（包含所有区段的所有点）----
x_global <- range(dat_all$x, na.rm = TRUE)
y_global <- range(dat_all$y, na.rm = TRUE)

# 稍微留点余量（2%），保证点不贴边
x_rng <- x_global + c(-1, 1) * diff(x_global) * 0.02
y_rng <- y_global + c(-1, 1) * diff(y_global) * 0.02

n_labels <- paste0(summary_df$Host_class, " (n=", summary_df$N_sequences, ")")
names(n_labels) <- summary_df$Host_class

grp_cols <- c(Human = "#E41A1C", Experimental_passage = "#681F72", Other_vertebrate = "#A65628",
              Bird = "#DDCE43", Mosquito = "#64A9E2", Other_arthropod = "#33A02C", Unknown = "#808080")

base_theme <- theme_classic(base_size = 10, base_family = fam) +
  theme(
    plot.title   = element_text(size = 24, face = "bold", hjust = 0.5, family = fam, margin = margin(b = 8)),
    axis.title.x = element_text(size = 24, face = "bold", family = fam, margin = margin(t = 10)),
    axis.title.y = element_text(size = 24, face = "bold", family = fam, margin = margin(r = 20)),
    axis.text    = element_text(size = 10, color = "black", family = fam),
    axis.line    = element_line(color = "black", linewidth = 0.5 * 0.352778),
    axis.ticks   = element_line(color = "black", linewidth = 0.5 * 0.352778),
    legend.title = element_text(size = 10, face = "bold", family = fam),
    legend.text  = element_text(size = 10, family = fam),
    plot.margin  = margin(8, 8, 8, 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5 * 0.352778)
  )

make_pr2 <- function(df, x_rng, y_rng, title, show_legend = FALSE) {
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_vline(xintercept = 0.5, color = "grey40", linewidth = 0.5, linetype = "dashed") +
    geom_hline(yintercept = 0.5, color = "grey40", linewidth = 0.5, linetype = "dashed") +
    geom_point(aes(color = Group), shape = 16, size = 2.4, alpha = 0.9) +
    scale_color_manual(values = grp_cols, name = "Host group", labels = n_labels) +
    scale_x_continuous(limits = x_rng, breaks = scales::breaks_pretty(n = 5), expand = c(0, 0)) +
    scale_y_continuous(limits = y_rng, breaks = scales::breaks_pretty(n = 5), expand = c(0, 0)) +
    coord_fixed(ratio = 1) +
    labs(title = title, x = "G3 / (G3 + C3)", y = "A3 / (A3 + U3)") +
    base_theme
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

plots <- list()
for (r in regions) {
  plots[[r]] <- make_pr2(dat_all %>% filter(Region == r), x_rng, y_rng, r,
                         show_legend = (r == "CDS"))
}

combined <- (plots[["CDS"]] | plots[["NSP"]] | plots[["SP"]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "PR2_all.pdf"), combined, width = 18, height = 6, device = "pdf")

# ---- Per-sequence tables ----
for (r in regions) {
  d <- dat_all %>% filter(Region == r) %>%
    arrange(meta_row) %>%
    select(Accession, Group, x, y) %>%
    mutate(x = sprintf("%.4f", x), y = sprintf("%.4f", y))

  mean_row <- data.frame(Accession = "Mean", Group = "",
                         x = sprintf("%.4f", mean(as.numeric(d$x))),
                         y = sprintf("%.4f", mean(as.numeric(d$y))))
  sd_row   <- data.frame(Accession = "SD", Group = "",
                         x = sprintf("%.4f", sd(as.numeric(d$x))),
                         y = sprintf("%.4f", sd(as.numeric(d$y))))

  out <- bind_rows(d, mean_row, sd_row)
  write.csv(out, file.path(out_tab, paste0("PR2_coords_", r, ".csv")),
            row.names = FALSE, quote = FALSE)
}

cat("Done.\n")
cat("Figure:", file.path(fig_dir, "PR2_all.pdf"), "\n")
cat("Tables:", out_tab, "\n")
