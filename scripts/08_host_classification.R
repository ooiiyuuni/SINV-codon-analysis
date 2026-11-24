# scripts/08_host_classification.R
# Purpose: Assign SINV sequences to seven host classes and output a summary table.

metadata_path <- "results/tables/metadata.csv"
out_file <- "results/tables/host_class_summary.csv"

meta <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("Accession", "Host", "Isolate") %in% colnames(meta)))

host <- trimws(meta$Host)
host[is.na(host)] <- ""
isolate <- trimws(meta$Isolate)
isolate[is.na(isolate)] <- ""

is_experimental <- grepl("^Treatment", isolate) &
  host %in% c("Cricetulus griseus", "Mesocricetus auratus")

is_mosquito <- grepl("^(Culex|Aedes|Ochlerotatus|Anopheles|Culiseta|Coquillettidia|Mansonia|Mimomyia)",
                     host) |
  host %in% c("mosquito", "Culicidae sp.")

is_bird <- host %in% c("Streptopelia turtur", "Nycticorax nycticorax",
                       "Columba palumbus", "Somateria mollissima",
                       "Bubulcus ibis", "Corvus corone cornix")

is_other_vertebrate <- host %in% c("Cricetus cricetus", "Pelophylax ridibundus",
                                   "Chinese soft-shelled turtle")

is_other_arthropod <- host %in% c("Culicoides sp.", "midge",
                                  "phlebotomine sand flies", "Hyalomma marginatum")

host_class <- rep(NA_character_, nrow(meta))
host_class[host == "Homo sapiens"] <- "Human"
host_class[is_experimental] <- "Experimental_passage"
host_class[is_other_vertebrate] <- "Other_vertebrate"
host_class[is_bird] <- "Bird"
host_class[is_mosquito] <- "Mosquito"
host_class[is_other_arthropod] <- "Other_arthropod"
host_class[host == ""] <- "Unknown"

unassigned <- which(is.na(host_class))
if (length(unassigned) > 0) {
  cat("WARNING: unassigned host values:\n")
  print(unique(host[unassigned]))
  stop("Unassigned hosts detected. Update mapping before proceeding.")
}

class_order <- c("Human", "Experimental_passage", "Other_vertebrate",
                 "Bird", "Mosquito", "Other_arthropod", "Unknown")
host_class <- factor(host_class, levels = class_order)

summary_list <- lapply(class_order, function(cl) {
  idx <- which(host_class == cl)
  sub_hosts <- sort(unique(host[idx]))
  sub_accs <- meta$Accession[idx]
  data.frame(
    Host_class = cl,
    N_sequences = length(idx),
    Host_values = paste(sub_hosts, collapse = "; "),
    All_accessions = paste(sub_accs, collapse = "; "),
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, out_file, row.names = FALSE)

cat("\nHost class counts:\n")
print(table(host_class))
cat("\nSummary table written to:", out_file, "\n")
