#!/usr/bin/env Rscript

# PDCoV NSP8 haplotype analysis
# Input: aligned NSP8 FASTA + metadata CSV (Sequence_ID, Country)
# Run from repository root:
# Rscript scripts/01_nsp8_haplotype_analysis.R

required_cran <- c("ape", "pegas", "dplyr", "tidyr", "ggplot2", "stringr")
required_bioc <- c("Biostrings")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(ape)
  library(pegas)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# ---- User configuration ----
fasta_file <- "data/nsp8_aligned.fasta"
metadata_file <- "data/metadata.csv"
outdir <- "results/nsp8_haplotype"

id_column <- "Sequence_ID"
country_column <- "Country"
minimum_metadata_match_fraction <- 0.95
# ----------------------------

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(fasta_file)) stop("FASTA file not found: ", fasta_file)
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)

dna <- ape::read.dna(fasta_file, format = "fasta")
if (nrow(dna) < 2) stop("At least two sequences are required.")
if (length(unique(nchar(as.character(dna)))) != 1) {
  stop("Sequences must have equal aligned length.")
}

seq_ids <- rownames(dna)
if (is.null(seq_ids) || anyDuplicated(seq_ids)) {
  stop("FASTA sequence identifiers must be present and unique.")
}

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c(id_column, country_column)
if (!all(required_cols %in% names(metadata))) {
  stop("Metadata must contain columns: ", paste(required_cols, collapse = ", "))
}
metadata <- metadata %>%
  transmute(
    Sequence_ID = as.character(.data[[id_column]]),
    Country = as.character(.data[[country_column]])
  ) %>%
  distinct(Sequence_ID, .keep_all = TRUE)

strip_version <- function(x) sub("\\.[0-9]+$", "", x)
match_idx <- match(seq_ids, metadata$Sequence_ID)
unmatched <- is.na(match_idx)
if (any(unmatched)) {
  match_idx[unmatched] <- match(
    strip_version(seq_ids[unmatched]),
    strip_version(metadata$Sequence_ID)
  )
}

match_fraction <- mean(!is.na(match_idx))
if (match_fraction < minimum_metadata_match_fraction) {
  stop(
    sprintf("Only %.1f%% of FASTA records matched metadata; required minimum is %.1f%%.",
            100 * match_fraction, 100 * minimum_metadata_match_fraction)
  )
}

country <- rep(NA_character_, length(seq_ids))
country[!is.na(match_idx)] <- metadata$Country[match_idx[!is.na(match_idx)]]
if (anyNA(country) || any(country == "")) {
  warning(sum(is.na(country) | country == ""),
          " sequence(s) have missing country metadata; these will be labeled 'Unknown'.")
  country[is.na(country) | country == ""] <- "Unknown"
}

# Infer haplotypes from aligned nucleotide sequences.
hap <- pegas::haplotype(dna)
hap_matrix <- as.matrix(hap)
hap_labels <- paste0("H", seq_len(nrow(hap_matrix)))
rownames(hap_matrix) <- hap_labels

# Map each input sequence to its inferred haplotype.
hap_index <- attr(hap, "index")
if (is.null(hap_index)) stop("Could not retrieve sequence-to-haplotype indices.")
sequence_to_haplotype <- data.frame(
  Sequence_ID = seq_ids,
  Haplotype = NA_character_,
  Country = country,
  stringsAsFactors = FALSE
)
for (i in seq_along(hap_index)) {
  sequence_to_haplotype$Haplotype[hap_index[[i]]] <- hap_labels[i]
}
if (anyNA(sequence_to_haplotype$Haplotype)) {
  stop("Some sequences could not be assigned to a haplotype.")
}

hap_summary <- sequence_to_haplotype %>%
  count(Haplotype, name = "Frequency") %>%
  arrange(desc(Frequency), Haplotype)

hap_country <- sequence_to_haplotype %>%
  count(Haplotype, Country, name = "Count") %>%
  arrange(Haplotype, Country)

hap_country_wide <- hap_country %>%
  pivot_wider(names_from = Country, values_from = Count, values_fill = 0)

country_counts <- sequence_to_haplotype %>%
  count(Country, name = "Sequence_count") %>%
  arrange(desc(Sequence_count))

# Haplotype diversity and nucleotide diversity.
n <- nrow(dna)
hap_freq <- table(sequence_to_haplotype$Haplotype)
hap_diversity <- if (n > 1) {
  1 - sum(hap_freq * (hap_freq - 1)) / (n * (n - 1))
} else {
  NA_real_
}
dist_matrix <- as.matrix(ape::dist.dna(
  dna, model = "N", pairwise.deletion = TRUE
))
upper_dist <- dist_matrix[upper.tri(dist_matrix)]
nucleotide_diversity <- if (length(upper_dist)) mean(upper_dist, na.rm = TRUE) else NA_real_

diversity_stats <- data.frame(
  Number_of_sequences = n,
  Number_of_haplotypes = length(hap_freq),
  Haplotype_diversity_Hd = hap_diversity,
  Mean_pairwise_p_distance = nucleotide_diversity,
  Metadata_match_fraction = match_fraction
)

# Construct haplotype network.
hap_dist <- pegas::dist.haplotype(hap)
net <- pegas::haploNet(hap)
country_levels <- sort(unique(country))
country_palette <- grDevices::hcl.colors(length(country_levels), palette = "Dark 3")
names(country_palette) <- country_levels

hap_country_matrix <- table(
  factor(sequence_to_haplotype$Haplotype, levels = hap_labels),
  factor(sequence_to_haplotype$Country, levels = country_levels)
)
pie_values <- lapply(seq_len(nrow(hap_country_matrix)), function(i) {
  as.numeric(hap_country_matrix[i, ])
})

png(file.path(outdir, "Haplotype_network_by_country.png"),
    width = 2400, height = 1800, res = 300)
plot(net, size = attr(net, "freq"), scale = TRUE,
     pie = pie_values, bg = country_palette,
     show.mutation = 1, labels = TRUE,
     main = "PDCoV NSP8 haplotype network")
legend("topleft", legend = country_levels, fill = country_palette,
       cex = 0.75, bty = "n", title = "Country")
dev.off()

pdf(file.path(outdir, "Haplotype_network_by_country.pdf"),
    width = 8, height = 6)
plot(net, size = attr(net, "freq"), scale = TRUE,
     pie = pie_values, bg = country_palette,
     show.mutation = 1, labels = TRUE,
     main = "PDCoV NSP8 haplotype network")
legend("topleft", legend = country_levels, fill = country_palette,
       cex = 0.75, bty = "n", title = "Country")
dev.off()

# Write inferred haplotype sequences in FASTA format.
hap_dna <- as.DNAbin(hap_matrix)
rownames(hap_dna) <- hap_labels
ape::write.dna(hap_dna, file.path(outdir, "Haplotype_sequences.fasta"),
               format = "fasta", nbcol = -1, colsep = "")

write.csv(hap_summary, file.path(outdir, "Haplotype_summary.csv"), row.names = FALSE)
write.csv(sequence_to_haplotype, file.path(outdir, "Sequence_to_haplotype.csv"), row.names = FALSE)
write.csv(hap_country, file.path(outdir, "Haplotype_by_country.csv"), row.names = FALSE)
write.csv(hap_country_wide, file.path(outdir, "Haplotype_by_country_wide.csv"), row.names = FALSE)
write.csv(country_counts, file.path(outdir, "Country_sample_counts.csv"), row.names = FALSE)
write.csv(diversity_stats, file.path(outdir, "Diversity_statistics.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(outdir, "R_sessionInfo.txt"))
cat("NSP8 haplotype analysis completed. Results saved to: ", outdir, "\n")
