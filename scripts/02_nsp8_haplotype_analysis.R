############################################################
# PDCoV nsp8 haplotype analysis
#
# Inputs:
#   1. An aligned nsp8 FASTA file (all sequences must have equal length)
#   2. A metadata table with at least:
#        Sequence_ID, Country
#
# Outputs:
#   - Haplotype_network_by_country.png / .pdf
#   - Haplotype_summary.csv
#   - Sequence_to_haplotype.csv
#   - Haplotype_by_country.csv
#   - Diversity_statistics.csv
#   - Haplotype_sequences.fasta
#
# Packages: ape, pegas, Biostrings, dplyr, tidyr, ggplot2
#
# IMPORTANT:
# This is a reproducible reconstruction/template. Confirm the input
# alignment, country labels, and filtering choices match the final
# manuscript dataset before using its outputs in the paper.
############################################################

# ==========================================================
# 0. Install and load packages
# ==========================================================

cran_pkgs <- c("ape", "pegas", "dplyr", "tidyr", "ggplot2", "stringr")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(ape)
  library(pegas)
  library(Biostrings)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# ==========================================================
# 1. User settings: EDIT THESE PATHS
# ==========================================================

fasta_file <- "data/nsp8_aligned.fasta"
metadata_file <- "data/metadata.csv"
outdir <- "results/nsp8_haplotype"

# Metadata column names
id_column <- "Sequence_ID"
country_column <- "Country"

# Minimum fraction of sequences that must match metadata.
# Set to 1 to require a complete match.
minimum_match_fraction <- 0.95

# Plot settings
plot_width <- 10
plot_height <- 8
plot_dpi <- 600

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(fasta_file)) {
  stop("Cannot find FASTA file: ", fasta_file)
}
if (!file.exists(metadata_file)) {
  stop("Cannot find metadata file: ", metadata_file)
}

# ==========================================================
# 2. Read aligned FASTA and validate sequences
# ==========================================================

dna_set <- Biostrings::readDNAStringSet(fasta_file, format = "fasta")
if (length(dna_set) < 3) {
  stop("At least three sequences are required.")
}

seq_ids <- names(dna_set)
seq_ids <- sub("^>", "", seq_ids)
seq_ids <- sub("[[:space:]].*$", "", seq_ids)

if (anyDuplicated(seq_ids)) {
  stop("Duplicate FASTA sequence IDs detected. Resolve them before analysis.")
}

seq_strings <- toupper(as.character(dna_set))
names(seq_strings) <- seq_ids

seq_lengths <- nchar(seq_strings)
if (length(unique(seq_lengths)) != 1) {
  stop(
    "Sequences are not aligned: sequence lengths differ. ",
    "Please align the nsp8 sequences before running this script."
  )
}

cat("Sequences in FASTA:", length(seq_strings), "\n")
cat("Alignment length:", unique(seq_lengths), "bp\n")

# ==========================================================
# 3. Read metadata and match sequence IDs
# ==========================================================

meta <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

required_cols <- c(id_column, country_column)
missing_cols <- setdiff(required_cols, names(meta))
if (length(missing_cols) > 0) {
  stop(
    "Metadata is missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

meta <- meta %>%
  transmute(
    Sequence_ID = trimws(as.character(.data[[id_column]])),
    Country = trimws(as.character(.data[[country_column]]))
  ) %>%
  filter(
    !is.na(Sequence_ID), Sequence_ID != "",
    !is.na(Country), Country != ""
  ) %>%
  distinct(Sequence_ID, .keep_all = TRUE)

# Exact accession match first, then accession without version suffix.
idx <- match(seq_ids, meta$Sequence_ID)
missing_idx <- which(is.na(idx))

if (length(missing_idx) > 0) {
  fasta_base <- sub("\\.[0-9]+$", "", seq_ids[missing_idx])
  meta_base <- sub("\\.[0-9]+$", "", meta$Sequence_ID)

  for (i in seq_along(missing_idx)) {
    hits <- which(meta_base == fasta_base[i])
    if (length(hits) == 1) idx[missing_idx[i]] <- hits
  }
}

matched <- !is.na(idx)
match_fraction <- mean(matched)

cat("Metadata match fraction:", round(match_fraction, 4), "\n")

if (match_fraction < minimum_match_fraction) {
  stop(
    "Only ", sum(matched), "/", length(seq_ids),
    " sequences matched metadata (", round(match_fraction * 100, 1),
    "%). Check accession IDs and metadata."
  )
}

if (any(!matched)) {
  warning(sum(!matched), " sequences lack metadata and will be excluded.")
}

seq_strings <- seq_strings[matched]
seq_ids <- seq_ids[matched]
meta_use <- meta[idx[matched], , drop = FALSE]
meta_use$Sequence_ID <- seq_ids

# Keep sequence and metadata order identical.
stopifnot(identical(meta_use$Sequence_ID, seq_ids))

# ==========================================================
# 4. Prepare aligned DNA data and identify haplotypes
# ==========================================================

# Convert the aligned FASTA to ape DNAbin format.
tmp_fasta <- tempfile(fileext = ".fasta")
writeXStringSet(
  DNAStringSet(seq_strings, use.names = TRUE),
  filepath = tmp_fasta,
  format = "fasta"
)
dna <- ape::read.dna(tmp_fasta, format = "fasta", as.character = FALSE)
unlink(tmp_fasta)

# Haplotype assignment using pegas.
haps <- pegas::haplotype(dna)

hap_names <- rownames(haps)
if (is.null(hap_names)) {
  hap_names <- paste0("Hap_", seq_len(nrow(haps)))
  rownames(haps) <- hap_names
}

# pegas stores the original sequence indices belonging to each haplotype.
hap_index <- attr(haps, "index")
if (is.null(hap_index) || length(hap_index) != nrow(haps)) {
  stop("Could not retrieve sequence membership for inferred haplotypes.")
}

# Normalize membership list and create sequence-to-haplotype mapping.
seq_to_hap <- bind_rows(lapply(seq_along(hap_index), function(i) {
  members <- as.integer(hap_index[[i]])
  data.frame(
    Sequence_ID = seq_ids[members],
    Haplotype = hap_names[i],
    stringsAsFactors = FALSE
  )
}))

if (nrow(seq_to_hap) != length(seq_ids)) {
  stop("Haplotype assignment did not account for every sequence.")
}

# Frequency table
hap_freq <- seq_to_hap %>%
  count(Haplotype, name = "Frequency") %>%
  arrange(desc(Frequency), Haplotype)

hap_freq$Haplotype <- factor(
  hap_freq$Haplotype,
  levels = hap_names
)

# ==========================================================
# 5. Haplotype diversity and nucleotide diversity
# ==========================================================

N <- length(seq_ids)
hap_counts <- table(seq_to_hap$Haplotype)
hap_proportions <- as.numeric(hap_counts) / N

haplotype_diversity <- if (N > 1) {
  (N / (N - 1)) * (1 - sum(hap_proportions^2))
} else {
  NA_real_
}

nucleotide_diversity <- tryCatch(
  as.numeric(ape::nuc.div(dna, pairwise.deletion = TRUE)),
  error = function(e) NA_real_
)

diversity_stats <- data.frame(
  Region = "nsp8",
  Number_of_sequences = N,
  Alignment_length_bp = unique(seq_lengths),
  Number_of_haplotypes = length(hap_counts),
  Haplotype_diversity_Hd = haplotype_diversity,
  Nucleotide_diversity_pi = nucleotide_diversity
)

# ==========================================================
# 6. Haplotype-by-country composition
# ==========================================================

hap_country <- seq_to_hap %>%
  left_join(meta_use, by = "Sequence_ID") %>%
  count(Haplotype, Country, name = "Count") %>%
  arrange(Haplotype, Country)

country_totals <- meta_use %>%
  count(Country, name = "Total_sequences") %>%
  arrange(desc(Total_sequences))

# Wide table for easy inspection
hap_country_wide <- hap_country %>%
  pivot_wider(
    names_from = Country,
    values_from = Count,
    values_fill = 0
  )

hap_summary <- hap_freq %>%
  mutate(Haplotype = as.character(Haplotype)) %>%
  left_join(
    hap_country %>%
      group_by(Haplotype) %>%
      summarise(
        Countries = paste(sort(unique(Country)), collapse = "; "),
        .groups = "drop"
      ),
    by = "Haplotype"
  )

# ==========================================================
# 7. Construct median-joining-style haplotype network
# ==========================================================

# pegas::haploNet constructs a statistical parsimony haplotype network.
# This is not a phylogenetic tree and should be interpreted as a
# representation of mutational connections among the observed haplotypes.
net <- pegas::haploNet(haps)

# Country counts per haplotype, ordered by the haplotype labels.
country_levels <- sort(unique(meta_use$Country))
pie_matrix <- matrix(
  0,
  nrow = length(hap_names),
  ncol = length(country_levels),
  dimnames = list(hap_names, country_levels)
)

for (r in seq_len(nrow(hap_country))) {
  pie_matrix[
    as.character(hap_country$Haplotype[r]),
    as.character(hap_country$Country[r])
  ] <- hap_country$Count[r]
}

# Save network as high-resolution PNG and vector PDF.
country_palette <- grDevices::hcl.colors(
  max(3, length(country_levels)),
  palette = "Dark 3"
)[seq_along(country_levels)]

png(
  filename = file.path(outdir, "Haplotype_network_by_country.png"),
  width = plot_width,
  height = plot_height,
  units = "in",
  res = plot_dpi
)

par(mar = c(1, 1, 3, 1))
plot(
  net,
  size = attr(haps, "freq"),
  scale.ratio = 1,
  pie = pie_matrix,
  bg = country_palette,
  labels = TRUE,
  show.mutation = 1,
  font = 2,
  cex = 0.9,
  main = "PDCoV nsp8 haplotype network"
)
legend(
  "topleft",
  legend = country_levels,
  fill = country_palette,
  border = NA,
  bty = "n",
  cex = 0.85,
  title = "Country"
)
dev.off()

pdf(
  file = file.path(outdir, "Haplotype_network_by_country.pdf"),
  width = plot_width,
  height = plot_height
)
par(mar = c(1, 1, 3, 1))
plot(
  net,
  size = attr(haps, "freq"),
  scale.ratio = 1,
  pie = pie_matrix,
  bg = country_palette,
  labels = TRUE,
  show.mutation = 1,
  font = 2,
  cex = 0.9,
  main = "PDCoV nsp8 haplotype network"
)
legend(
  "topleft",
  legend = country_levels,
  fill = country_palette,
  border = NA,
  bty = "n",
  cex = 0.85,
  title = "Country"
)
dev.off()

# ==========================================================
# 8. Export haplotype sequences and result tables
# ==========================================================

# Convert inferred haplotypes to a FASTA file.
hap_strings <- as.character(haps)
names(hap_strings) <- hap_names
hap_dna <- Biostrings::DNAStringSet(hap_strings)
names(hap_dna) <- hap_names

Biostrings::writeXStringSet(
  hap_dna,
  filepath = file.path(outdir, "Haplotype_sequences.fasta"),
  format = "fasta"
)

write.csv(
  seq_to_hap,
  file.path(outdir, "Sequence_to_haplotype.csv"),
  row.names = FALSE
)

write.csv(
  hap_summary,
  file.path(outdir, "Haplotype_summary.csv"),
  row.names = FALSE
)

write.csv(
  hap_country,
  file.path(outdir, "Haplotype_by_country.csv"),
  row.names = FALSE
)

write.csv(
  hap_country_wide,
  file.path(outdir, "Haplotype_by_country_wide.csv"),
  row.names = FALSE
)

write.csv(
  country_totals,
  file.path(outdir, "Country_sample_counts.csv"),
  row.names = FALSE
)

write.csv(
  diversity_stats,
  file.path(outdir, "Diversity_statistics.csv"),
  row.names = FALSE
)

# Record package/session information for reproducibility.
capture.output(
  sessionInfo(),
  file = file.path(outdir, "R_sessionInfo.txt")
)

cat("\nAnalysis completed.\n")
cat("Sequences analyzed:", N, "\n")
cat("Haplotypes identified:", length(hap_counts), "\n")
cat("Haplotype diversity (Hd):", round(haplotype_diversity, 4), "\n")
cat("Nucleotide diversity (pi):", round(nucleotide_diversity, 6), "\n")
cat("Output directory:", normalizePath(outdir), "\n")
