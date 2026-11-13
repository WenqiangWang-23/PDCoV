#!/usr/bin/env Rscript

# Summarize per-genome EMBOSS compseq outputs and perform temporal CpG analyses.
# Usage:
#   Rscript scripts/02_summarize_compseq.R results/compseq results/compseq_summary.tsv data/Genome_year.txt
#
# The parser expects the standard EMBOSS compseq tabular section with columns
# corresponding to Word, Obs Count, Obs Frequency, Exp Count, Exp Frequency,
# and Obs/Exp. It stops with a diagnostic if the installed EMBOSS output format
# differs, rather than silently generating potentially incorrect values.

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
compseq_dir <- if (length(args) >= 1) args[[1]] else "results/compseq"
out_tsv <- if (length(args) >= 2) args[[2]] else file.path(compseq_dir, "compseq_summary.tsv")
year_file <- if (length(args) >= 3) args[[3]] else "data/Genome_year.txt"
plot_dir <- file.path(dirname(out_tsv), "temporal_analysis")
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

manifest_file <- file.path(compseq_dir, "compseq_manifest.tsv")
if (!file.exists(manifest_file)) stop("Manifest not found: ", manifest_file)
manifest <- read.delim(manifest_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("sequence_id", "compseq_output") %in% names(manifest))) {
  stop("Manifest must contain sequence_id and compseq_output columns.")
}

normalize_header <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

parse_one_compseq <- function(path, sequence_id) {
  if (!file.exists(path)) stop("Missing compseq output: ", path)
  lines <- readLines(path, warn = FALSE)
  header_i <- which(grepl("Word.*Obs.*Count.*Obs.*Frequency", lines, ignore.case = TRUE))[1]
  if (is.na(header_i)) {
    stop(
      "Could not identify the dinucleotide table header in: ", path,
      "\nInspect this file and adapt parse_one_compseq() to the installed EMBOSS output."
    )
  }

  tab_lines <- lines[(header_i + 1):length(lines)]
  tab_lines <- tab_lines[nzchar(trimws(tab_lines))]
  tab_lines <- tab_lines[grepl("^[ACGT]{2}[[:space:]]", tab_lines)]
  if (length(tab_lines) == 0) stop("No dinucleotide rows parsed from: ", path)

  dat <- read.table(text = tab_lines, header = FALSE, fill = TRUE,
                    stringsAsFactors = FALSE, comment.char = "")
  if (ncol(dat) < 6) {
    stop("Expected at least six columns in compseq table, found ", ncol(dat),
         " in: ", path)
  }
  names(dat)[1:6] <- c("Dinucleotide", "Observed_count",
                       "Observed_frequency", "Expected_count",
                       "Expected_frequency", "OE")
  dat$sequence_id <- sequence_id
  dat
}

parsed <- lapply(seq_len(nrow(manifest)), function(i) {
  parse_one_compseq(manifest$compseq_output[i], manifest$sequence_id[i])
})
all_dinuc <- do.call(rbind, parsed)
write.table(all_dinuc, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

# Keep all 16 dinucleotides in the complete output, then extract CpG (CG).
cpg <- all_dinuc[all_dinuc$Dinucleotide == "CG",
                 c("sequence_id", "OE", "Observed_count",
                   "Observed_frequency", "Expected_frequency")]
names(cpg)[names(cpg) == "OE"] <- "CpG_OE"

if (!file.exists(year_file)) {
  message("Dinucleotide summary written to: ", out_tsv)
  message("Year metadata not supplied/found; temporal analyses were skipped.")
  quit(status = 0)
}

# Year metadata: tab- or comma-delimited file with a sequence/accession column
# and a sampling-year column. Recommended headers: sequence_id, year.
read_metadata <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  sep <- if (grepl("\t", first)) "\t" else if (grepl(",", first)) "," else ""
  x <- read.table(path, header = TRUE, sep = sep, check.names = FALSE,
                  stringsAsFactors = FALSE, comment.char = "",
                  quote = "\"", fill = TRUE)
  nm <- tolower(gsub("[^a-z0-9]", "", names(x)))
  id_candidates <- which(nm %in% c("sequenceid", "accession", "id", "strain", "genomeid"))
  year_candidates <- which(nm %in% c("year", "samplingyear", "collectionyear", "date"))
  if (length(id_candidates) == 0 || length(year_candidates) == 0) {
    stop("Metadata needs recognizable ID and year columns. Use headers sequence_id and year.")
  }
  id <- as.character(x[[id_candidates[1]]])
  yr_raw <- as.character(x[[year_candidates[1]]])
  yr <- suppressWarnings(as.numeric(substr(yr_raw, 1, 4)))
  data.frame(sequence_id = id, year = yr, stringsAsFactors = FALSE)
}
meta <- read_metadata(year_file)

# Exact ID matching first; a versionless accession fallback is provided.
versionless <- function(x) sub("\\.[0-9]+$", "", x)
cpg$key <- versionless(trimws(cpg$sequence_id))
meta$key <- versionless(trimws(meta$sequence_id))
joined <- merge(cpg, meta[, c("key", "year")], by = "key")
joined <- joined[is.finite(joined$year) & is.finite(joined$CpG_OE), ]
if (nrow(joined) < 10) {
  stop("Fewer than 10 genomes matched year metadata. Check sequence identifiers.")
}
if (length(unique(joined$year)) < 3) {
  stop("At least three distinct sampling years are required for temporal analyses.")
}
write.table(joined, file.path(plot_dir, "CpG_OE_with_year.csv"),
            sep = ",", quote = TRUE, row.names = FALSE)

# Spearman correlation and percentile bootstrap CI.
set.seed(2026)
spearman <- cor.test(joined$year, joined$CpG_OE, method = "spearman", exact = FALSE)
boot_rho <- replicate(5000, {
  ix <- sample(seq_len(nrow(joined)), replace = TRUE)
  suppressWarnings(cor(joined$year[ix], joined$CpG_OE[ix],
                       method = "spearman", use = "complete.obs"))
})
boot_rho <- boot_rho[is.finite(boot_rho)]
rho_ci <- quantile(boot_rho, c(0.025, 0.975), na.rm = TRUE)

linear_fit <- lm(CpG_OE ~ year, data = joined)
linear_summary <- summary(linear_fit)

annual <- aggregate(CpG_OE ~ year, data = joined, function(x) {
  c(median = median(x), q1 = unname(quantile(x, 0.25)),
    q3 = unname(quantile(x, 0.75)), n = length(x))
})
annual_out <- data.frame(
  year = annual$year,
  median = annual$CpG_OE[, "median"],
  Q1 = annual$CpG_OE[, "q1"],
  Q3 = annual$CpG_OE[, "q3"],
  n = annual$CpG_OE[, "n"]
)
write.table(annual_out, file.path(plot_dir, "annual_CpG_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Leave-one-year-out sensitivity analysis.
years <- sort(unique(joined$year))
loo <- do.call(rbind, lapply(years, function(y) {
  d <- joined[joined$year != y, ]
  if (length(unique(d$year)) < 2 || nrow(d) < 4) return(NULL)
  ct <- suppressWarnings(cor.test(d$year, d$CpG_OE,
                                  method = "spearman", exact = FALSE))
  data.frame(excluded_year = y, rho = unname(ct$estimate),
             p_value = ct$p.value)
}))
write.table(loo, file.path(plot_dir, "leave_one_year_out.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# GAM sensitivity across k values. mgcv is part of the recommended R setup.
if (!requireNamespace("mgcv", quietly = TRUE)) {
  warning("Package 'mgcv' not installed; GAM analyses skipped.")
  gam_results <- data.frame()
} else {
  gam_results <- do.call(rbind, lapply(c(4, 6, 8, 10), function(k_value) {
    fit <- mgcv::gam(CpG_OE ~ s(year, k = k_value), data = joined, method = "REML")
    sm <- summary(fit)
    data.frame(k = k_value,
               EDF = unname(sm$s.table[1, "edf"]),
               adjusted_R2 = unname(sm$r.sq),
               deviance_explained_percent = 100 * unname(sm$dev.expl))
  }))
}
write.table(gam_results, file.path(plot_dir, "GAM_model_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p_year <- ggplot(joined, aes(x = year, y = CpG_OE)) +
  geom_point(alpha = 0.45, size = 1.5) +
  stat_summary(fun = median, geom = "line", aes(group = 1), linewidth = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  labs(x = "Sampling year", y = "CpG observed/expected ratio") +
  theme_classic()
ggsave(file.path(plot_dir, "CpG_temporal_trend.png"), p_year,
       width = 7, height = 5, dpi = 300)

p_annual <- ggplot(annual_out, aes(x = year, y = median)) +
  geom_ribbon(aes(ymin = Q1, ymax = Q3), alpha = 0.2) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(x = "Sampling year", y = "Annual median CpG O/E") +
  theme_classic()
ggsave(file.path(plot_dir, "CpG_annual_median_IQR.png"), p_annual,
       width = 7, height = 5, dpi = 300)

if (nrow(loo) > 0) {
  p_loo <- ggplot(loo, aes(x = factor(excluded_year), y = rho)) +
    geom_hline(yintercept = unname(spearman$estimate), linetype = 2) +
    geom_point(size = 2) +
    labs(x = "Excluded sampling year", y = "Spearman's rho") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(plot_dir, "CpG_leave_one_year_out.png"), p_loo,
         width = 8, height = 5, dpi = 300)
}

if (requireNamespace("mgcv", quietly = TRUE) && nrow(gam_results) > 0) {
  fit_final <- mgcv::gam(CpG_OE ~ s(year, k = 4), data = joined, method = "REML")
  pred_year <- data.frame(year = seq(min(joined$year), max(joined$year), length.out = 200))
  pr <- predict(fit_final, newdata = pred_year, se.fit = TRUE)
  pred <- data.frame(year = pred_year$year, fit = pr$fit,
                     lower = pr$fit - 1.96 * pr$se.fit,
                     upper = pr$fit + 1.96 * pr$se.fit)
  p_gam <- ggplot(joined, aes(x = year, y = CpG_OE)) +
    geom_point(alpha = 0.35) +
    geom_ribbon(data = pred, aes(x = year, ymin = lower, ymax = upper),
                inherit.aes = FALSE, alpha = 0.2) +
    geom_line(data = pred, aes(x = year, y = fit),
              inherit.aes = FALSE, linewidth = 0.9) +
    labs(x = "Sampling year", y = "CpG observed/expected ratio") +
    theme_classic()
  ggsave(file.path(plot_dir, "CpG_GAM_k4.png"), p_gam,
         width = 7, height = 5, dpi = 300)
}

report <- c(
  "CpG temporal analysis summary",
  paste0("Matched genomes: ", nrow(joined)),
  paste0("Sampling years: ", min(joined$year), "-", max(joined$year)),
  paste0("Spearman rho: ", round(unname(spearman$estimate), 4)),
  paste0("Spearman p-value: ", signif(spearman$p.value, 4)),
  paste0("Bootstrap 95% CI: [", round(rho_ci[1], 4), ", ", round(rho_ci[2], 4), "]"),
  paste0("Linear model slope: ", signif(coef(linear_fit)[["year"]], 5)),
  paste0("Linear model adjusted R-squared: ", round(linear_summary$adj.r.squared, 4)),
  "",
  "R session information:",
  capture.output(sessionInfo())
)
writeLines(report, file.path(plot_dir, "analysis_report.txt"))

message("Completed. Dinucleotide table: ", out_tsv)
message("Temporal outputs: ", plot_dir)
