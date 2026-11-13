#!/usr/bin/env Rscript

# PDCoV NSP8 sequence clustering
# Run from repository root:
# Rscript scripts/02_nsp8_clustering_analysis.R

required_packages <- c("ape", "cluster", "igraph")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  library(ape)
  library(cluster)
  library(igraph)
})

# ---- User configuration ----
fasta_file <- "data/nsp8_aligned.fasta"
outdir <- "results/nsp8_clustering"
threshold_bases <- 2
k_values <- 2:8
# ----------------------------

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(fasta_file)) stop("FASTA file not found: ", fasta_file)

dna <- ape::read.dna(fasta_file, format = "fasta")
if (nrow(dna) < 3) stop("At least three sequences are required for this clustering analysis.")
seq_ids <- rownames(dna)
if (is.null(seq_ids) || anyDuplicated(seq_ids)) {
  stop("FASTA sequence identifiers must be present and unique.")
}

alignment_length <- ncol(dna)
if (length(unique(nchar(as.character(dna)))) != 1) {
  stop("Sequences must have equal aligned length.")
}

# ape::dist.dna(model = 'N') returns p-distance (a proportion).
distmat <- as.matrix(ape::dist.dna(
  dna, model = "N", pairwise.deletion = TRUE
))
if (any(!is.finite(distmat))) stop("Distance matrix contains non-finite values; check gaps and ambiguous bases.")

# ---- Method A: hierarchical clustering + silhouette analysis ----
hc <- hclust(as.dist(distmat), method = "average")
valid_k <- k_values[k_values >= 2 & k_values < nrow(dna)]
if (length(valid_k) == 0) stop("No valid k values for silhouette analysis.")

silhouette_means <- sapply(valid_k, function(k) {
  cl <- cutree(hc, k = k)
  mean(cluster::silhouette(cl, dist = as.dist(distmat))[, "sil_width"])
})
best_k <- valid_k[which.max(silhouette_means)]
groups_tree <- cutree(hc, k = best_k)

silhouette_df <- data.frame(
  k = valid_k,
  Mean_silhouette_width = as.numeric(silhouette_means)
)
write.csv(silhouette_df, file.path(outdir, "silhouette_by_k.csv"), row.names = FALSE)

# ---- Method B: threshold graph + Louvain community detection ----
# Convert the requested base-difference threshold to a p-distance cutoff.
# With pairwise deletion, this is an alignment-length-scaled approximation.
threshold_p <- threshold_bases / alignment_length
adj <- distmat <= threshold_p
diag(adj) <- FALSE

g <- igraph::graph_from_adjacency_matrix(
  adj, mode = "undirected", diag = FALSE
)
if (igraph::ecount(g) == 0) {
  stop("The threshold graph has no edges. Increase threshold_bases or review the alignment.")
}
if (igraph::components(g)$no == igraph::vcount(g)) {
  warning("All sequences are isolated under the selected threshold.")
}
comm <- igraph::cluster_louvain(g)
groups_net <- igraph::membership(comm)

cluster_df <- data.frame(
  Sequence_ID = seq_ids,
  TreeCluster = as.integer(groups_tree),
  NetCluster = as.integer(groups_net),
  stringsAsFactors = FALSE
)

write.csv(cluster_df, file.path(outdir, "hap_clusters.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "R_sessionInfo.txt"))

cat("Best number of hierarchical clusters (silhouette): ", best_k, "\n")
cat("Louvain communities: ", length(unique(groups_net)), "\n")
cat("Clustering results saved to: ", outdir, "\n")
