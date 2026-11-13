# PDCoV manuscript reproducibility scripts

This repository contains selected custom scripts supporting the manuscript
"Global Evolutionary Dynamics and Phylogeography of Porcine Deltacoronavirus".

## Repository structure

```text
PDCoV_reproducibility_scripts/
├── README.md
└── scripts/
    ├── 01_nsp8_haplotype_analysis.R
    └── 02_nsp8_clustering_analysis.R
```

The scripts perform two complementary analyses of aligned PDCoV NSP8
nucleotide sequences:

1. Haplotype inference, haplotype-network visualization, and summary statistics.
2. Sequence-level grouping using hierarchical clustering with silhouette
   evaluation and Louvain community detection.

The scripts do not perform sequence alignment. Use an aligned NSP8 FASTA
containing sequences from the final non-recombinant dataset (265 genome
records in the current manuscript analysis). Confirm that the FASTA and
metadata correspond to the exact dataset used for the manuscript.

## Requirements

R (version 4.0 or later recommended) and the following packages:

- ape
- pegas
- Biostrings
- dplyr
- tidyr
- ggplot2
- stringr
- cluster
- igraph

The scripts install missing packages from CRAN/Bioconductor where possible.
For reproducible publication analyses, record the R and package versions used.

## Input files

Place the input files in a `data/` directory under the repository root:

- `data/nsp8_aligned.fasta`: aligned NSP8 nucleotide sequences in FASTA format.
- `data/metadata.csv`: comma-separated metadata with at least these columns:
  - `Sequence_ID`: sequence identifiers matching the FASTA headers.
  - `Country`: sampling country for each sequence.

The scripts first attempt exact sequence-ID matching and then attempt matching
after removing GenBank version suffixes (for example, `.1`). They do not infer
missing countries, sampling dates, or other metadata.

If your metadata uses different column names, edit `id_column` and
`country_column` near the beginning of `01_nsp8_haplotype_analysis.R`.
The clustering script does not require metadata.

## 1. NSP8 haplotype analysis

Run from the repository root:

```bash
Rscript scripts/01_nsp8_haplotype_analysis.R
```

The script infers unique haplotypes from the aligned sequences, generates a
haplotype network colored by country, and exports sequence-to-haplotype and
country-frequency summaries.

Default output directory: `results/nsp8_haplotype/`

Main outputs:

- `Haplotype_network_by_country.png`
- `Haplotype_network_by_country.pdf`
- `Haplotype_summary.csv`
- `Sequence_to_haplotype.csv`
- `Haplotype_by_country.csv`
- `Haplotype_by_country_wide.csv`
- `Country_sample_counts.csv`
- `Diversity_statistics.csv`
- `Haplotype_sequences.fasta`
- `R_sessionInfo.txt`

The network is generated with `pegas::haploNet`. It represents mutational
connections among inferred haplotypes and should not be interpreted as a
phylogenetic tree or as evidence of transmission direction.

## 2. NSP8 sequence clustering

Run from the repository root:

```bash
Rscript scripts/02_nsp8_clustering_analysis.R
```

The script calculates pairwise nucleotide p-distances, evaluates hierarchical
clusters for k = 2 to 8 using average silhouette width, and performs Louvain
community detection on a thresholded sequence-similarity graph.

Default input: `data/nsp8_aligned.fasta`
Default output: `results/nsp8_clustering/hap_clusters.csv`

The distance returned by `ape::dist.dna(model = "N")` is a proportion
(p-distance), not an integer number of nucleotide differences. Therefore, the
default `threshold_bases <- 2` is converted to a proportion of the alignment
length before constructing the graph. Because pairwise deletion is enabled,
the effective number of compared sites can vary among sequence pairs; the
threshold should therefore be interpreted as an approximate alignment-length
scaled cutoff. Change `threshold_bases` only with a clear biological or
analytical justification.

The hierarchical-clustering and Louvain labels are algorithm-derived groups.
They are not automatically equivalent to the manuscript's manually defined
haplotype groups (A-D), and the number or membership of groups may differ
between methods.

Default output directory: `results/nsp8_clustering/`

Main outputs:

- `hap_clusters.csv`: sequence ID, hierarchical-clustering assignment, and
  Louvain community assignment.
- `silhouette_by_k.csv`: mean silhouette width for each tested k.
- `R_sessionInfo.txt`

## Reproducibility notes

- Both scripts require aligned sequences of equal length; neither script
  performs alignment.
- Check sequence identifiers and metadata matching before interpreting
  country-level summaries.
- The haplotype analysis is region-specific and does not replace whole-genome
  phylogenetic inference.
- The clustering analysis is exploratory. Its results should not be described
  as validated biological lineages without independent support.
- Do not upload restricted or unpublished sequence data unless permitted.
  Public accession numbers and metadata can be provided separately.

## Citation and versioning

For manuscript submission, create a tagged GitHub release or archive the
repository with Zenodo and cite the exact release/DOI in the manuscript.
Record the software and package versions used for the final analyses.
