# PDCoV manuscript reproducibility scripts

This repository contains custom scripts for selected analyses in the manuscript
"Global Spread and Evolutionary Dynamics of Porcine Deltacoronavirus".

## Repository structure

```text
PDCoV_reproducibility_scripts/
├── README.md
├── scripts/
│   ├── 01_split_multifasta.sh
│   └── 02_nsp8_haplotype_analysis.R
└── data_template/
    └── metadata_template.csv
```

The manuscript's Figure 3 dinucleotide/CpG R script can be added as
`scripts/03_dinucleotide_CpG_analysis.R` after checking that the public version
uses relative paths and matches the final manuscript dataset.

## 1. Split a multi-FASTA file

Requirements: Bash and awk. On Windows, Git Bash or WSL can be used.

```bash
bash scripts/01_split_multifasta.sh input.fasta split_sequences/
```

The script creates one FASTA file per record. The original FASTA header is
preserved in each output file. Output filenames are based on the first token
of the sequence identifier.

## 2. nsp8 haplotype analysis

### Required input files

1. An aligned nsp8 FASTA file. All sequences must have the same alignment
   length. The script does not perform sequence alignment.
2. A comma-separated metadata file containing at least two columns:
   `Sequence_ID` and `Country`.

A metadata template is provided in `data_template/metadata_template.csv`.

### R packages

- ape
- pegas
- Biostrings
- dplyr
- tidyr
- ggplot2
- stringr

The script installs missing packages when necessary and writes `sessionInfo()`
to the output directory.

### Configure paths

Open `scripts/02_nsp8_haplotype_analysis.R` and edit:

```r
fasta_file <- "data/nsp8_aligned.fasta"
metadata_file <- "data/metadata.csv"
outdir <- "results/nsp8_haplotype"
```

If your metadata columns have different names, change `id_column` and
`country_column` accordingly.

### Run

From the repository root:

```bash
Rscript scripts/02_nsp8_haplotype_analysis.R
```

### Outputs

- `Haplotype_network_by_country.png` and `.pdf`
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
phylogenetic tree.

## Important reproducibility notes

- Confirm that the nsp8 alignment and metadata correspond to the final
  manuscript dataset before using outputs in the paper.
- The script matches sequence identifiers exactly first, then attempts a
  secondary match after removing GenBank version suffixes.
- Sequences without matching country metadata are excluded only if the
  match fraction is at least the configured threshold. The default threshold
  is 95%; set it to 1 to require every sequence to match.
- The script does not infer missing sampling locations or dates.
- The nsp8 analysis is region-specific and does not replace whole-genome
  phylogenetic inference.
- Do not upload restricted or unpublished sequence data unless permitted.
  Public accession numbers and metadata can be provided separately.

## Citation and versioning

For manuscript submission, create a tagged GitHub release or archive the
repository with Zenodo and cite the exact release/DOI in the manuscript.
Record the software and package versions used for the final analyses.
