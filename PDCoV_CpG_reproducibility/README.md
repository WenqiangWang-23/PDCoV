# PDCoV CpG and dinucleotide reproducibility scripts

This package contains **only the whole-genome dinucleotide composition and CpG temporal-analysis workflow**. It does not include NSP8 haplotype or clustering analyses.

## Directory structure

```text
PDCoV_CpG_reproducibility/
├── README.md
├── scripts/
│   ├── 01_run_compseq.sh
│   └── 02_summarize_compseq.R
├── data/
│   └── README.txt
└── manuscript_reference/
    └── CpG_methods_and_results.md
```

## Input files

Place the following files in `data/`:

1. `PDCoV_265.fas`: multi-FASTA containing the 265 non-recombinant PDCoV whole-genome records.
2. `Genome_year.txt`: sampling metadata with a sequence identifier column and a sampling-year column. Recommended tab-delimited format:

```text
sequence_id	year
ACCESSION_A	2010
ACCESSION_B	2012
```

Identifiers should correspond to the FASTA headers. Versionless accession matching is also attempted (for example, `ABC123.1` can match `ABC123`).

## Software

- EMBOSS 6.6.0, with `compseq` available on PATH
- R 4.3.0
- R packages: `ggplot2`, `mgcv`

Install the R packages if needed:

```r
install.packages(c("ggplot2", "mgcv"))
```

## Run

Run from the package root directory.

### Step 1. Calculate dinucleotide composition with EMBOSS

```bash
bash scripts/01_run_compseq.sh data/PDCoV_265.fas results/compseq
```

The script splits the multi-FASTA into one FASTA file per genome and runs:

```bash
compseq -sequence individual_genome.fasta -word 2 -outfile output.compseq -calcfreq
```

It preserves each individual compseq output and writes a manifest.

### Step 2. Summarize all 16 dinucleotides and perform CpG temporal analyses

```bash
Rscript scripts/02_summarize_compseq.R \
  results/compseq \
  results/compseq_summary.tsv \
  data/Genome_year.txt
```

## Main outputs

- `results/compseq/split_fasta/`: one FASTA file per genome
- `results/compseq/per_genome/`: raw EMBOSS compseq output for each genome
- `results/compseq/compseq_manifest.tsv`: sequence-to-output mapping
- `results/compseq_summary.tsv`: all 16 dinucleotides for each genome, including observed/expected statistics and O/E
- `results/temporal_analysis/CpG_OE_with_year.csv`: matched CpG O/E values and sampling years
- `annual_CpG_summary.tsv`: annual median, Q1, Q3, and sample size
- `leave_one_year_out.tsv`: Spearman correlation after excluding each year in turn
- `GAM_model_summary.tsv`: GAM EDF, adjusted R-squared, and deviance explained for k = 4, 6, 8, and 10
- `CpG_temporal_trend.png`, `CpG_annual_median_IQR.png`, `CpG_leave_one_year_out.png`, `CpG_GAM_k4.png`
- `analysis_report.txt`: statistical summary and R session information

## Important reproducibility notes

1. The temporal analysis includes genomes with unambiguous sampling years that successfully match the FASTA identifiers.
2. The R script uses a fixed random seed and 5,000 bootstrap replicates for the percentile confidence interval.
3. The parser checks for the EMBOSS compseq table header and stops if it cannot recognize the installed output format. Inspect a raw `.compseq` file before adapting the parser.
4. This workflow implements the method stated in the manuscript: EMBOSS `compseq`, word size 2, and `-calcfreq`. Do not replace the EMBOSS values with an independently calculated R O/E formula without validating equivalence.
5. The script reports the results obtained from the supplied input files. Values quoted in the manuscript should be updated only after confirming that the same sequence set, metadata, software versions, and filtering rules were used.
