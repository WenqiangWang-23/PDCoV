#!/usr/bin/env bash
set -euo pipefail

# Run EMBOSS compseq (v6.6.0) separately for each genome in a multi-FASTA file.
# Usage:
#   bash scripts/01_run_compseq.sh data/PDCoV_265.fas results/compseq
#
# Requirements:
#   EMBOSS compseq available on PATH (tested workflow target: v6.6.0).

FASTA="${1:-data/PDCoV_265.fas}"
OUTDIR="${2:-results/compseq}"
WORD_SIZE=2

if ! command -v compseq >/dev/null 2>&1; then
  echo "ERROR: compseq was not found. Install EMBOSS and ensure compseq is on PATH." >&2
  exit 1
fi
if [[ ! -s "$FASTA" ]]; then
  echo "ERROR: FASTA file does not exist or is empty: $FASTA" >&2
  exit 1
fi

mkdir -p "$OUTDIR/per_genome" "$OUTDIR/split_fasta"

# Split FASTA records without changing sequence identifiers or nucleotide content.
awk -v splitdir="$OUTDIR/split_fasta" '
  /^>/ {
    if (out != "") close(out)
    id = substr($0, 2)
    # Use the first non-whitespace token as a filesystem-safe record name.
    split(id, fields, /[ \t]/)
    name = fields[1]
    gsub(/[^A-Za-z0-9_.-]/, "_", name)
    if (name == "") name = "sequence_" ++n
    file = splitdir "/" name ".fasta"
    # Avoid overwriting if two headers share the same first token.
    if (system("test -e \"" file "\"") == 0) {
      file = splitdir "/" name "_" ++n ".fasta"
    }
    out = file
    print $0 > out
    next
  }
  {
    if (out != "") print $0 >> out
  }
' "$FASTA"

shopt -s nullglob
files=("$OUTDIR"/split_fasta/*.fasta)
if (( ${#files[@]} == 0 )); then
  echo "ERROR: No FASTA records were split from $FASTA" >&2
  exit 1
fi

: > "$OUTDIR/compseq_manifest.tsv"
printf "file\tsequence_id\tcompseq_output\n" >> "$OUTDIR/compseq_manifest.tsv"

for fasta in "${files[@]}"; do
  stem="$(basename "$fasta" .fasta)"
  out="$OUTDIR/per_genome/${stem}.compseq"
  compseq -sequence "$fasta" -word "$WORD_SIZE" -outfile "$out" -calcfreq
  header="$(grep -m1 '^>' "$fasta" | sed 's/^>//')"
  printf "%s\t%s\t%s\n" "$(basename "$fasta")" "$header" "$out" >> "$OUTDIR/compseq_manifest.tsv"
done

echo "Finished compseq for ${#files[@]} genomes."
echo "Raw outputs: $OUTDIR/per_genome/"
echo "Manifest: $OUTDIR/compseq_manifest.tsv"
echo "Next: Rscript scripts/02_summarize_compseq.R $OUTDIR $OUTDIR/compseq_summary.tsv"
