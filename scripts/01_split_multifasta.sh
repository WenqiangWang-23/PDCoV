#!/usr/bin/env bash
# Split a multi-FASTA file into one FASTA file per sequence.
# Usage:
#   bash 01_split_multifasta.sh input.fasta output_directory
#
# Requirements: Bash and awk (available in Git Bash, Linux, and macOS).

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: bash $0 <input.fasta> <output_directory>" >&2
    exit 1
fi

INPUT="$1"
OUTDIR="$2"

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: Input FASTA file not found: $INPUT" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

# Use awk to preserve sequence lines and write one record per file.
# Header text is retained in each output FASTA. The filename is based on
# the first non-whitespace token after '>', with unsafe characters replaced.
awk -v outdir="$OUTDIR" '
function safe_name(s, t) {
    t = s
    sub(/^>/, "", t)
    sub(/[[:space:]].*$/, "", t)
    gsub(/[^A-Za-z0-9_.-]/, "_", t)
    if (t == "") t = "sequence"
    return t
}
function close_current() {
    if (outfile != "") close(outfile)
}
BEGIN {
    outfile = ""
    count = 0
}
 /^>/ {
    close_current()
    header = $0
    id = safe_name(header)
    outfile = outdir "/" id ".fasta"

    # Avoid silently overwriting files when IDs collide.
    if (system("test -e \"" outfile "\"") == 0) {
        suffix = 2
        while (system("test -e \"" outdir "/" id "_" suffix ".fasta\"") == 0) {
            suffix++
        }
        outfile = outdir "/" id "_" suffix ".fasta"
    }

    print header > outfile
    count++
    next
}
{
    if (outfile == "") {
        if ($0 ~ /[^[:space:]]/) {
            print "ERROR: Sequence data encountered before first FASTA header." > "/dev/stderr"
            exit 2
        }
    } else {
        print $0 >> outfile
    }
}
END {
    close_current()
    if (count == 0) {
        print "ERROR: No FASTA records found." > "/dev/stderr"
        exit 3
    }
    print "Successfully split " count " FASTA records into: " outdir
}
' "$INPUT"
