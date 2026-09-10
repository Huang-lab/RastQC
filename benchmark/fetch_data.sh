#!/usr/bin/env bash
set -euo pipefail

# Download public sequencing data for benchmarking.
#
# The dataset committed under benchmark/data/ is deliberately tiny so the
# repository stays small, but 200k reads finish in well under a second — which
# is process startup, not throughput, and tells you nothing about how any of
# these tools behave on a real run. Fetch a real one before drawing
# conclusions.
#
#   ./benchmark/fetch_data.sh            # ~870 MB, 18.7M NextSeq 500 reads
#   ./benchmark/fetch_data.sh small      # ~56 MB, 1.2M NextSeq 500 reads
#
# Then: ./benchmark/run_benchmark.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATADIR="${DATADIR:-$SCRIPT_DIR/data}"
SET="${1:-default}"

# ENA read files. Both are Illumina NextSeq 500 WGS runs — the instrument and
# library type most short-read QC is actually run on.
case "$SET" in
    small)
        # DRR048760: 1.2M reads, 76 bp
        URLS=("https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR048/DRR048760/DRR048760.fastq.gz")
        ;;
    default)
        # DRR045135: 18.7M reads, 144 bp, paired — R1 only is plenty
        URLS=("https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR045/DRR045135/DRR045135_1.fastq.gz")
        ;;
    *)
        echo "Unknown dataset '$SET' (expected 'default' or 'small')" >&2
        exit 1
        ;;
esac

mkdir -p "$DATADIR"
for url in "${URLS[@]}"; do
    out="$DATADIR/$(basename "$url")"
    if [ -s "$out" ]; then
        echo "already present: $out"
        continue
    fi
    echo "downloading $(basename "$url") -> $DATADIR"
    # -C - resumes a partial download, so an interrupted fetch can be retried.
    curl -fL -C - -o "$out" "$url"
done

echo ""
echo "Data in $DATADIR:"
ls -lh "$DATADIR"/*.fastq.gz 2>/dev/null || true
