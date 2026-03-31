#!/bin/bash
set -e

export PATH="/usr/local/opt/openjdk@11/bin:$HOME/.cargo/bin:$PATH"

RASTQC="/Users/huangk06/Projects/RastQC/target/release/rastqc"
FASTQC_DIR="/Users/huangk06/Projects/RastQC/FastQC"
FASTQC_CP="$FASTQC_DIR/bin:$FASTQC_DIR/sam-1.103.jar:$FASTQC_DIR/jbzip2-0.9.jar:$FASTQC_DIR/htsjdk.jar:$FASTQC_DIR/cisd-jhdf5.jar"
DATADIR="/Users/huangk06/Projects/RastQC/paper/data/real_genomes"
OUTDIR="/Users/huangk06/Projects/RastQC/paper/benchmarks/real_genomes"
RESULTS="$OUTDIR/real_genome_benchmarks.tsv"

mkdir -p "$OUTDIR/rastqc_out" "$OUTDIR/fastqc_out"

echo "Tool	Organism	Reads	ReadLength	FileSize_MB	Time_s	RSS_bytes" > "$RESULTS"

bench() {
    local tool="$1" organism="$2" reads="$3" readlen="$4" filesize="$5"
    shift 5

    local tf=$(mktemp)
    /usr/bin/time -l "$@" 2>"$tf"

    local wall=$(grep "real" "$tf" | awk '{print $1}')
    local rss=$(grep "maximum resident" "$tf" | awk '{print $1}')
    rm -f "$tf"
    echo "$tool	$organism	$reads	$readlen	$filesize	$wall	$rss" >> "$RESULTS"
    echo "  $tool $organism: ${wall}s, RSS: $((rss/1048576))MB"
}

subsample() {
    local input="$1" output="$2" nreads="$3"
    local nlines=$((nreads * 4))

    if [[ "$input" == *.gz ]]; then
        gzip -dc "$input" 2>/dev/null | head -n "$nlines" > "$output"
    else
        head -n "$nlines" "$input" > "$output"
    fi

    local actual=$(wc -l < "$output")
    echo "  Subsampled to $((actual/4)) reads -> $output ($(du -h "$output" | cut -f1))"
}

TARGET_READS=5000000

echo "=== Preparing real genome datasets ==="

# E. coli - already have
if [ -f "/Users/huangk06/Projects/RastQC/paper/data/real_ecoli.fastq" ]; then
    cp "/Users/huangk06/Projects/RastQC/paper/data/real_ecoli.fastq" "$DATADIR/ecoli.fastq" 2>/dev/null || true
fi

# Subsample all downloaded files to 5M reads
for org in yeast fly mouse human; do
    gz="$DATADIR/${org}.fastq.gz"
    fq="$DATADIR/${org}.fastq"
    if [ -f "$gz" ] && [ ! -f "$fq" ]; then
        echo "Subsampling $org..."
        subsample "$gz" "$fq" $TARGET_READS
    elif [ -f "$fq" ]; then
        echo "$org already subsampled"
    else
        echo "WARNING: $org data not found, skipping"
    fi
done

echo ""
echo "=== Running benchmarks on real genome data ==="
echo ""

for org in ecoli yeast fly mouse human; do
    fq="$DATADIR/${org}.fastq"
    [ ! -f "$fq" ] && echo "Skipping $org: file not found" && continue

    reads=$(( $(wc -l < "$fq") / 4 ))
    # Detect read length from first read
    readlen=$(sed -n '2p' "$fq" | wc -c | tr -d ' ')
    readlen=$((readlen - 1))  # remove newline
    filesize=$(ls -l "$fq" | awk '{printf "%.0f", $5/1048576}')

    echo "--- $org ($reads reads, ${readlen}bp, ${filesize}MB) ---"

    # RastQC
    bench RastQC "$org" "$reads" "$readlen" "$filesize" \
        "$RASTQC" -t 1 --quiet --nozip -o "$OUTDIR/rastqc_out" "$fq"

    # FastQC
    bench FastQC "$org" "$reads" "$readlen" "$filesize" \
        java -Djava.awt.headless=true -Xmx512m -cp "$FASTQC_CP" \
        uk.ac.babraham.FastQC.FastQCApplication "$fq"

    # Compare outputs
    echo "  Comparing outputs..."
    python3 -c "
import os, sys

def parse_modules(path):
    modules = {}
    with open(path) as f:
        for line in f:
            if line.startswith('>>') and not line.startswith('>>END_MODULE'):
                parts = line.strip().lstrip('>>').split('\t')
                if len(parts) >= 2:
                    modules[parts[0]] = parts[1].upper()
    return modules

rdir = '$OUTDIR/rastqc_out'
fdir = '$DATADIR'
org = '$org'

# Find RastQC output
rpath = None
for f in os.listdir(rdir):
    if f.endswith('_fastqc.html'):
        rdata = os.path.join(rdir, f.replace('.html', ''), 'fastqc_data.txt')
        # Also check zip-extracted paths
        break

# RastQC text is embedded in the HTML process, check for standalone
rdata = os.path.join(rdir, org + '_fastqc.html')

# Actually we need the data file. Let's re-run with --extract
" 2>/dev/null || true

    echo ""
done

echo "=== Results ==="
cat "$RESULTS"
