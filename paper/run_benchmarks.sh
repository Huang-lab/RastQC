#!/bin/bash
set -e

export PATH="/usr/local/opt/openjdk@11/bin:$HOME/.cargo/bin:$PATH"

RASTQC="/Users/huangk06/Projects/FastQC2/target/release/rastqc"
FASTQC_DIR="/Users/huangk06/Projects/FastQC2/FastQC"
FASTQC_CP="$FASTQC_DIR/bin:$FASTQC_DIR/sam-1.103.jar:$FASTQC_DIR/jbzip2-0.9.jar:$FASTQC_DIR/htsjdk.jar:$FASTQC_DIR/cisd-jhdf5.jar"
DATADIR="/Users/huangk06/Projects/FastQC2/paper/data"
OUTDIR="/Users/huangk06/Projects/FastQC2/paper/benchmarks"

mkdir -p "$OUTDIR/rastqc_out" "$OUTDIR/fastqc_out"

RESULTS="$OUTDIR/benchmark_results.tsv"
echo "Tool	Dataset	Reads	Threads	Median_s	Min_s	Max_s	MaxRSS_KB" > "$RESULTS"

bench() {
    local tool="$1" dataset="$2" reads="$3" threads="$4"
    shift 4
    local cmd=("$@")

    echo "  $tool $dataset t=$threads ..."

    # Use /usr/bin/time for memory, hyperfine for timing
    local rss_file=$(mktemp)

    # Run once with /usr/bin/time to get peak RSS
    /usr/bin/time -l "${cmd[@]}" 2>"$rss_file" >/dev/null || true
    local rss=$(grep "maximum resident" "$rss_file" | awk '{print $1}')
    local rss_kb=$((${rss:-0} / 1024))
    rm -f "$rss_file"

    # Run 3 times with hyperfine for timing
    local hf_json=$(mktemp)
    hyperfine --warmup 0 --runs 3 --export-json "$hf_json" -- "${cmd[*]}" 2>/dev/null

    local median=$(python3 -c "import json; d=json.load(open('$hf_json')); print(f'{d[\"results\"][0][\"median\"]:.3f}')" 2>/dev/null || echo "0")
    local mintime=$(python3 -c "import json; d=json.load(open('$hf_json')); print(f'{d[\"results\"][0][\"min\"]:.3f}')" 2>/dev/null || echo "0")
    local maxtime=$(python3 -c "import json; d=json.load(open('$hf_json')); print(f'{d[\"results\"][0][\"max\"]:.3f}')" 2>/dev/null || echo "0")
    rm -f "$hf_json"

    echo "    -> ${median}s (RSS: ${rss_kb}KB)"
    echo "$tool	$dataset	$reads	$threads	$median	$mintime	$maxtime	$rss_kb" >> "$RESULTS"
}

echo "=== Performance Benchmarks: RastQC vs FastQC ==="
echo ""

# Clean outputs
rm -rf "$OUTDIR/rastqc_out"/* "$OUTDIR/fastqc_out"/* 2>/dev/null

# ---- 100K reads ----
echo "--- 100K synthetic reads (33 MB) ---"
bench RastQC 100K 100000 1 \
    "$RASTQC" -t 1 --quiet --nozip -o "$OUTDIR/rastqc_out" "$DATADIR/synthetic_100K.fastq"
bench FastQC 100K 100000 1 \
    java -Djava.awt.headless=true -Xmx512m -cp "$FASTQC_CP" \
    uk.ac.babraham.FastQC.FastQCApplication --threads 1 --outdir "$OUTDIR/fastqc_out" "$DATADIR/synthetic_100K.fastq"

# ---- 1M reads ----
echo "--- 1M synthetic reads (325 MB) ---"
bench RastQC 1M 1000000 1 \
    "$RASTQC" -t 1 --quiet --nozip -o "$OUTDIR/rastqc_out" "$DATADIR/synthetic_1M.fastq"
bench FastQC 1M 1000000 1 \
    java -Djava.awt.headless=true -Xmx512m -cp "$FASTQC_CP" \
    uk.ac.babraham.FastQC.FastQCApplication --threads 1 --outdir "$OUTDIR/fastqc_out" "$DATADIR/synthetic_1M.fastq"

# ---- 1M gzip ----
echo "--- 1M synthetic reads gzipped (149 MB) ---"
bench RastQC 1M_gz 1000000 1 \
    "$RASTQC" -t 1 --quiet --nozip -o "$OUTDIR/rastqc_out" "$DATADIR/synthetic_1M.fastq.gz"
bench FastQC 1M_gz 1000000 1 \
    java -Djava.awt.headless=true -Xmx512m -cp "$FASTQC_CP" \
    uk.ac.babraham.FastQC.FastQCApplication --threads 1 --outdir "$OUTDIR/fastqc_out" "$DATADIR/synthetic_1M.fastq.gz"

# ---- 10M reads ----
echo "--- 10M synthetic reads (3.2 GB) ---"
bench RastQC 10M 10000000 1 \
    "$RASTQC" -t 1 --quiet --nozip -o "$OUTDIR/rastqc_out" "$DATADIR/synthetic_10M.fastq"
bench FastQC 10M 10000000 1 \
    java -Djava.awt.headless=true -Xmx512m -cp "$FASTQC_CP" \
    uk.ac.babraham.FastQC.FastQCApplication --threads 1 --outdir "$OUTDIR/fastqc_out" "$DATADIR/synthetic_10M.fastq"

# ---- Real E. coli ----
if [ -f "$DATADIR/real_ecoli.fastq" ]; then
    echo "--- Real E. coli data (5M reads, 1.2 GB) ---"
    bench RastQC real_ecoli 5000000 1 \
        "$RASTQC" -t 1 --quiet --nozip -o "$OUTDIR/rastqc_out" "$DATADIR/real_ecoli.fastq"
    bench FastQC real_ecoli 5000000 1 \
        java -Djava.awt.headless=true -Xmx512m -cp "$FASTQC_CP" \
        uk.ac.babraham.FastQC.FastQCApplication --threads 1 --outdir "$OUTDIR/fastqc_out" "$DATADIR/real_ecoli.fastq"
fi

echo ""
echo "=== Startup overhead test (empty-ish run) ==="
# Measure JVM startup vs native startup
printf "@r\nA\n+\nI\n" > /tmp/tiny.fastq
bench RastQC startup 1 1 \
    "$RASTQC" -t 1 --quiet --nozip -o "$OUTDIR/rastqc_out" /tmp/tiny.fastq
bench FastQC startup 1 1 \
    java -Djava.awt.headless=true -Xmx512m -cp "$FASTQC_CP" \
    uk.ac.babraham.FastQC.FastQCApplication --threads 1 --outdir "$OUTDIR/fastqc_out" /tmp/tiny.fastq
rm -f /tmp/tiny.fastq

echo ""
echo "=== Binary size comparison ==="
echo "RastQC binary: $(ls -lh "$RASTQC" | awk '{print $5}')"
echo "FastQC JARs + classes: $(du -sh "$FASTQC_DIR/bin" "$FASTQC_DIR"/*.jar 2>/dev/null | tail -1 | awk '{print $1}')"

echo ""
echo "=== Complete ==="
echo "Results saved to: $RESULTS"
echo ""
cat "$RESULTS"
