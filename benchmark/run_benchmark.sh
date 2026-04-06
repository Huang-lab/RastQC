#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# RastQC Benchmark: Short-Read + Long-Read vs FastQC
# ============================================================

RASTQC="/Users/kuan-lin.huang/Projects/RastQC/target/release/rastqc"
DATADIR="/Users/kuan-lin.huang/Projects/RastQC/benchmark/data"
RESULTSDIR="/Users/kuan-lin.huang/Projects/RastQC/benchmark/results"
THREADS=4

mkdir -p "$RESULTSDIR/fastqc" "$RESULTSDIR/rastqc"

# ---- Classify files ----
SHORT_FILES=()
LONG_FILES=()
for f in "$DATADIR"/*.fastq.gz; do
    [ -f "$f" ] || continue
    fname=$(basename "$f")
    if [[ "$fname" == *_ont_* ]] || [[ "$fname" == *_pacbio_* ]]; then
        LONG_FILES+=("$f")
    else
        SHORT_FILES+=("$f")
    fi
done

TOTAL=$((${#SHORT_FILES[@]} + ${#LONG_FILES[@]}))
if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: No FASTQ files found in $DATADIR"
    exit 1
fi

echo "=============================================="
echo "  RastQC Benchmark"
echo "  Date: $(date)"
echo "  Threads: $THREADS"
echo "  Short-read files: ${#SHORT_FILES[@]}"
echo "  Long-read files:  ${#LONG_FILES[@]}"
echo "=============================================="
echo ""

# ---- CSV output ----
CSV="$RESULTSDIR/benchmark_results.csv"
echo "tool,file,type,size_mb,reads,real_sec,user_sec,sys_sec,max_rss_mb" > "$CSV"

# ---- Benchmark helper ----
bench() {
    local tool="$1"
    local label="$2"
    local ftype="$3"
    shift 3
    local cmd=("$@")

    echo "  Running $tool..."

    local timefile=$(mktemp)
    local start=$(python3 -c "import time; print(f'{time.time():.3f}')")

    /usr/bin/time -l "${cmd[@]}" > /dev/null 2> "$timefile"

    local end=$(python3 -c "import time; print(f'{time.time():.3f}')")
    local wall=$(python3 -c "print(f'{$end - $start:.2f}')")

    local user_sec=$(grep "user" "$timefile" | head -1 | awk '{print $1}')
    local sys_sec=$(grep "sys" "$timefile" | head -1 | awk '{print $1}')
    local max_rss=$(grep "maximum resident set size" "$timefile" | awk '{print $1}')
    local rss_mb=$((max_rss / 1048576))

    echo "    Wall: ${wall}s | MaxRSS: ${rss_mb}MB"

    local size_mb=""
    if [ -f "$DATADIR/$label" ]; then
        size_mb=$(( $(stat -f%z "$DATADIR/$label") / 1024 / 1024 ))
    fi

    echo "$tool,$label,$ftype,$size_mb,,$wall,$user_sec,$sys_sec,$rss_mb" >> "$CSV"
    rm -f "$timefile"
}

# ---- Per-file benchmarks ----
run_file() {
    local f="$1"
    local ftype="$2"
    local fname=$(basename "$f")

    echo ""
    echo "=== $fname ($ftype) ==="

    # Count reads
    echo "  Counting reads..."
    local nreads=$(gzip -dc "$f" 2>/dev/null | wc -l | awk '{print int($1/4)}')
    echo "  Reads: $nreads"

    # FastQC
    bench "fastqc" "$fname" "$ftype" fastqc -t "$THREADS" -o "$RESULTSDIR/fastqc" --quiet "$f"

    # RastQC (short-read mode)
    bench "rastqc" "$fname" "$ftype" "$RASTQC" -t "$THREADS" -o "$RESULTSDIR/rastqc" -q --time "$f"

    # RastQC with --long-read (only for long-read files)
    if [ "$ftype" = "long" ]; then
        bench "rastqc_lr" "$fname" "$ftype" "$RASTQC" --long-read -t "$THREADS" -o "$RESULTSDIR/rastqc" -q --time "$f"
    fi

    # Patch reads count into CSV
    sed -i '' "s/,$fname,$ftype,[^,]*,,/,$fname,$ftype,$(( $(stat -f%z "$f") / 1024 / 1024 )),$nreads,/g" "$CSV"
}

# Short-read files
if [ ${#SHORT_FILES[@]} -gt 0 ]; then
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  SHORT-READ BENCHMARKS"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    for f in "${SHORT_FILES[@]}"; do
        run_file "$f" "short"
    done

    # All short-read files together
    echo ""
    echo "=== ALL SHORT-READ FILES TOGETHER ==="
    bench "fastqc" "ALL_SHORT" "short" fastqc -t "$THREADS" -o "$RESULTSDIR/fastqc" --quiet "${SHORT_FILES[@]}"
    bench "rastqc" "ALL_SHORT" "short" "$RASTQC" -t "$THREADS" -o "$RESULTSDIR/rastqc" -q --time "${SHORT_FILES[@]}"
fi

# Long-read files
if [ ${#LONG_FILES[@]} -gt 0 ]; then
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  LONG-READ BENCHMARKS"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    for f in "${LONG_FILES[@]}"; do
        run_file "$f" "long"
    done
fi

echo ""
echo "=============================================="
echo "  Results saved to: $CSV"
echo "=============================================="

# ---- Summary ----
python3 << 'PYEOF'
import csv

results = {}
with open("/Users/kuan-lin.huang/Projects/RastQC/benchmark/results/benchmark_results.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        key = row["file"]
        tool = row["tool"]
        if key not in results:
            results[key] = {"type": row.get("type", ""), "reads": row.get("reads", "")}
        results[key][tool] = {
            "wall": float(row["real_sec"]),
            "size": row.get("size_mb", ""),
            "rss": row.get("max_rss_mb", ""),
        }

# Short-read table
print()
print("━━━ SHORT-READ RESULTS ━━━")
print(f"{'File':<30} {'Size':>6} {'Reads':>10} {'FastQC':>9} {'RastQC':>9} {'Speedup':>8} {'FQC RSS':>8} {'RQC RSS':>8}")
print("─" * 95)
for fname, d in results.items():
    if d["type"] != "short" or "fastqc" not in d or "rastqc" not in d:
        continue
    fqc = d["fastqc"]["wall"]
    rqc = d["rastqc"]["wall"]
    speedup = fqc / rqc if rqc > 0 else 0
    size = d.get("fastqc", d.get("rastqc", {})).get("size", "")
    reads = d.get("reads", "")
    fmem = d["fastqc"].get("rss", "")
    rmem = d["rastqc"].get("rss", "")
    print(f"{fname:<30} {size:>5}M {reads:>10} {fqc:>8.2f}s {rqc:>8.2f}s {speedup:>7.1f}x {fmem:>7}M {rmem:>7}M")

# Long-read table
long_entries = {k: v for k, v in results.items() if v["type"] == "long" and "fastqc" in v}
if long_entries:
    print()
    print("━━━ LONG-READ RESULTS ━━━")
    print(f"{'File':<35} {'Size':>6} {'Reads':>8} {'FastQC':>9} {'RastQC':>9} {'RQC+LR':>9} {'Speedup':>8}")
    print("─" * 95)
    for fname, d in long_entries.items():
        fqc = d["fastqc"]["wall"]
        rqc = d["rastqc"]["wall"]
        rqc_lr = d.get("rastqc_lr", {}).get("wall", 0)
        speedup = fqc / rqc if rqc > 0 else 0
        size = d.get("fastqc", d.get("rastqc", {})).get("size", "")
        reads = d.get("reads", "")
        lr_str = f"{rqc_lr:>8.2f}s" if rqc_lr > 0 else "       -"
        print(f"{fname:<35} {size:>5}M {reads:>8} {fqc:>8.2f}s {rqc:>8.2f}s {lr_str} {speedup:>7.1f}x")
PYEOF
