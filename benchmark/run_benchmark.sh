#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# RastQC benchmark: vs FastQC and Falco, short-read + long-read
# ============================================================
#
# Measures wall time and peak resident memory for each tool on every FASTQ in
# the data directory, repeating each measurement and reporting the median.
#
# FastQC and Falco are optional: whichever are found on PATH (or pointed at by
# $FASTQC / $FALCO) are included, and the rest are skipped. Only RastQC is
# required.
#
#   ./benchmark/run_benchmark.sh
#   THREADS=4 REPS=5 ./benchmark/run_benchmark.sh
#   DATADIR=/path/to/fastqs ./benchmark/run_benchmark.sh
#
# Note on threads: Falco is single-threaded (its -t flag is documented as "NOT
# YET IMPLEMENTED"), and FastQC's -t parallelizes across files rather than
# within one. So a single-file row compares one Falco core against $THREADS
# RastQC cores; the "rastqc -t 1" row is there for a same-core comparison.

# Resolve paths from this script's own location so the benchmark runs from any
# checkout. Each may be overridden from the environment.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

RASTQC="${RASTQC:-$REPO_ROOT/target/release/rastqc}"
FASTQC="${FASTQC:-$(command -v fastqc || true)}"
FALCO="${FALCO:-$(command -v falco || true)}"
DATADIR="${DATADIR:-$REPO_ROOT/benchmark/data}"
RESULTSDIR="${RESULTSDIR:-$REPO_ROOT/benchmark/results}"
THREADS="${THREADS:-4}"
REPS="${REPS:-3}"

if [ ! -x "$RASTQC" ]; then
    echo "ERROR: rastqc binary not found at $RASTQC" >&2
    echo "Build it first: cargo build --release" >&2
    exit 1
fi

mkdir -p "$RESULTSDIR"/{fastqc,falco,rastqc}

# ---- Portable peak-RSS measurement ----------------------------------------
# BSD/macOS `time -l` reports "maximum resident set size" in bytes; GNU
# `time -v` reports "Maximum resident set size (kbytes)" in KB.
case "$(uname -s)" in
    Darwin) TIME_FLAG="-l"; RSS_DIVISOR=1048576 ;;
    *)      TIME_FLAG="-v"; RSS_DIVISOR=1024 ;;
esac
if ! /usr/bin/time $TIME_FLAG true >/dev/null 2>&1; then
    echo "WARNING: /usr/bin/time $TIME_FLAG unavailable; memory will report as 0" >&2
    TIME_FLAG=""
fi

# ---- Classify files ----
SHORT_FILES=()
LONG_FILES=()
shopt -s nullglob
for f in "$DATADIR"/*.fastq "$DATADIR"/*.fq "$DATADIR"/*.fastq.gz "$DATADIR"/*.fq.gz; do
    [ -f "$f" ] || continue
    fname=$(basename "$f")
    if [[ "$fname" == *_ont_* ]] || [[ "$fname" == *_pacbio_* ]]; then
        LONG_FILES+=("$f")
    else
        SHORT_FILES+=("$f")
    fi
done
shopt -u nullglob

TOTAL=$((${#SHORT_FILES[@]} + ${#LONG_FILES[@]}))
if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: No FASTQ files found in $DATADIR" >&2
    exit 1
fi

echo "=============================================="
echo "  RastQC benchmark"
echo "  Date:    $(date)"
echo "  Host:    $(uname -sm)"
echo "  Threads: $THREADS   Repetitions: $REPS"
echo "  RastQC:  $($RASTQC --version 2>&1 | head -1)"
echo "  FastQC:  ${FASTQC:-(not found, skipped)}"
echo "  Falco:   ${FALCO:-(not found, skipped)}"
echo "  Files:   ${#SHORT_FILES[@]} short-read, ${#LONG_FILES[@]} long-read"
echo "=============================================="

CSV="$RESULTSDIR/benchmark_results.csv"
echo "tool,file,type,size_mb,reads,real_sec,max_rss_mb,reps" > "$CSV"

# ---- Benchmark helper ----
# bench <tool-label> <file-label> <type> -- <command...>
bench() {
    local tool="$1" label="$2" ftype="$3"
    shift 4  # drop the "--" separator too
    local cmd=("$@")

    local walls=() rsss=()
    for _ in $(seq 1 "$REPS"); do
        local timefile
        timefile=$(mktemp)
        local start end
        start=$(date +%s.%N 2>/dev/null || python3 -c 'import time;print(time.time())')
        if [ -n "$TIME_FLAG" ]; then
            /usr/bin/time $TIME_FLAG "${cmd[@]}" >/dev/null 2>"$timefile" || true
        else
            "${cmd[@]}" >/dev/null 2>"$timefile" || true
        fi
        end=$(date +%s.%N 2>/dev/null || python3 -c 'import time;print(time.time())')

        walls+=("$(awk -v s="$start" -v e="$end" 'BEGIN{printf "%.2f", e-s}')")
        local rss_raw
        rss_raw=$(grep -i "maximum resident set size" "$timefile" 2>/dev/null \
                  | grep -oE '[0-9]+' | head -1 || true)
        rsss+=("$(awk -v r="${rss_raw:-0}" -v d="$RSS_DIVISOR" 'BEGIN{printf "%.0f", r/d}')")
        rm -f "$timefile"
    done

    local wall rss
    wall=$(printf '%s\n' "${walls[@]}" | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')
    rss=$(printf '%s\n' "${rsss[@]}" | sort -n | awk 'END{print $1}')

    local size_mb=""
    if [ -f "$DATADIR/$label" ]; then
        size_mb=$(awk -v b="$(wc -c < "$DATADIR/$label")" 'BEGIN{printf "%.0f", b/1048576}')
    fi

    printf '  %-24s %8ss  %6s MB\n' "$tool" "$wall" "$rss"
    echo "$tool,$label,$ftype,$size_mb,,$wall,$rss,$REPS" >> "$CSV"
}

# ---- Per-file benchmarks ----
run_file() {
    local f="$1" ftype="$2"
    local fname
    fname=$(basename "$f")

    echo ""
    echo "=== $fname ($ftype) ==="

    local nreads
    case "$fname" in
        *.gz) nreads=$(gzip -dc "$f" | wc -l | awk '{print int($1/4)}') ;;
        *)    nreads=$(wc -l < "$f" | awk '{print int($1/4)}') ;;
    esac
    echo "  reads: $nreads"

    [ -n "$FASTQC" ] && bench "fastqc" "$fname" "$ftype" -- \
        "$FASTQC" -t "$THREADS" -o "$RESULTSDIR/fastqc" --quiet "$f"
    [ -n "$FALCO" ] && bench "falco" "$fname" "$ftype" -- \
        "$FALCO" -o "$RESULTSDIR/falco" "$f"

    bench "rastqc -t 1" "$fname" "$ftype" -- \
        "$RASTQC" -t 1 -o "$RESULTSDIR/rastqc" -q "$f"
    bench "rastqc -t $THREADS" "$fname" "$ftype" -- \
        "$RASTQC" -t "$THREADS" -o "$RESULTSDIR/rastqc" -q "$f"

    if [ "$ftype" = "long" ]; then
        bench "rastqc --long-read" "$fname" "$ftype" -- \
            "$RASTQC" --long-read -t "$THREADS" -o "$RESULTSDIR/rastqc" -q "$f"
    fi

    # Patch the read count in for every row of this file.
    local size_mb
    size_mb=$(awk -v b="$(wc -c < "$f")" 'BEGIN{printf "%.0f", b/1048576}')
    awk -F, -v OFS=, -v fn="$fname" -v n="$nreads" -v sz="$size_mb" \
        '$2==fn { $4=sz; $5=n } { print }' "$CSV" > "$CSV.tmp" && mv "$CSV.tmp" "$CSV"
}

if [ ${#SHORT_FILES[@]} -gt 0 ]; then
    echo ""
    echo "━━━ SHORT-READ BENCHMARKS ━━━"
    for f in "${SHORT_FILES[@]}"; do run_file "$f" "short"; done

    if [ ${#SHORT_FILES[@]} -gt 1 ]; then
        echo ""
        echo "=== ALL SHORT-READ FILES TOGETHER ==="
        [ -n "$FASTQC" ] && bench "fastqc" "ALL_SHORT" "short" -- \
            "$FASTQC" -t "$THREADS" -o "$RESULTSDIR/fastqc" --quiet "${SHORT_FILES[@]}"
        [ -n "$FALCO" ] && bench "falco" "ALL_SHORT" "short" -- \
            "$FALCO" -o "$RESULTSDIR/falco" "${SHORT_FILES[@]}"
        bench "rastqc -t $THREADS" "ALL_SHORT" "short" -- \
            "$RASTQC" -t "$THREADS" -o "$RESULTSDIR/rastqc" -q "${SHORT_FILES[@]}"
    fi
fi

if [ ${#LONG_FILES[@]} -gt 0 ]; then
    echo ""
    echo "━━━ LONG-READ BENCHMARKS ━━━"
    for f in "${LONG_FILES[@]}"; do run_file "$f" "long"; done
fi

echo ""
echo "=============================================="
echo "  Results: $CSV"
echo "=============================================="

# ---- Summary ----
CSV="$CSV" python3 << 'PYEOF'
import csv, os, collections

rows = collections.defaultdict(dict)
meta = {}
with open(os.environ["CSV"]) as fh:
    for row in csv.DictReader(fh):
        rows[row["file"]][row["tool"]] = {
            "wall": float(row["real_sec"]),
            "rss": row["max_rss_mb"],
        }
        meta[row["file"]] = (row["type"], row["size_mb"], row["reads"])

ref_tools = ["fastqc", "falco"]
for ftype, heading in (("short", "SHORT-READ RESULTS"), ("long", "LONG-READ RESULTS")):
    entries = {f: d for f, d in rows.items() if meta[f][0] == ftype}
    if not entries:
        continue
    print()
    print(f"━━━ {heading} ━━━")
    for fname, tools in entries.items():
        _, size, reads = meta[fname]
        print(f"\n{fname}  ({size} MB, {reads} reads)")
        print(f"  {'tool':<22} {'wall':>9} {'peak RSS':>10}  {'vs fastqc':>10} {'vs falco':>9}")
        print("  " + "─" * 66)
        for tool, m in tools.items():
            cmp = []
            for ref in ref_tools:
                if ref in tools and tool != ref and m["wall"] > 0:
                    cmp.append(f"{tools[ref]['wall'] / m['wall']:.1f}x")
                else:
                    cmp.append("-")
            print(f"  {tool:<22} {m['wall']:>8.2f}s {m['rss']:>8} MB  "
                  f"{cmp[0]:>10} {cmp[1]:>9}")
PYEOF
