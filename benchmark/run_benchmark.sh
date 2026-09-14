#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# RastQC benchmark: vs FastQC and Falco, short-read + long-read
# ============================================================
#
# Measures wall time and peak resident memory for each tool on every FASTQ in
# the data directory, repeating each measurement and reporting the fastest wall
# time and the peak RSS across repetitions (see the note in bench() for why
# the fastest rather than the median).
#
# FastQC and Falco are optional: whichever are found on PATH (or pointed at by
# $FASTQC / $FALCO) are included, and the rest are skipped. Only RastQC is
# required.
#
#   ./benchmark/fetch_data.sh all        # get real data first
#   ./benchmark/run_benchmark.sh
#   THREADS=4 REPS=5 ./benchmark/run_benchmark.sh
#   DATADIR=/path/to/fastqs ./benchmark/run_benchmark.sh
#
# Writes benchmark/results/benchmark_results.csv, which is the single source
# for both benchmark/RESULTS.md and the paper's figures — regenerate those with
# `python3 paper/analyze_benchmarks.py`.
#
# Note on threads: Falco is single-threaded (its -t flag is documented in its
# own --help as "NOT YET IMPLEMENTED IN FALCO"), and FastQC's -t parallelizes
# across files rather than within one, so on a single file it also runs one
# core. The rastqc_t1 rows are therefore the like-for-like one-core comparison;
# the rastqc_tN rows show what RastQC does when given more cores.

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

# Mean read length above which a file is treated as long-read. ONT and PacBio
# runs are thousands of bases; no Illumina platform exceeds a few hundred. The
# gap is wide enough that this needs no per-platform configuration, and unlike
# a filename convention it cannot be defeated by how ENA happens to name a file.
LONG_READ_MIN_MEAN_LEN="${LONG_READ_MIN_MEAN_LEN:-1000}"

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

# ---- Per-file stats, computed once and cached ------------------------------
# Both the read count and the mean read length need a full pass over the file,
# which for a 1.5 GB gzipped input costs more than some of the benchmark runs
# it annotates. Cache them next to the data, keyed on the file's own size so a
# re-downloaded or truncated file recomputes rather than reusing stale numbers.
STATS_CACHE="$DATADIR/.benchmark_stats"
touch "$STATS_CACHE"

file_stats() {
    local f="$1"
    local bytes key hit
    bytes=$(wc -c <"$f" | tr -d ' ')
    key="$(basename "$f"):$bytes"

    hit=$(grep -F "$key " "$STATS_CACHE" 2>/dev/null | head -1 || true)
    if [ -n "$hit" ]; then
        echo "${hit#"$key" }"
        return
    fi

    local stats
    case "$f" in
        *.gz) stats=$(gzip -dc "$f" | awk 'NR%4==2 {n++; b+=length($0)} END{printf "%d %.0f", n, (n?b/n:0)}') ;;
        *)    stats=$(awk 'NR%4==2 {n++; b+=length($0)} END{printf "%d %.0f", n, (n?b/n:0)}' "$f") ;;
    esac
    echo "$key $stats" >> "$STATS_CACHE"
    echo "$stats"
}

# ---- Classify files by measured read length ----
SHORT_FILES=()
LONG_FILES=()

# Per-file stats as a newline-delimited "name reads meanlen size_mb" table.
# macOS ships bash 3.2, which has no associative arrays, and these scripts are
# meant to run on whatever bash the machine already has.
FILE_STATS=""

stat_of() {  # stat_of <name> <1=reads|2=meanlen|3=size_mb>
    printf '%s\n' "$FILE_STATS" | awk -v n="$1" -v c="$2" '$1==n {print $(c+1); exit}'
}

shopt -s nullglob
CANDIDATES=("$DATADIR"/*.fastq "$DATADIR"/*.fq "$DATADIR"/*.fastq.gz "$DATADIR"/*.fq.gz)
shopt -u nullglob

if [ ${#CANDIDATES[@]} -eq 0 ]; then
    echo "ERROR: No FASTQ files found in $DATADIR" >&2
    echo "Fetch real data first: ./benchmark/fetch_data.sh all" >&2
    exit 1
fi

echo "Scanning ${#CANDIDATES[@]} input file(s) for read count and mean length..."
for f in "${CANDIDATES[@]}"; do
    [ -f "$f" ] || continue
    fname=$(basename "$f")
    read -r nreads meanlen <<<"$(file_stats "$f")"
    size_mb=$(awk -v b="$(wc -c <"$f")" 'BEGIN{printf "%.0f", b/1048576}')
    FILE_STATS="$FILE_STATS$fname $nreads $meanlen $size_mb
"

    if [ "$meanlen" -ge "$LONG_READ_MIN_MEAN_LEN" ]; then
        LONG_FILES+=("$f")
        kind=long
    else
        SHORT_FILES+=("$f")
        kind=short
    fi
    printf '  %-34s %10s reads  %6s bp mean  %5s MB  -> %s\n' \
        "$fname" "$nreads" "$meanlen" "$size_mb" "$kind"
done

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
echo "tool,file,type,size_mb,reads,mean_read_len,real_sec,max_rss_mb,reps,threads" > "$CSV"

# ---- Benchmark helper ----
# bench <tool-id> <threads> <file-label> <type> -- <command...>
#
# <tool-id> is the stable identifier written to the CSV (fastqc, falco,
# rastqc_t1, rastqc_tN, rastqc_lr). paper/analyze_benchmarks.py matches on
# these, so they are part of this script's contract — renaming one silently
# drops that series from the figures.
bench() {
    local tool="$1" nthreads="$2" label="$3" ftype="$4"
    shift 5  # drop the "--" separator too
    local cmd=("$@")

    local walls=() rsss=()
    for _ in $(seq 1 "$REPS"); do
        local timefile
        timefile=$(mktemp)
        local start end
        start=$(date +%s.%N 2>/dev/null || python3 -c 'import time;print(time.time())')
        local rc=0
        if [ -n "$TIME_FLAG" ]; then
            /usr/bin/time $TIME_FLAG "${cmd[@]}" >/dev/null 2>"$timefile" || rc=$?
        else
            "${cmd[@]}" >/dev/null 2>"$timefile" || rc=$?
        fi
        end=$(date +%s.%N 2>/dev/null || python3 -c 'import time;print(time.time())')

        # A tool that rejected the input exits fast; recording that as a wall
        # time would report a failure as a speedup.
        if [ "$rc" -ne 0 ]; then
            printf '  %-24s FAILED (exit %s) — not recorded\n' "$tool" "$rc"
            sed -n '1,3p' "$timefile" | sed 's/^/      /'
            rm -f "$timefile"
            return
        fi

        walls+=("$(awk -v s="$start" -v e="$end" 'BEGIN{printf "%.2f", e-s}')")
        local rss_raw
        rss_raw=$(grep -i "maximum resident set size" "$timefile" 2>/dev/null \
                  | grep -oE '[0-9]+' | head -1 || true)
        rsss+=("$(awk -v r="${rss_raw:-0}" -v d="$RSS_DIVISOR" 'BEGIN{printf "%.0f", r/d}')")
        rm -f "$timefile"
    done

    # Fastest wall time, peak RSS across repetitions.
    #
    # The minimum, not the median. Contention is one-sided: another process
    # competing for CPU or I/O can only ever make a run slower, never faster,
    # so the fastest observed run is the best estimate of what the tool costs
    # uncontended, and it is the figure least sensitive to whatever else the
    # machine happened to be doing. On this project's own hardware a median of
    # 3 produced a 5x outlier for one tool on one file and a 40x outlier on
    # another, because the interference outlasted two of the three repetitions
    # — a median cannot survive that, a minimum can. Compare across two
    # independent runs before publishing regardless.
    local wall rss
    wall=$(printf '%s\n' "${walls[@]}" | sort -n | head -1)
    rss=$(printf '%s\n' "${rsss[@]}" | sort -n | awk 'END{print $1}')

    printf '  %-24s %8ss  %6s MB\n' "$tool" "$wall" "$rss"
    echo "$tool,$label,$ftype,$(stat_of "$label" 3),$(stat_of "$label" 1),$(stat_of "$label" 2),$wall,$rss,$REPS,$nthreads" >> "$CSV"
}

# ---- Per-file benchmarks ----
run_file() {
    local f="$1" ftype="$2"
    local fname
    fname=$(basename "$f")

    echo ""
    echo "=== $fname ($ftype, $(stat_of "$fname" 1) reads, $(stat_of "$fname" 2) bp mean) ==="

    [ -n "$FASTQC" ] && bench "fastqc" 1 "$fname" "$ftype" -- \
        "$FASTQC" -t "$THREADS" -o "$RESULTSDIR/fastqc" --quiet "$f"
    [ -n "$FALCO" ] && bench "falco" 1 "$fname" "$ftype" -- \
        "$FALCO" -o "$RESULTSDIR/falco" "$f"

    bench "rastqc_t1" 1 "$fname" "$ftype" -- \
        "$RASTQC" -t 1 -o "$RESULTSDIR/rastqc" -q "$f"
    bench "rastqc_tN" "$THREADS" "$fname" "$ftype" -- \
        "$RASTQC" -t "$THREADS" -o "$RESULTSDIR/rastqc" -q "$f"

    if [ "$ftype" = "long" ]; then
        bench "rastqc_lr" "$THREADS" "$fname" "$ftype" -- \
            "$RASTQC" --long-read -t "$THREADS" -o "$RESULTSDIR/rastqc" -q "$f"
    fi
}

# A multi-file row is the case FastQC's -t actually parallelizes, and the case
# where RastQC's -t has to be split across two levels, so it is the one that
# caught the 0.1.0 memory blowup. Aggregate rows carry no per-file stats.
run_group() {
    local label="$1" ftype="$2"
    shift 2
    local files=("$@")
    [ ${#files[@]} -gt 1 ] || return 0

    local total_mb=0 total_reads=0
    for f in "${files[@]}"; do
        local b
        b=$(basename "$f")
        total_mb=$((total_mb + $(stat_of "$b" 3)))
        total_reads=$((total_reads + $(stat_of "$b" 1)))
    done
    FILE_STATS="$FILE_STATS$label $total_reads 0 $total_mb
"

    echo ""
    echo "=== $label — ${#files[@]} files together ($total_mb MB, $total_reads reads) ==="
    [ -n "$FASTQC" ] && bench "fastqc" "$THREADS" "$label" "$ftype" -- \
        "$FASTQC" -t "$THREADS" -o "$RESULTSDIR/fastqc" --quiet "${files[@]}"
    [ -n "$FALCO" ] && bench "falco" 1 "$label" "$ftype" -- \
        "$FALCO" -o "$RESULTSDIR/falco" "${files[@]}"
    bench "rastqc_tN" "$THREADS" "$label" "$ftype" -- \
        "$RASTQC" -t "$THREADS" -o "$RESULTSDIR/rastqc" -q "${files[@]}"
}

if [ ${#SHORT_FILES[@]} -gt 0 ]; then
    echo ""
    echo "━━━ SHORT-READ BENCHMARKS ━━━"
    for f in "${SHORT_FILES[@]}"; do run_file "$f" "short"; done
    run_group "ALL_SHORT" "short" "${SHORT_FILES[@]}"
fi

if [ ${#LONG_FILES[@]} -gt 0 ]; then
    echo ""
    echo "━━━ LONG-READ BENCHMARKS ━━━"
    for f in "${LONG_FILES[@]}"; do run_file "$f" "long"; done
fi

echo ""
echo "=============================================="
echo "  Results: $CSV"
echo "  Figures: python3 paper/analyze_benchmarks.py"
echo "=============================================="

# ---- Summary ----
CSV="$CSV" python3 << 'PYEOF'
import csv, os, collections

LABELS = {
    "fastqc": "FastQC",
    "falco": "Falco",
    "rastqc_t1": "RastQC -t 1",
    "rastqc_tN": "RastQC -t N",
    "rastqc_lr": "RastQC --long-read",
}

rows = collections.defaultdict(dict)
meta = {}
order = []
with open(os.environ["CSV"]) as fh:
    for row in csv.DictReader(fh):
        if row["file"] not in rows:
            order.append(row["file"])
        rows[row["file"]][row["tool"]] = {
            "wall": float(row["real_sec"]),
            "rss": row["max_rss_mb"],
        }
        meta[row["file"]] = (row["type"], row["size_mb"], row["reads"])

for ftype, heading in (("short", "SHORT-READ RESULTS"), ("long", "LONG-READ RESULTS")):
    entries = [f for f in order if meta[f][0] == ftype]
    if not entries:
        continue
    print()
    print(f"━━━ {heading} ━━━")
    for fname in entries:
        tools = rows[fname]
        _, size, reads = meta[fname]
        print(f"\n{fname}  ({size} MB, {reads} reads)")
        print(f"  {'tool':<22} {'wall':>9} {'peak RSS':>10}  {'vs FastQC':>10} {'vs Falco':>9}")
        print("  " + "─" * 66)
        for tool, m in tools.items():
            cmp = []
            for ref in ("fastqc", "falco"):
                if ref in tools and tool != ref and m["wall"] > 0:
                    cmp.append(f"{tools[ref]['wall'] / m['wall']:.1f}x")
                else:
                    cmp.append("-")
            print(f"  {LABELS.get(tool, tool):<22} {m['wall']:>8.2f}s {m['rss']:>8} MB  "
                  f"{cmp[0]:>10} {cmp[1]:>9}")
PYEOF
