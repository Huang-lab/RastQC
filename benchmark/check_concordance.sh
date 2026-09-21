#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# RastQC vs FastQC module-call concordance
# ============================================================
#
# For every FASTQ in the data directory, runs FastQC as the reference and
# compares RastQC's per-module PASS/WARN/FAIL calls against it — and Falco's
# too, when Falco is installed, since it is the other FastQC reimplementation
# and the relevant question is whether RastQC matches FastQC at least as
# closely. This is the check behind the paper's concordance claim, and it is
# the one that matters most: RastQC is a drop-in replacement, so a run that is
# fast but calls a module differently from FastQC is a regression, not an
# optimization.
#
#   ./benchmark/fetch_data.sh all
#   ./benchmark/check_concordance.sh
#   DATADIR=/path/to/fastqs ./benchmark/check_concordance.sh
#
# Requires fastqc on PATH (or $FASTQC); Falco is optional. Exits non-zero if
# any shared module disagrees with FastQC, so it can be used as a gate.
#
# Only modules both tools emit are compared. RastQC's long-read modules have
# no FastQC counterpart, and FastQC 0.12 ships with Kmer Content disabled, so
# those are skipped rather than counted as disagreements.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

RASTQC="${RASTQC:-$REPO_ROOT/target/release/rastqc}"
FASTQC="${FASTQC:-$(command -v fastqc || true)}"
FALCO="${FALCO:-$(command -v falco || true)}"
DATADIR="${DATADIR:-$REPO_ROOT/benchmark/data}"
OUTDIR="${OUTDIR:-$REPO_ROOT/benchmark/results/concordance}"
THREADS="${THREADS:-4}"

if [ ! -x "$RASTQC" ]; then
    echo "ERROR: rastqc binary not found at $RASTQC (cargo build --release)" >&2
    exit 1
fi
if [ -z "$FASTQC" ]; then
    echo "ERROR: fastqc not found on PATH; set \$FASTQC" >&2
    exit 1
fi

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR/rastqc" "$OUTDIR/fastqc" "$OUTDIR/falco"

shopt -s nullglob
FILES=("$DATADIR"/*.fastq "$DATADIR"/*.fq "$DATADIR"/*.fastq.gz "$DATADIR"/*.fq.gz)
shopt -u nullglob
if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: no FASTQ files in $DATADIR — run ./benchmark/fetch_data.sh all" >&2
    exit 1
fi

echo "=============================================="
echo "  RastQC vs FastQC concordance"
echo "  Reference: $($FASTQC --version 2>&1 | head -1)"
echo "  RastQC:    $($RASTQC --version 2>&1 | head -1)"
echo "  Falco:     ${FALCO:+$($FALCO --version 2>&1 | head -1)}${FALCO:-(not found, skipped)}"
echo "  Files:     ${#FILES[@]}"
echo "=============================================="

TOTAL=0
AGREE=0
FALCO_TOTAL=0
FALCO_AGREE=0
DISAGREE_LOG="$OUTDIR/disagreements.txt"
: > "$DISAGREE_LOG"

for f in "${FILES[@]}"; do
    base=$(basename "$f")
    stem="${base%.gz}"; stem="${stem%.fastq}"; stem="${stem%.fq}"

    echo ""
    echo "=== $base ==="

    "$RASTQC" -t "$THREADS" -q --extract -o "$OUTDIR/rastqc" "$f" >/dev/null 2>&1
    "$FASTQC" -t "$THREADS" --quiet --extract -o "$OUTDIR/fastqc" "$f" >/dev/null 2>&1
    if [ -n "$FALCO" ]; then
        mkdir -p "$OUTDIR/falco/$stem"
        "$FALCO" -o "$OUTDIR/falco/$stem" "$f" >/dev/null 2>&1 || true
    fi

    r_sum=$(find "$OUTDIR/rastqc" -name summary.txt -path "*${stem}*" | head -1)
    f_sum=$(find "$OUTDIR/fastqc" -name summary.txt -path "*${stem}*" | head -1)
    # Falco writes summary.txt (and, for multiple inputs, a prefixed variant).
    fa_sum=$(find "$OUTDIR/falco/$stem" -name "*summary.txt" 2>/dev/null | head -1)

    if [ -z "$r_sum" ] || [ -z "$f_sum" ]; then
        echo "  SKIPPED — could not locate both summaries" >&2
        continue
    fi

    # summary.txt is "STATUS<TAB>Module name<TAB>filename" for both tools, so
    # join on the module name and compare the status column.
    while IFS=$'\t' read -r module r_status f_status; do
        TOTAL=$((TOTAL + 1))
        if [ "$r_status" = "$f_status" ]; then
            AGREE=$((AGREE + 1))
            printf '  %-34s %-5s == %-5s\n' "$module" "$r_status" "$f_status"
        else
            printf '  %-34s %-5s != %-5s   <-- DISAGREE\n' "$module" "$r_status" "$f_status"
            echo "$base: $module: rastqc=$r_status fastqc=$f_status" >> "$DISAGREE_LOG"
        fi
    done < <(
        join -t$'\t' -1 1 -2 1 \
            <(awk -F'\t' '{print $2"\t"$1}' "$r_sum" | sort -t$'\t' -k1,1) \
            <(awk -F'\t' '{print $2"\t"$1}' "$f_sum" | sort -t$'\t' -k1,1)
    )

    # Same comparison for Falco, reported separately rather than gating.
    if [ -n "$fa_sum" ]; then
        while IFS=$'\t' read -r module fa_status f_status; do
            FALCO_TOTAL=$((FALCO_TOTAL + 1))
            if [ "$fa_status" = "$f_status" ]; then
                FALCO_AGREE=$((FALCO_AGREE + 1))
            else
                echo "$base: $module: falco=$fa_status fastqc=$f_status" >> "$OUTDIR/falco_disagreements.txt"
            fi
        done < <(
            join -t$'\t' -1 1 -2 1 \
                <(awk -F'\t' '{print $2"\t"$1}' "$fa_sum" | sort -t$'\t' -k1,1) \
                <(awk -F'\t' '{print $2"\t"$1}' "$f_sum" | sort -t$'\t' -k1,1)
        )
    fi
done

DISAGREE=$((TOTAL - AGREE))
echo ""
echo "=============================================="
if [ "$TOTAL" -eq 0 ]; then
    echo "  No comparable modules found." >&2
    exit 1
fi
printf '  RastQC vs FastQC: %d/%d module calls identical (%.1f%%)\n' \
    "$AGREE" "$TOTAL" "$(awk -v a="$AGREE" -v t="$TOTAL" 'BEGIN{print 100*a/t}')"
if [ "$FALCO_TOTAL" -gt 0 ]; then
    printf '  Falco  vs FastQC: %d/%d module calls identical (%.1f%%)\n' \
        "$FALCO_AGREE" "$FALCO_TOTAL" \
        "$(awk -v a="$FALCO_AGREE" -v t="$FALCO_TOTAL" 'BEGIN{print 100*a/t}')"
fi
echo "=============================================="

if [ "$DISAGREE" -ne 0 ]; then
    echo ""
    echo "$DISAGREE disagreement(s), listed in $DISAGREE_LOG:" >&2
    cat "$DISAGREE_LOG" >&2
    exit 1
fi
