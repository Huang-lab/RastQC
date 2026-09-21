#!/usr/bin/env bash
set -euo pipefail

# Download the public sequencing data used by every RastQC benchmark.
#
# The dataset committed under benchmark/data/ is deliberately tiny so the
# repository stays small, but 200k reads finish in well under a second — which
# is process startup, not throughput, and tells you nothing about how any of
# these tools behave on a real run. Fetch real data before drawing conclusions.
#
#   ./benchmark/fetch_data.sh short   # ~3.0 GB, 7 Illumina runs
#   ./benchmark/fetch_data.sh long    # ~0.7 GB, ONT + PacBio
#   ./benchmark/fetch_data.sh human   # 12 GB, one full 30X human WGS run
#   ./benchmark/fetch_data.sh nextseq # just the issue #12 reproducer (828 MB)
#   ./benchmark/fetch_data.sh all     # everything below (~16 GB)
#
# The `human` group is a single 12 GB file and is the only one that takes
# hours rather than minutes, both to fetch and to benchmark. `short long`
# covers every platform and three orders of magnitude of size without it.
#
# Then: ./benchmark/run_benchmark.sh
#
# ENA throttles hard per connection — a single stream off the FTP mirror runs
# around 200 KB/s, which is over five hours for the full set — but concurrent
# range requests do aggregate. Large files are therefore fetched as $PARALLEL
# byte ranges at once and reassembled, which measured ~2 MB/s at PARALLEL=12.
# Set PARALLEL=1 to fall back to a single plain stream.
#
# Every file is a public ENA read run, fetched over HTTPS from the ENA FTP
# mirror. The expected byte size of each comes from the ENA filereport API and
# is checked after download, because a truncated FASTQ does not announce
# itself: it parses fine and simply reports fewer reads, so a partial fetch
# silently becomes a wrong benchmark. Re-query sizes with:
#
#   curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=DRR045135\
#&result=read_run&fields=fastq_bytes,fastq_ftp&format=tsv"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATADIR="${DATADIR:-$SCRIPT_DIR/data}"
SET="${1:-all}"

# Concurrent range requests per file, and the size below which a file is just
# fetched in one stream (the per-request overhead is not worth splitting 20 MB).
PARALLEL="${PARALLEL:-12}"
CHUNK_MIN_BYTES="${CHUNK_MIN_BYTES:-33554432}"   # 32 MB
CHUNK_RETRIES="${CHUNK_RETRIES:-12}"

# group | filename | expected bytes | url | description
#
# The short-read set spans three orders of magnitude of file size and four
# Illumina instruments, so a speedup that only holds at one size or one read
# length shows up as an outlier rather than hiding in an average. DRR045135_1
# is the run from issue #12 and is the headline figure in RESULTS.md.
DATASETS=(
  "short|DRR609229_1.fastq.gz|23446394|https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR609/DRR609229/DRR609229_1.fastq.gz|iSeq 100, 720k reads — smallest short-read case"
  "short|DRR609229_2.fastq.gz|24564824|https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR609/DRR609229/DRR609229_2.fastq.gz|iSeq 100, 720k reads — R2 of the above"
  "short|ERR5897746_1.fastq.gz|335746460|https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR589/006/ERR5897746/ERR5897746_1.fastq.gz|HiSeq 1500, 4.3M reads — mid-size"
  "short|ERR5897746_2.fastq.gz|343437023|https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR589/006/ERR5897746/ERR5897746_2.fastq.gz|HiSeq 1500, 4.3M reads — R2 of the above"
  "short|DRR013000_1.fastq.gz|1499679049|https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR013/DRR013000/DRR013000_1.fastq.gz|Genome Analyzer IIx, 24.8M reads — largest short-read case"
  "nextseq|DRR045135_1.fastq.gz|868618666|https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR045/DRR045135/DRR045135_1.fastq.gz|NextSeq 500, 18.7M reads — the run from issue #12"
  "nextseq|DRR048760.fastq.gz|56324213|https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR048/DRR048760/DRR048760.fastq.gz|NextSeq 500, 1.2M reads — small single-end run"
  "long|DRR242198_1.fastq.gz|425771068|https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR242/DRR242198/DRR242198_1.fastq.gz|ONT MinION, 76k reads, ~5.9 kb mean"
  "long|DRR723651_subreads.fastq.gz|295266620|https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR723/DRR723651/DRR723651_subreads.fastq.gz|PacBio Revio, 42k reads, ~17.6 kb mean"
  "human|ERR3239334_1.fastq.gz|13021167175|https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR323/004/ERR3239334/ERR3239334_1.fastq.gz|NovaSeq 6000, 30X human WGS (1000 Genomes) — full-scale case, 12 GB"
)

case "$SET" in
    short)   WANT=("short") ;;
    long)    WANT=("long") ;;
    human)   WANT=("human") ;;
    nextseq) WANT=("nextseq") ;;
    all)     WANT=("short" "nextseq" "long" "human") ;;
    *)
        echo "Unknown dataset group '$SET' (expected short, long, human, nextseq or all)" >&2
        exit 1
        ;;
esac

wanted() {
    local g
    for g in "${WANT[@]}"; do [ "$g" = "$1" ] && return 0; done
    return 1
}

# fetch_file <url> <dest> <expected-bytes>
#
# Only moves <dest> into place once its size matches exactly, so an
# interrupted fetch can never be mistaken for a complete one.
fetch_file() {
    local url="$1" out="$2" bytes="$3"

    if [ "$bytes" -lt "$CHUNK_MIN_BYTES" ] || [ "$PARALLEL" -le 1 ]; then
        # -C - resumes a .part left by a previous interrupted run.
        curl -fL -C - --connect-timeout 30 --speed-limit 4096 --speed-time 60 \
             -o "$out.part" "$url" || return 1
    else
        local chunkdir="$out.chunks"
        local layout="$bytes:$PARALLEL"

        # Chunk boundaries are derived from $PARALLEL, but the partial chunks
        # already on disk were written under whatever $PARALLEL the interrupted
        # run used. Resuming with a different value would append bytes from the
        # new layout onto bytes from the old one: every chunk would still reach
        # its expected length, the total would still equal $bytes, the size
        # check below would pass, and the file would be silently corrupt in the
        # middle. Record the layout each set of chunks was written under and
        # start the file over if it no longer matches.
        if [ -f "$chunkdir/.layout" ]; then
            local had
            had=$(cat "$chunkdir/.layout")
            if [ "$had" != "$layout" ]; then
                echo "  chunk layout changed ($had -> $layout); refetching this file"
                rm -rf "$chunkdir"
            fi
        elif [ -d "$chunkdir" ]; then
            # A chunk dir from before the layout was recorded. Which layout
            # produced it cannot be recovered from the chunks: a chunk that was
            # complete under a larger $PARALLEL is *shorter* than its range
            # under a smaller one, so it is indistinguishable from a partial
            # chunk and would be resumed from the wrong offset. Stamping the
            # current layout on it would bless exactly the corruption the
            # recorded layout exists to prevent. Refuse rather than guess —
            # and rather than delete, since the bytes are still usable to
            # whoever knows how they were fetched.
            echo "  $name has partial chunks with no recorded layout." >&2
            echo "  Resuming them under PARALLEL=$PARALLEL could silently corrupt the file." >&2
            echo "  Discard them:" >&2
            echo "      rm -rf '$chunkdir'" >&2
            echo "  or, if you know they were fetched with PARALLEL=$PARALLEL, adopt them:" >&2
            echo "      printf '%s' '$layout' > '$chunkdir/.layout'" >&2
            return 1
        fi
        mkdir -p "$chunkdir"
        printf '%s' "$layout" > "$chunkdir/.layout"

        local csize=$(( (bytes + PARALLEL - 1) / PARALLEL ))

        local attempt
        for attempt in $(seq 1 "$CHUNK_RETRIES"); do
            local pids=() i
            for i in $(seq 0 $((PARALLEL - 1))); do
                local s=$((i * csize)) e=$(( i * csize + csize - 1 ))
                [ "$e" -ge "$bytes" ] && e=$((bytes - 1))
                [ "$s" -gt "$e" ] && continue
                local want=$((e - s + 1)) cf="$chunkdir/$i"

                local have=0
                [ -f "$cf" ] && have=$(wc -c <"$cf" | tr -d ' ')
                # Skip chunks a previous attempt already completed.
                [ "$have" -eq "$want" ] && continue
                # A chunk somehow longer than its range is corrupt; start over.
                if [ "$have" -gt "$want" ]; then
                    rm -f "$cf"; have=0
                fi

                # Resume mid-chunk by asking for the remainder of the range and
                # appending. ENA stalls connections under sustained load, and a
                # 1.4 GB file split twelve ways means each chunk is large enough
                # that restarting one from zero wastes real time.
                # --speed-time/--speed-limit abort a connection that has gone
                # quiet so the retry loop can replace it; without them a stalled
                # transfer hangs forever with no output.
                curl -fsSL -r "$((s + have))-$e" \
                     --connect-timeout 30 --speed-limit 4096 --speed-time 60 \
                     "$url" >> "$cf" &
                pids+=($!)
            done
            # Wait for every range fetch. A curl that failed or stalled is
            # not fatal here: each chunk is re-verified below and the retry
            # loop replaces whatever is still short. `"${pids[@]:-}"` yields a
            # single empty element when no chunk needed fetching — an unguarded
            # expansion of an empty array is an error under `set -u` on the
            # bash 3.2 macOS ships — so skip that element rather than waiting
            # on it. Written as three lines because `A && B || C` reads as
            # if-then-else and is not one (shellcheck SC2015).
            for i in "${pids[@]:-}"; do
                [ -n "$i" ] || continue
                wait "$i" || true
            done

            # Re-verify every chunk; retry only what is still short.
            local ok=1
            for i in $(seq 0 $((PARALLEL - 1))); do
                local s=$((i * csize)) e=$(( i * csize + csize - 1 ))
                [ "$e" -ge "$bytes" ] && e=$((bytes - 1))
                [ "$s" -gt "$e" ] && continue
                local want=$((e - s + 1)) cf="$chunkdir/$i"
                if [ ! -f "$cf" ] || [ "$(wc -c <"$cf" | tr -d ' ')" != "$want" ]; then
                    ok=0
                fi
            done
            [ "$ok" = "1" ] && break
            if [ "$attempt" = "$CHUNK_RETRIES" ]; then
                echo "  chunks still incomplete after $CHUNK_RETRIES attempts" >&2
                return 1
            fi
            echo "  retrying incomplete chunks (attempt $((attempt + 1))/$CHUNK_RETRIES)"
        done

        # Concatenate in index order. seq keeps this numeric rather than the
        # lexicographic order a glob would give (which puts 10 before 2).
        : > "$out.part"
        for i in $(seq 0 $((PARALLEL - 1))); do
            [ -f "$chunkdir/$i" ] && cat "$chunkdir/$i" >> "$out.part"
        done
        rm -rf "$chunkdir"
    fi

    local have
    have=$(wc -c <"$out.part" | tr -d ' ')
    if [ "$have" != "$bytes" ]; then
        echo "  SIZE MISMATCH: got $have, expected $bytes — left at $out.part" >&2
        return 1
    fi
    mv "$out.part" "$out"
}

mkdir -p "$DATADIR"
failed=0

for row in "${DATASETS[@]}"; do
    IFS='|' read -r group name bytes url desc <<<"$row"
    wanted "$group" || continue

    out="$DATADIR/$name"

    # Treat a file of the right size as done, but a file of the wrong size as
    # a previous truncated fetch to be resumed rather than trusted.
    if [ -f "$out" ]; then
        have=$(wc -c <"$out" | tr -d ' ')
        if [ "$have" = "$bytes" ]; then
            echo "ok       $name  ($desc)"
            continue
        fi
        if [ "$bytes" -lt "$CHUNK_MIN_BYTES" ] || [ "$PARALLEL" -le 1 ]; then
            echo "resuming $name — have $have bytes, expected $bytes"
            mv "$out" "$out.part"
        else
            # The chunked path resumes from $out.chunks, not from $out.part.
            echo "discarding wrong-sized $name — have $have bytes, expected $bytes"
            rm -f "$out"
        fi
    fi

    echo "fetching $name  ($desc)"
    if ! fetch_file "$url" "$out" "$bytes"; then
        echo "  FAILED to download $name" >&2
        failed=1
        continue
    fi
done

echo ""
echo "Data in $DATADIR:"
ls -lh "$DATADIR"/*.fastq.gz 2>/dev/null || true

if [ "$failed" -ne 0 ]; then
    echo ""
    echo "One or more downloads did not complete. Re-run to resume." >&2
    exit 1
fi
