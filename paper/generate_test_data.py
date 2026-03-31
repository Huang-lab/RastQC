#!/usr/bin/env python3
"""Generate synthetic FASTQ files for benchmarking RastQC vs FastQC."""

import random
import gzip
import os
import sys

READ_LENGTH = 150
BASES = b"ACGT"

# Illumina-like quality profile: high quality at start, declining toward 3' end
MEAN_QUALS = [36] * 20 + [35] * 30 + [33] * 30 + [30] * 30 + [27] * 20 + [22] * 10 + [18] * 10

# Illumina Universal Adapter
ADAPTER = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

ADAPTER_RATE = 0.05  # 5% of reads get adapter contamination
ADAPTER_INSERT_SIZES = list(range(80, 145))  # variable insert sizes

def generate_quality(length):
    """Generate a realistic quality string."""
    quals = bytearray(length)
    for i in range(length):
        base_qual = MEAN_QUALS[i] if i < len(MEAN_QUALS) else 15
        q = max(2, min(40, base_qual + random.randint(-5, 5)))
        quals[i] = q + 33  # Phred+33
    return quals

def generate_sequence(length, rng):
    """Generate a random DNA sequence."""
    seq = bytearray(length)
    for i in range(length):
        seq[i] = BASES[rng.randint(0, 3)]
    return seq

def generate_fastq(outpath, num_reads, gz=False):
    """Generate a FASTQ file with num_reads reads."""
    rng = random.Random(42)  # reproducible

    opener = gzip.open if gz else open
    mode = 'wb'

    with opener(outpath, mode) as f:
        for i in range(num_reads):
            # Header
            tile = 1100 + (i % 100)
            x = 1000 + (i % 10000)
            y = 2000 + (i // 10000) % 10000
            header = f"@SYNTH:1:{tile}:{x}:{y} 1:N:0:ATCACG\n".encode()

            # Sequence
            seq = generate_sequence(READ_LENGTH, rng)

            # Add adapter contamination to some reads
            if rng.random() < ADAPTER_RATE:
                insert_size = rng.choice(ADAPTER_INSERT_SIZES)
                if insert_size < READ_LENGTH:
                    adapter_len = min(len(ADAPTER), READ_LENGTH - insert_size)
                    seq[insert_size:insert_size + adapter_len] = ADAPTER[:adapter_len]

            # Quality
            qual = generate_quality(READ_LENGTH)

            f.write(header)
            f.write(seq)
            f.write(b"\n+\n")
            f.write(qual)
            f.write(b"\n")

            if (i + 1) % 1000000 == 0:
                print(f"  {outpath}: {(i+1)//1000000}M reads...", file=sys.stderr)

    size_mb = os.path.getsize(outpath) / (1024 * 1024)
    print(f"  Created {outpath} ({size_mb:.1f} MB, {num_reads} reads)", file=sys.stderr)

if __name__ == "__main__":
    outdir = sys.argv[1] if len(sys.argv) > 1 else "paper/data"
    os.makedirs(outdir, exist_ok=True)

    datasets = [
        ("100K", 100_000),
        ("1M", 1_000_000),
        ("10M", 10_000_000),
        ("50M", 50_000_000),
    ]

    for label, n_reads in datasets:
        path = os.path.join(outdir, f"synthetic_{label}.fastq")
        print(f"Generating {label} reads...", file=sys.stderr)
        generate_fastq(path, n_reads)

        # Also create gzipped version for 1M dataset
        if label == "1M":
            gz_path = os.path.join(outdir, f"synthetic_{label}.fastq.gz")
            print(f"Generating {label} gzipped...", file=sys.stderr)
            generate_fastq(gz_path, n_reads, gz=True)

    print("Done.", file=sys.stderr)
