# Dataset provenance

Every FASTQ used to benchmark or validate RastQC, where it came from, and what
it actually contains. No sequencing data is committed to this repository —
these directories are in `.gitignore` and are populated by
[`fetch_data.sh`](fetch_data.sh) or, for the legacy set, by hand.

Organism and platform below are from the ENA filereport API, keyed on the run
accession. Read counts and mean read lengths are **measured** from the files
themselves rather than taken from metadata, because ENA's `read_count` and
`base_count` are per *spot* for some runs and per *file* for others — dividing
one by the other gives the wrong read length for paired runs. Re-derive any
row with:

```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=DRR045135\
&result=read_run&fields=run_accession,scientific_name,instrument_model,library_strategy&format=tsv"

gzip -dc <file> | awk 'NR%4==2 {n++; b+=length($0)} END{print n" reads", b/n" bp mean"}'
```

## Canonical benchmark set — `benchmark/data/`

Fetched and byte-size-verified by `./benchmark/fetch_data.sh all`. These are
the datasets behind [`RESULTS.md`](RESULTS.md), the paper's tables, and
`check_concordance.sh`.

| File | Accession | Organism | Platform | Strategy | Reads | Mean len | Size |
|---|---|---|---|---|---|---|---|
| `DRR609229_1.fastq.gz` | DRR609229 | *Homo sapiens* | Illumina iSeq 100 | WXS | 719,828 | 83 bp | 22 MB |
| `DRR609229_2.fastq.gz` | DRR609229 | *Homo sapiens* | Illumina iSeq 100 | WXS | 719,828 | 84 bp | 23 MB |
| `ERR5897746_1.fastq.gz` | ERR5897746 | *Homo sapiens* | Illumina HiSeq 1500 | WXS | 4,252,217 | 126 bp | 320 MB |
| `ERR5897746_2.fastq.gz` | ERR5897746 | *Homo sapiens* | Illumina HiSeq 1500 | WXS | 4,252,217 | 126 bp | 328 MB |
| `DRR013000_1.fastq.gz` | DRR013000 | *Homo sapiens* | Genome Analyzer IIx | WXS | 24,778,423 | 76 bp | 1430 MB |
| `DRR045135_1.fastq.gz` | DRR045135 | *Dioscorea cayenensis* | NextSeq 500 | WGS | 18,731,118 | 72 bp | 828 MB |
| `DRR048760.fastq.gz` | DRR048760 | *Cricetulus griseus* | NextSeq 500 | WGS | 1,215,541 | 67 bp | 54 MB |
| `DRR242198_1.fastq.gz` | DRR242198 | *Escherichia coli* | ONT MinION | WGS | 75,766 | 5920 bp | 406 MB |
| `DRR723651_subreads.fastq.gz` | DRR723651 | *Escherichia coli* | PacBio Revio | WGS | 41,996 | 17,609 bp | 282 MB |

Notes on individual rows that are easy to describe wrongly:

- **DRR045135 is a yam (*Dioscorea*) WGS run, not a human one.** It is the run
  from [#12](https://github.com/Huang-lab/RastQC/issues/12) and is used because
  it is a large real NextSeq 500 run, which is the instrument and library type
  most short-read QC is run on — not because of its species. Describe it as a
  NextSeq run, not a human one.
- **DRR045135_1 reads are 72 bp, not 144 bp.** ENA reports `base_count`
  2,697,280,992 over 18,731,118 reads, which divides to 144 — but that base
  count covers both mates, so R1 alone is 72 bp. Dividing ENA's two fields is
  how the 144 bp figure in earlier versions of the README was arrived at.
- **DRR609229 reads are variable length**, mean 83/84 bp with a 151 bp maximum,
  so a single fixed read length does not describe it.

### `yeast_200k.fastq.gz` is neither yeast nor committed

A 200k-read subsample of **SRR19072702**, which is *SARS-CoV-2* amplicon
sequencing (Illumina NovaSeq 6000, Broad Institute patient sequencing) — not
*Saccharomyces cerevisiae*. The filename is wrong, and it has never been
committed to this repository in any revision despite older text in `RESULTS.md`
calling it "the committed sample data".

It is also too small to measure anything: 200k reads finish in under half a
second for every tool tested, which measures process startup rather than
throughput. Use `fetch_data.sh` before drawing any conclusion about
performance. It is kept only as a tiny smoke-test input.

## Legacy set — `paper/data/`

Assembled by hand for the 0.1.0 paper, never recorded anywhere, and superseded
by the canonical set above. Nothing in this repository reads from it. Its
contents are reconstructed here only because the 0.1.0 paper's concordance
claim rested on it and three of its files are not what their names say:

| File | Accession | Name claims | Actually is |
|---|---|---|---|
| `real_genomes/fly.fastq[.gz]` | ERR1942264 | fly | ***Mus musculus***, RNA-Seq — a second mouse dataset |
| `real_genomes/yeast.fastq[.gz]` | SRR19072702 | yeast | ***SARS-CoV-2***, amplicon |
| `real_genomes/zebrafish.fastq.gz` | — | zebrafish | **broken** — 196 bytes, not gzip; a failed download |
| `real_genomes/ecoli.fastq` | ERR022075 | E. coli | *E. coli* K-12 MG1655 ✓, but duplicated byte-for-byte at `real_ecoli.fastq` |
| `real_genomes/mouse.fastq[.gz]` | ERR3085830 | mouse | *Mus musculus* ✓ |
| `real_genomes/human.fastq[.gz]` | ERR3239334 | human | *Homo sapiens* ✓ (1000 Genomes 30X) |
| `synthetic_*.fastq` | — | synthetic | generated, `@SYNTH:...` headers |

So the "five model organisms" are three distinct ones — *E. coli*, *Mus
musculus* and *Homo sapiens* — with mouse counted twice under two names and a
coronavirus standing in for yeast. Whether RastQC's module calls match
FastQC's does not depend on the species, so the concordance *result* was never
affected; only the description of the validation set was wrong. Concordance is
now measured on the canonical set above by
[`check_concordance.sh`](check_concordance.sh), so nothing depends on these
files any more.
