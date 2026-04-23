# RastQC Bioconda recipes

RastQC is distributed through Bioconda as two sibling packages that share the
same upstream source tag but differ in build features and runtime deps:

- **`rastqc`** (this directory) — core build (short-read QC). No HDF5/Arrow deps.
  Installs the `rastqc` binary.
- **`rastqc-nanopore`** (`../rastqc-nanopore/`) — full build with
  `--features nanopore` (Fast5 + POD5 readers). Depends on `hdf5` and `libarrow`.
  Installs the `rastqc-nanopore` binary (renamed from `rastqc` so it can
  coexist with the core package in the same env).

We keep them as two separate recipes rather than a single multi-output recipe
because Bioconda's CI renders recipes without downloading source, and
multi-output Rust recipes with shared `host:` deps don't render cleanly in
that mode.

## Local build + test

```bash
# 1. Install conda-build (conda-forge rust toolchain is pulled in automatically).
conda create -n cbuild -c conda-forge conda-build
conda activate cbuild

# 2. Build the core package.
conda build -c conda-forge -c bioconda recipes/rastqc

# 3. Build the nanopore package.
conda build -c conda-forge -c bioconda recipes/rastqc-nanopore

# 4. Install and smoke-test the core package.
conda create -n rqc-test -c local -c bioconda -c conda-forge rastqc
conda activate rqc-test
rastqc --version
printf '@r1\nACGT\n+\n!!!!\n' | rastqc --stdin -o /tmp/rqc-out
ls /tmp/rqc-out

# 5. Install and smoke-test the nanopore variant.
conda create -n rqc-nano-test -c local -c bioconda -c conda-forge rastqc-nanopore
conda activate rqc-nano-test
rastqc-nanopore --long-read --help
```

## Before submitting to bioconda-recipes

1. Push a release tag upstream: `git tag v0.1.0 && git push --tags`.
2. Fill `source.sha256` in both `meta.yaml` files:
   ```bash
   curl -sL https://github.com/Huang-lab/RastQC/archive/refs/tags/v0.1.0.tar.gz | sha256sum
   ```
3. Confirm the maintainer GitHub handle in `extra.recipe-maintainers`.

## Submission

Fork [`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes),
copy both `rastqc/` and `rastqc-nanopore/` directories into `recipes/` in the
fork, open a PR, and address bioconda CI feedback (linter + bulk build).
