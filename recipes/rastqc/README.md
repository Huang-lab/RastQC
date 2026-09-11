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

Both recipes build from the release tarball with `--locked`, so the tag must
carry a `Cargo.lock` that resolves against its own `Cargo.toml`. Check that
against the tarball rather than the working tree — a lockfile that is stale
only in the tag is invisible locally:

```bash
VERSION=0.2.0

# 1. Push the release tag upstream first; the recipe source URL points at it.
git tag "v${VERSION}" && git push --tags

# 2. Set source.sha256 in both meta.yaml files (same tarball, same hash).
curl -sL "https://github.com/Huang-lab/RastQC/archive/refs/tags/v${VERSION}.tar.gz" \
  | shasum -a 256   # sha256sum on Linux

# 3. Confirm the tag builds the way bioconda builds it.
curl -sL "https://github.com/Huang-lab/RastQC/archive/refs/tags/v${VERSION}.tar.gz" | tar xz
cd "RastQC-${VERSION}"
cargo install --no-track --locked --root /tmp/rqc-prefix --path . --bin rastqc
/tmp/rqc-prefix/bin/rastqc --version   # must print ${VERSION}
```

Then confirm the maintainer GitHub handle in `extra.recipe-maintainers`.

## Submission

For the **first** submission of a package, fork
[`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes) and
copy both `rastqc/` and `rastqc-nanopore/` directories into `recipes/` in the
fork.

For a **version bump** of a package already on bioconda, diff against the live
recipe before copying, because bioconda's linter and its maintainers amend
merged recipes in place and those edits are not mirrored back here:

```bash
for pkg in rastqc rastqc-nanopore; do
  curl -sL "https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/recipes/${pkg}/meta.yaml" \
    | diff -u - "recipes/${pkg}/meta.yaml"
done
```

Every difference should be one you intended. In particular `extra.additional-platforms`
(`linux-aarch64`, `osx-arm64`) is added on the bioconda side during review;
dropping it silently stops publishing ARM builds, which is not something CI
flags. Reset `build.number` to `0` for a new version.

Open a PR and address bioconda CI feedback (linter + bulk build). Merging it
also rebuilds the BioContainer image for the new version.
