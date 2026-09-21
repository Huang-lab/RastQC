# RastQC Bioconda recipe

RastQC ships through Bioconda as two packages that share one upstream source
tag but differ in build features and runtime deps:

- **`rastqc`** — core build (short-read QC). No HDF5/Arrow deps. Installs the
  `rastqc` binary.
- **`rastqc-nanopore`** — full build with `--features nanopore` (Fast5 + POD5
  readers). Depends on `hdf5` and `libarrow`. Installs the `rastqc-nanopore`
  binary (renamed from `rastqc` so it can coexist with the core package in the
  same env).

Both are **outputs of this single recipe**, because they are the same program
from the same tag and are always released together. Carrying them as two
recipes meant two version bumps and two `sha256` edits per release, and 0.2.0
showed the failure mode: two autobump PRs plus one maintainer PR for one tag.
They were merged into one recipe in
[bioconda-recipes#69452](https://github.com/bioconda/bioconda-recipes/pull/69452).

## Why the directory is `rastqc-meta`

Three Bioconda policy rules fix the layout, and they interact:

- `outputs_name_same_as_package_name` — an output may not share the top-level
  package name, so the top level cannot be called `rastqc`.
- `folder_and_package_name_must_match` — so the folder is named for the
  top-level package, `rastqc-meta`, following `openms-meta`.
- CEP 0014 — the top level builds nothing, so `requirements`, `run_exports`
  and the tests belong to the outputs, not to the top level.

`extra.skip-lints: [missing_run_exports]` is a deliberate, temporary entry.
`run_exports` is declared per output, which is what CEP 0014 asks for and what
`bioconda-utils` checks from v5.0.0 on — but the deployed linter is v4.5.0,
whose check reads the top-level `build:` section exclusively. A top-level copy
satisfies v4.5.0 and is itself an error from v5.0.0
(`unnecessary_run_exports_in_main_build_section_with_multiple_outputs`), so the
two generations cannot both be satisfied by placement. **Drop the skip once
Bioconda CI is on v5.0.0**; a v5-era linter already reads this recipe clean
with the skip inert.

## Local build + test

```bash
# 1. Install conda-build (the conda-forge rust toolchain comes with it).
conda create -n cbuild -c conda-forge conda-build
conda activate cbuild

# 2. Build both outputs in one pass.
conda build -c conda-forge -c bioconda recipes/rastqc-meta

# 3. Install both into one env and smoke-test them together — the point of the
#    binary rename is that they coexist.
conda create -n rqc-test -c local -c bioconda -c conda-forge rastqc rastqc-nanopore
conda activate rqc-test
rastqc --version
printf '@r1\nACGT\n+\n!!!!\n' | rastqc --stdin -o /tmp/rqc-out
ls /tmp/rqc-out
rastqc-nanopore --long-read --help
```

To check that a change is structural only, render it and diff the metadata
rather than trusting the build to be equivalent:

```bash
conda render -c conda-forge -c bioconda recipes/rastqc-meta
```

The `build`, `host` and `run` requirement sets, the `test` commands, the
`run_exports` pins and the computed build strings should all be unchanged for
both packages, across every `libarrow` variant. A changed build string means
already published artifacts no longer match and a rebuild is implied.

## Before submitting to bioconda-recipes

The recipe builds from the release tarball with `--locked`, so the tag must
carry a `Cargo.lock` that resolves against its own `Cargo.toml`. Check that
against the tarball rather than the working tree — a lockfile that is stale
only in the tag builds fine locally and fails in Bioconda CI:

```bash
VERSION=0.2.0

# 1. Push the release tag upstream first; the recipe source URL points at it.
git tag "v${VERSION}" && git push --tags

# 2. Set source.sha256 (one place now, not two).
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

**Diff against the live recipe before copying anything.** Bioconda's linter and
its maintainers amend merged recipes in place, and those edits are not mirrored
back here automatically:

```bash
for f in meta.yaml build-rastqc.sh build-rastqc-nanopore.sh; do
  curl -sL "https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/recipes/rastqc-meta/${f}" \
    | diff -u - "recipes/rastqc-meta/${f}" && echo "same: ${f}"
done
```

Every difference should be one you intended. Two that have bitten this project:

- `extra.additional-platforms` (`linux-aarch64`, `osx-arm64`) is added on the
  Bioconda side during review. Dropping it silently stops publishing ARM
  builds — including Apple Silicon — and nothing in Bioconda CI flags a
  platform that simply stopped being built.
- `extra.skip-lints` likewise only exists on the Bioconda side for linter
  reasons, and is expected to be removed there (see above) rather than here.

Then copy `rastqc-meta/` into `recipes/` in a fork of
[`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes),
reset `build.number` to `0` for a new version, open a PR, and address CI
feedback (linter + bulk build). Merging also rebuilds the BioContainer image.

One caveat when reading CI on a structural change: if the version and build
number already exist in the channel, `bioconda-utils` reports `not building
recipe recipes/rastqc-meta because the same number of builds are in
channel(s)` and the build jobs are a **no-op**. Green checks in that case do
not demonstrate that the recipe builds; verify locally as above.
