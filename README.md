# ErythroWorkbench

R companion to the ErythroWorkbench web server (https://erythromics.cn). It
profiles the erythroid state of a sample **before** any modelling: how much of
the signal is erythroid, how far a unit has drifted from fresh, which disease
state its cells resemble, and which pathways are perturbed.

Every reference quantity is measured in house and ships with the package, so a
score computed here reproduces the web server bit for bit. The test suite checks
that against the Python engine's own output, not against hand-written values.

## Install

```r
# straight from GitHub
install.packages("remotes")
remotes::install_github("OuLab-BioMed/ErythroWorkbench")

# or from a checkout of this repository
install.packages("ErythroWorkbench", repos = NULL, type = "source")
```

The package has **no dependencies outside base R**; `remotes` is only used to
fetch it, never at runtime. The reference bundle lives in `inst/extdata` (about
640 KB compressed), so installation needs no network access beyond GitHub.

Verify the installation against the shipped reference bundle:

```r
library(ErythroWorkbench)
ery_provenance()      # engine version, export timestamp, md5 of every source file
```

## Versioning

The package version tracks the web server engine, not the website. The engine
commit that produced the shipped reference bundle is recorded in `DESCRIPTION`
as `Config/ErythroWorkbench/engineVersion` and reported by `ery_provenance()`.
A reference bundle exported from a different engine version is a different
package version — the parity tests would fail otherwise.

## Use

```r
library(ErythroWorkbench)

m <- read.delim("matrix.tsv", row.names = 1)   # features x samples

## one call, all modules
p <- ery_profile(m)
p

## or module by module
ery_purity(m)                       # erythroid identity and admixture
ery_transfer(m)                     # which reference identity this looks like
ery_clock(m)                        # storage days, pre-trained projection
ery_pathways(m, mode = "within")    # pathway perturbation fingerprint
ery_attribution(m)                  # which condition class these profiles resemble
ery_age(lipid_matrix)               # position on the reticulocyte-to-senescent axis
```

## What the numbers mean — and do not mean

* **The primary identity reading is a rank-agreement score, not cell-count
  purity.** Mature red cells are enucleated, so what a matrix measures is
  erythroid-derived RNA. The score is capped at the 95th percentile of in-house
  self-agreement, which is why it never reaches 100.
* **Label transfer is a projection, not a classifier.** Geometrically every
  sample lands next to some reference group, so two gates (out-of-distribution
  and vote margin) decide whether a label is allowed at all. Failing either
  gives **no label**, not a low-confidence one.
* **The storage clock is pre-trained.** Annotation days are used to validate the
  predictions, never to fit them. On gene or protein input there is no same-layer
  clock anywhere in the reference, so those go through a cross-layer module
  bridge and are reported as an extrapolation.
* **The within-sample pathway fingerprint answers only "high relative to the
  other pathways in this sample"**, not "high relative to other samples".

## Reference provenance

```r
ery_provenance()
```

returns the engine version that produced the bundle, when it was exported, and
the md5 of every source file. Regenerate it with `data-raw/export_refs.R`
(needs a checkout of the web server); regenerate the test fixtures with
`data-raw/make_fixtures.py`.

## Module status

| Module | Function | Verified against the engine |
|---|---|---|
| Erythroid identity and admixture | `ery_purity()` | yes, bit for bit |
| Reference identity transfer | `ery_transfer()` | yes, bit for bit |
| Storage clock (same-layer + module bridge) | `ery_clock()` | yes, bit for bit |
| Pathway perturbation fingerprint | `ery_pathways()` | yes, bit for bit |
| Condition-profile matching | `ery_attribution()` | yes, bit for bit |
| Erythroid maturation axis | `ery_age()` | yes, bit for bit |

All six modules reproduce the web server. "Bit for bit" is meant literally: the
parity tests assert equality at 1e-9 on unrounded values (the rounded figures
are compared too, but rounding at a .0005 boundary can differ by one unit in the
last digit for reasons that have nothing to do with the algorithm).

Three things are *not* covered by that claim, and each is stated where it
matters:

* **Permutation p-values** are the one quantity that cannot match across
  languages, because numpy's PCG64 and R's Mersenne-Twister are different
  random streams. Every function that permutes accepts a `permutations=`
  argument; the parity tests feed the permutations the engine actually used, so
  the *statistic* is verified even though the default RNG differs.
* **The identity gate's RNA-fraction branch is not ported.** On gene input the
  engine can gate on an erythroid-derived-RNA fraction from a two-endpoint
  mixture model; that model is not in this package, so a gene matrix is reported
  as `unverifiable` rather than pretending it was checked. Lipid input — the
  ageing module's actual use case — goes through the annotation/sample-name
  branch, which *is* ported and matches.
* **Upstream reference profiles are frozen.** The 84 condition centroids, the 381
  reference samples and the 45-species lipid panel are a snapshot. Re-export with
  `data-raw/export_refs.R` when the resource changes.

## Layout

```
R/                 module code, one file per module
inst/extdata/      frozen reference bundle (RDS + curated TSVs)
data-raw/          scripts that regenerate the bundle and the fixtures
tests/testthat/    cross-language parity tests against the Python engine
```

## Note on Windows paths

Two Windows-specific traps are documented here because both produce silent
wrong answers rather than errors:

* **`matrix - vector` is column-major recycling in R**, not row broadcasting.
  Every per-feature standardisation in this package uses `sweep(..., 1L, ...)`;
  writing `(P - mu) / sd` on a matrix silently standardises with the wrong
  numbers (it produced correlation coefficients above 1 during development).
* **Chinese string literals must be written as `\uXXXX` escapes.** R parses
  source bytes as native encoding, so a literal 血红蛋白 loads with
  `Encoding == "unknown"`, prints as `<U+8840>`, and — worse — silently fails to
  match in `grepl(..., perl = TRUE)`. Identity classification depends on Chinese
  keywords, so a mis-encoded literal would mis-classify samples without any
  error. `data-raw/escape_cjk.py --check R tests/testthat` enforces this.

R on Windows fails on non-ASCII **absolute** paths (`normalizePath()` raises a
file-name conversion error, and `devtools::load_all()` can segfault). This is why
the regeneration scripts use paths relative to the package root, and why the
module code never calls `normalizePath()` on user input.
