# ErythroWorkbench 4.1.0

First public release (2026-11-18). Companion to the ErythroWorkbench web server
(https://erythromics.cn).

## Modules

Six modules, each reproducing the Python engine's output bit for bit:

* `ery_purity()` — erythroid identity and admixture (rank-agreement score).
* `ery_transfer()` — which reference identity a query profile looks like, with
  out-of-distribution and vote-margin gates that can return no label.
* `ery_clock()` — storage days from a pre-trained projection, same-layer or
  through the cross-layer module bridge.
* `ery_pathways()` — within-sample pathway perturbation fingerprint.
* `ery_attribution()` — which condition class a profile resembles, with a
  permutation null over the query's feature ranks.
* `ery_age()` — position on the reticulocyte-to-senescent maturation axis from
  a lipid matrix.

## Parity with the engine

* 1,544 assertions across the six modules, compared at `1e-9` on unrounded
  values; rounded figures are checked separately. `R CMD check` reports 0
  failures.
* Permutation p-values cannot match across languages (numpy PCG64 vs R
  Mersenne-Twister), so every permuting function takes `permutations=`. The
  parity tests feed the permutations the engine actually used, which verifies
  the statistic even though the default RNG differs.
* R package name, function names and argument names are stable; see
  `?ErythroWorkbench` for the module index.

## Defects found by the parity work

Cross-language comparison of unrounded values surfaced three real bugs that
neither implementation's own tests could see:

* Engine: `_rank_rows()` assigned group ids in the sorted space but averaged
  ranks taken from the original positions. Ties were given wrong ranks, which
  shifted the permutation null downwards (26 of 152 reference samples produced a
  negative 95th percentile, which is impossible by construction). Fixed by
  scattering the averaged ranks back to the original positions.
* Engine: `parse_lipid()` summed the chains of a multi-chain lipid only when no
  total was given, and its ether-link regex required `e` to precede `)`. Both
  were silent: `TG(18:2_17:1_18:2)` parsed as 18:2 instead of 53:5, and
  `PC(34:1)e` was classified as a diacyl species. Fixed, while keeping the
  declared total authoritative when the record states one
  (`PC:38:5 (18:0-20:4w6)` remains 38:5 — the chain annotation is not allowed to
  rewrite the measurement).
* R: `ery_id_space()` classified `L-alanine` as a gene symbol. Fixed by
  aligning the decision order and thresholds with the engine.

## Reference bundle

* Frozen snapshot: 381 reference samples, 84 condition centroids, a 45-species
  lipid panel, and the in-house marker sets used for identity calls.
* `ery_provenance()` reports the engine version that produced the bundle, the
  export timestamp, and the md5 of every source file.
* Regenerate with `data-raw/export_refs.R` and `data-raw/make_fixtures.py`
  (both need a checkout of the web server). These scripts are excluded from the
  installed package but kept in the repository.

## Known limits (stated, not hidden)

* The RNA-fraction branch of the identity gate is not ported; a gene matrix is
  reported as `unverifiable` rather than silently passed.
* Upstream profiles are a snapshot, not a live query.
* Scores are rank agreement with an in-house reference, not cell-count purity:
  mature red cells are enucleated.
