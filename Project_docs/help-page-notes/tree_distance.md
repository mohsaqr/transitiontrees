# Implementation notes removed from `?tree_distance` (2026-11-17)

Still true of R/compare.R; kept for maintainers.

- Contexts: union of `names(tree_a$nodes)` and `names(tree_b$nodes)`.
- A context absent from a tree is looked up with `.ct_match_context()` (the
  longest fitted suffix) in `.pt_lookup_prob()`; its count in that tree is 0
  (`.pt_lookup_count()`).
- Per-context divergence: `.ct_kl()` in nats; `0 log 0 := 0`; `Inf` when
  `p > 0` and `q == 0`.
- Weights: `count_a + count_b`; contexts with zero weight are masked out
  before the weighted sum; a total weight of 0 returns 0.
- Units differ across the package: `tree_distance()` / `compare_trees()` in
  nats; `tree_pathways()` divergence, `bootstrap_pathways()` divergence and
  `compare_groups()` JSD in bits.
