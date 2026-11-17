# Implementation notes removed from `?compare_trees` and its methods (2026-11-17)

Still true of R/compare.R; kept for maintainers.

- A 2-group `transitiontrees_group` is unpacked as `grp[[1]]`, `grp[[2]]`
  (group order = key order of the object).
- Permutation: `set.seed(seed)` on the global RNG (not restored); each
  iteration draws `sample.int(n_total)`, the first `n_a` indices form group A.
  Pooled data are `c(tree_a$data, tree_b$data)`.
- Refits use each original tree's `max_depth`, `nmin` (stored field name),
  `smoothing`, and the shared alphabet; when a template tree is pruned,
  `.pt_reapply_pruning()` re-applies its stored `criterion`, `alpha`
  (default 0.05) and `threshold` (default 0.005). Missing pruning metadata
  is a `transitiontrees_bad_input` error.
- p-value: `(1 + sum(null >= observed)) / (iter + 1)`; no tolerance on the
  comparison (unlike `compare_groups()`, which uses `>= obs - 1e-12`).
- Object fields: `pdist`, `null_dist`, `p_value`, `pathways` (breakdown
  sorted by `divergence_sym` descending).
- plot(): histogram fill `#0072B2`, observed line `#D55E00`, linewidth 1.1.
- summary(): `null_q95` is `quantile(type = 7)`; `null_sd` is `NA` for one
  permutation.
