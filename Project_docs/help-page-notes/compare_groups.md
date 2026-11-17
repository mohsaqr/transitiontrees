# Implementation notes removed from `?compare_groups` and its methods (2026-11-17)

Still true of R/compare_groups.R; kept for maintainers.

- Contexts tested: union of the member trees' node names, plus the root if
  absent. Counts are recomputed from raw sequences with `.ct_count_table()`
  (the routine used by the fit), so observed and permuted values agree.
- JSD (`.cg_jsd_bits()`): rows with zero total dropped; weights = row totals /
  grand total; `H(mixture) - sum(w * H(row))`, floored at 0; 0 when fewer
  than two groups have the context.
- Usage G^2 (`.cg_usage_g2()`): K x 2 table `[count_c,g ; n_g - count_c,g]`
  where `n_g` is the number of windows of the context's own depth in group
  g (sum of all same-depth context counts), so prevalence is comparable
  across groups with different sequence lengths. Expected counts from
  `.cg_expected()` (shared with `plot_difference()` residuals).
- Root: usage is structurally 0 / p = 1; reported as `NA` and excluded from
  the usage BH family.
- Omnibus: behavioral = `sum(n_c * jsd_c)`; usage = `sum(usage_g2_c)`.
- Permutation counts use `stat_perm >= stat_obs - 1e-12`.
  `seed` non-NULL calls `set.seed()` on the global RNG.
- Block: when `block = NULL`, each member's `"block"` attribute (from
  `context_tree(block = )`) is filtered by the `.ct_traj` survivor index
  `attr(trajs, "idx")`; used only if every member has one. Stratified
  shuffle permutes label values within each block (loop over blocks).
- `modal_<g>`: `which.max` (ties to alphabet order); `NA` when the group has
  zero count at the context. `flips` ignores NA modals.
- `min_count` filter is applied after BH adjustment, so the FDR family
  includes contexts later dropped from the table.
- Table order: `jsd_bits` descending, then `usage_g2` descending.
- `distance_matrix`: K x K `tree_distance()` (symmetric KL, nats); exposed as
  `as.data.frame(x, what = "distance")`.
- plot(): divergence style uses `#D55E00` (significant) / `#bbbbbb`
  (n.s.), segment `grey75`; matrix style `scale_fill_gradient("#deebf7",
  "#08519c")`.
- summary(): an `NA` adjusted p-value is not counted as significant; the
  `alpha` check is an unclassed `stopifnot()`.
