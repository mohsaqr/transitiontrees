# Implementation notes removed from `?tree_pathways`, `?common_pathways`, `?divergent_pathways` and `?sharp_pathways` (2026-11-17)

Still true of R/pathways.R; kept for maintainers.

- `tree_pathways()` was renamed from `pathways()` to avoid a naming
  collision with a sibling package.
- The root is stored under the internal sentinel `"<root>"` and displayed
  as `"(start)"`.
- **divergence** is computed in bits (`log(base = 2)`), whereas the
  internal `.ct_kl()` used by `prune_tree(criterion = "KL")` works in nats.
  The two unit conventions are intentional. `NA` also when the parent
  context is absent from the node list.
- **Sorting**: `sort_by = "divergence"` treats `NA` as `-Inf` (the root
  sorts last when decreasing). Ties keep node-list order (`order()` is
  stable).
- `common_pathways(depth = )` filters after sorting by count, then
  `head(top)`.
- `divergent_pathways()` drops rows with `NA` divergence (root, missing
  parent); `flips_only` drops `NA` `changes_prediction`.
- `sharp_pathways()` sorts by count first, then re-sorts by
  `-next_probability` (stable), so ties in probability are broken by
  count.
