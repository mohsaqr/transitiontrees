# Implementation notes removed from `?query_pathway`, `?pathway_exists`, `?subtree` and `?state_frequencies` (2026-11-17)

Still true of R/query.R; kept for maintainers.

- **Pathway normalisation** (`.pt_normalise_pathway()`): `NULL` or a
  zero-length vector, `"(start)"`, the legacy `"(root)"` and the internal
  sentinel `"<root>"` all map to the root. A single string is `trimws()`ed
  and used as is (arrow notation must use `" -> "` exactly); a vector is
  pasted with `" -> "`.
- `query_pathway(exact = FALSE)` uses `.ct_match_context()` (the same
  longest-suffix walk as `predict()`), capped at the tree's `max_depth`.
  With `exact = TRUE` a non-node returns a vector of `NA` (named by
  state), or `NA` for one `next_state`.
- `subtree()` keeps node names absolute and stores the cut point in
  `attr(x, "local_root")`, which `print()` reports as "subtree of:". The
  printed text tree is drawn from the root, so a subtree (which has no
  root) prints only its header; `summary()`/`as.data.frame()` list its
  contexts. The internal `.pt_children_of()` recursion collects
  descendants.
- `state_frequencies()` reads the root counts; `proportion` is `NA` when
  the total is 0. Sort is `order(-count, alphabet position)`.
