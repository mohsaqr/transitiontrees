# Implementation notes removed from `?tree_dependence` (2026-11-17)

Still true of R/path_dependence.R; kept for maintainers.

- Renamed from `path_dependence()` to avoid a naming collision with a
  sibling package.
- Contexts whose parent is not in the node list are skipped (no row).
- `pathway` holds the internal context string; the root never appears.
- KL returns `Inf` when the context gives mass to a state its parent gives
  probability 0; entropy uses `0 log 0 = 0`.
- `base` validation is a named `stopifnot()` (unclassed error).
- Sorting is by `-column` (decreasing), stable for ties.
- **Removed claim**: the old page said "the mean of `n * KL` across rows
  recovers, up to constants, the chain-level mutual-information gain from
  the variable-depth model over the order-1 model". Not verified; dropped
  from the help page. (The count-weighted KL sum relates to the
  log-likelihood ratio of a context against its parent, which is what the
  G2 statistic measures, but the "mean" and "order-1" wording is not
  exact.)
- Cover & Thomas (2006) checked on Crossref (Wiley, DOI
  10.1002/047174882X; Crossref lists it as the online edition dated
  2005); kept without DOI.
