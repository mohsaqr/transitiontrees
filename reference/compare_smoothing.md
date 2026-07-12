# Compare Smoothing Schemes on One Dataset

Fits a context tree under several smoothing schemes — holding
`max_depth`, `nmin` and every other argument fixed — and returns a tidy
one-row-per-scheme comparison of tree size and in-sample perplexity. A
convenience wrapper over repeated
[`context_tree`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
calls that collapses the usual five-line
[`lapply()`](https://rdrr.io/r/base/lapply.html) loop into a single
call.

## Usage

``` r
compare_smoothing(
  data,
  smoothing = c("floor", "laplace", "kneser_ney", "witten_bell", "jelinek_mercer"),
  ...
)
```

## Arguments

- data:

  Either sequence data in any form accepted by
  [`context_tree`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
  (wide matrix / data.frame, list of character vectors, an `stslist`, or
  a transition/network object) — fitted afresh under each scheme —
  **or** an already-fitted `transitiontrees`, which is *re-smoothed*
  under each scheme (topology frozen, no re-count; e.g. to sweep
  smoothers on a pruned tree).

- smoothing:

  Character vector of smoothing-method names to compare. Defaults to all
  five: `"floor"`, `"laplace"`, `"kneser_ney"`, `"witten_bell"`,
  `"jelinek_mercer"`.

- ...:

  Further arguments passed to
  [`context_tree`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
  (e.g. `max_depth`, `nmin`, `alphabet`), held fixed across every
  scheme. Ignored when `data` is a fitted tree.

## Value

A `data.frame` with one row per scheme (in the order given by
`smoothing`) and columns `smoothing` (method name), `n_nodes` (tree
size) and `perplexity` (in-sample).

## Details

The perplexity reported is **in-sample** (computed on the fitting data),
so it rewards memorisation and must *not* be used to pick a smoother —
use
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
for out-of-sample selection. The point of this table is the side-by-side
view and the invariance of `n_nodes` across schemes: smoothing changes
the *probabilities* inside the tree, never *which* contexts exist
(topology is set by `nmin`, not by the smoother).

## See also

[`smooth_tree`](https://pak.dynasite.org/transitiontrees/reference/smooth_tree.md)
to re-smooth a fitted tree without re-counting;
[`tune_tree`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
for cross-validated selection.

## Examples

``` r
# \donttest{
set.seed(1)
seqs <- replicate(50, sample(c("A", "B", "C"), 12, replace = TRUE),
                  simplify = FALSE)
compare_smoothing(seqs, max_depth = 3L, min_count = 5L)
#>        smoothing n_nodes perplexity
#> 1          floor      40   2.831017
#> 2        laplace      40   2.838534
#> 3     kneser_ney      40   2.834433
#> 4    witten_bell      40   2.835931
#> 5 jelinek_mercer      40   2.867182
compare_smoothing(seqs, smoothing = c("floor", "kneser_ney"),
                  max_depth = 2L)
#>    smoothing n_nodes perplexity
#> 1      floor      13   2.938122
#> 2 kneser_ney      13   2.938189
# }
```
