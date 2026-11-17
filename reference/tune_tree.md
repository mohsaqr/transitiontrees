# Select context tree settings by cross-validated perplexity

`tune_tree()` selects the settings of a variable-order context tree by
cross-validated perplexity. Four settings decide how much history the
tree uses and how far it trusts sparse histories: the maximum depth, the
minimum count a context needs to be kept, the smoothing of the
next-state distributions, and whether the tree is pruned. A deeper, less
smoothed tree reproduces the sequences it was fitted to more closely,
but it does not necessarily predict new sequences better. Held-out
perplexity measures the second property: it is the exponential of the
mean negative log-likelihood per scored transition in sequences the
model has not seen, and can be read as the effective number of next
states the model is choosing between.

It fits a tree for every combination of the supplied values and scores
it by k-fold cross-validation. Folds are formed from whole sequences, so
a held-out sequence contributes no transitions to the fit, and every
configuration is scored on the same folds. The selected configuration is
the one with the lowest held-out perplexity among those that scored on
every fold; a configuration with a failed fold is kept in the grid but
is never selected. It returns the grid as a `transitiontrees_tune` data
frame with one row per configuration, lowest held-out perplexity first.

The grid compares settings by out-of-sample prediction. To describe the
fit of one tree to its own data use
[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md)
or
[`perplexity()`](https://pak.dynasite.org/transitiontrees/reference/perplexity.md);
to compare smoothing methods or pruning rules on one data set without
cross-validation use
[`compare_smoothing()`](https://pak.dynasite.org/transitiontrees/reference/compare_smoothing.md)
and
[`compare_pruning()`](https://pak.dynasite.org/transitiontrees/reference/compare_pruning.md).

## Usage

``` r
tune_tree(
  data,
  max_depth = 2L:5L,
  min_count = c(3L, 5L, 10L),
  smoothing = "floor",
  prune = c(FALSE, TRUE),
  alpha = 0.05,
  folds = 5L,
  seed = 1L,
  actor = NULL,
  time = NULL,
  action = NULL,
  order = NULL,
  session = NULL,
  time_threshold = 900
)
```

## Arguments

- data:

  Sequences in any format accepted by
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

- max_depth:

  Maximum depths to try. Default `2:5`.

- min_count:

  Minimum context counts to try. Default `c(3, 5, 10)`.

- smoothing:

  Smoothing methods to try: a character vector of method names, each
  used with its default settings (default `"floor"`), or a list of
  specifications such as
  `list(list("floor", ymin = 0.001), list("floor", ymin = 0.005))` to
  vary the settings of one method.

- prune:

  Whether to prune with the \\G^2\\ criterion; both values are tried by
  default (`c(FALSE, TRUE)`).

- alpha:

  Significance level of the \\G^2\\ pruning test. Default 0.05.

- folds:

  Number of cross-validation folds. Default 5.

- seed:

  Seed for the assignment of sequences to folds. Default 1. The call
  sets the session's random-number state.

- actor, time, action, order, session, time_threshold:

  Column names and session gap for a long event log, as in
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md);
  the log is reshaped with
  [`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md)
  before tuning. Default `NULL` (data already one row per sequence) and
  900 seconds.

## Value

An object of class `transitiontrees_tune`: a data frame with one row per
configuration, sorted by held-out perplexity (lowest first), and columns
`max_depth`, `min_count`, `smoothing`, `prune`, `logLik` (held-out
log-likelihood summed over the folds that scored; `NA` when none did),
`n_scored` (held-out transitions scored), `perplexity`, `n_nodes_avg`
(mean number of contexts of the fitted trees) and `folds_failed` (folds
that could not be fitted or scored). `as.data.frame(x, what = "best")`
returns the selected configuration,
[`summary()`](https://rdrr.io/r/base/summary.html) a one-row overview,
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) the
perplexity surface.

## Conditions

- `transitiontrees_bad_input`:

  (error) fewer sequences than `folds`; `smoothing` neither a character
  vector nor a list; `data` not in a sequence format.

- `transitiontrees_fold_failed`:

  (warning) at least one configuration had a fold that could not be
  fitted or scored; the message gives the first error.

- `transitiontrees_no_best`:

  (warning) no configuration scored on every fold, so none is selected.

## See also

[`summary.transitiontrees_tune`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_tune.md),
[`as.data.frame.transitiontrees_tune`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees_tune.md)
and
[`plot.transitiontrees_tune`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees_tune.md)
for the result;
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
and
[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md)
to fit the selected configuration;
[`perplexity()`](https://pak.dynasite.org/transitiontrees/reference/perplexity.md)
and
[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md)
for the fit of one tree.

## Examples

``` r
tuned <- tune_tree(trajectories, max_depth = 1:3, min_count = c(5, 10),
                   prune = c(FALSE, TRUE), folds = 3)
tuned
#> <transitiontrees_tune>  12 configurations
#>  max_depth min_count                           smoothing prune    logLik
#>          2         5 floor(ymin=0.001, rule=interpolate) FALSE -1592.852
#>          2        10 floor(ymin=0.001, rule=interpolate) FALSE -1592.852
#>          3        10 floor(ymin=0.001, rule=interpolate) FALSE -1597.931
#>          3        10 floor(ymin=0.001, rule=interpolate)  TRUE -1602.659
#>          2         5 floor(ymin=0.001, rule=interpolate)  TRUE -1603.322
#>          2        10 floor(ymin=0.001, rule=interpolate)  TRUE -1603.322
#>          3         5 floor(ymin=0.001, rule=interpolate)  TRUE -1608.892
#>          3         5 floor(ymin=0.001, rule=interpolate) FALSE -1610.315
#>          1         5 floor(ymin=0.001, rule=interpolate) FALSE -1659.904
#>          1        10 floor(ymin=0.001, rule=interpolate) FALSE -1659.904
#>  n_scored perplexity n_nodes_avg folds_failed
#>      1870   2.343845    13.00000            0
#>      1870   2.343845    13.00000            0
#>      1870   2.350220    30.00000            0
#>      1870   2.356169    14.00000            0
#>      1870   2.357004    10.33333            0
#>      1870   2.357004    10.33333            0
#>      1870   2.364036    14.33333            0
#>      1870   2.365835    36.00000            0
#>      1870   2.429412     4.00000            0
#>      1870   2.429412     4.00000            0
#> 
#> best (min perplexity):
#>  max_depth min_count                           smoothing prune    logLik
#>          2         5 floor(ymin=0.001, rule=interpolate) FALSE -1592.852
#>  n_scored perplexity n_nodes_avg folds_failed
#>      1870   2.343845          13            0
summary(tuned)
#> <transitiontrees_tune summary>  12 configurations, 0 with a failed fold
#> selected (minimum held-out perplexity):
#>  max_depth min_count                           smoothing prune perplexity
#>          2         5 floor(ymin=0.001, rule=interpolate) FALSE      2.344
#>  n_nodes_avg
#>           13
```
