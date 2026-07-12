# Cross-Validated Hyperparameter Tuning for context trees

Runs k-fold cross-validation over a grid of fitting and pruning
hyperparameters. Returns a data.frame ranked by held-out perplexity. The
configuration with minimum perplexity is exposed via
`attr(result, "best")`.

Folds are at the sequence level (each fold holds out whole sequences,
not positions within sequences).

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

  Sequence data; format accepted by
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

- max_depth:

  Integer vector. Grid values for tree depth. Default `2:5`.

- min_count:

  Integer vector. Grid for minimum-count threshold. Default
  `c(3L, 5L, 10L)`.

- smoothing:

  Smoothing grid. A character vector of method names (e.g.
  `c("floor", "kneser_ney")`) — each method is tried with its default
  hyperparameters — or a list of explicit specs (e.g.
  `list(list("floor", ymin = 0.001), list("floor", ymin = 0.005))` for a
  hyperparameter sweep within one method).

- prune:

  Logical vector. Whether to apply G^2 pruning. Default
  `c(FALSE, TRUE)`.

- alpha:

  Numeric. Significance level for G^2 pruning when `prune = TRUE`.
  Default `0.05`.

- folds:

  Integer. Number of CV folds. Default 5.

- seed:

  Integer. RNG seed for reproducible folds. Default 1.

- actor, time, action, order, session, time_threshold:

  Long-format reshaping, forwarded to
  [`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md)
  when `action` is named (exactly as in
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)),
  so a long event log is reshaped before tuning rather than read as one
  row per sequence. Default `NULL`/900 (data already in sequence shape).

## Value

A `transitiontrees_tune` object: a data.frame with one row per grid
point and columns `max_depth`, `nmin`, `smoothing`, `prune`, `logLik`,
`n_scored`, `perplexity`, `n_nodes_avg`, and `folds_failed` (the number
of CV folds that errored for that configuration), sorted by `perplexity`
ascending. `attr(result, "best")` carries the minimum-perplexity row
**among configurations whose every fold scored**; it is `NULL` if none
did. A warning is issued when any configuration had a failed fold.

## Examples

``` r
# \donttest{
set.seed(1)
m <- matrix(sample(c("A","B","C"), 30 * 12, replace = TRUE), 30, 12)
tune_tree(m, max_depth = 1:3,
              smoothing = c("floor", "kneser_ney"),
              prune = FALSE, folds = 4)
#> <transitiontrees_tune>  18 configurations
#>  max_depth nmin                           smoothing prune    logLik n_scored
#>          1    3 floor(ymin=0.001, rule=interpolate) FALSE -400.6736      360
#>          1    5 floor(ymin=0.001, rule=interpolate) FALSE -400.6736      360
#>          1   10 floor(ymin=0.001, rule=interpolate) FALSE -400.6736      360
#>          1    3           kneser_ney(discount=0.75) FALSE -400.7109      360
#>          1    5           kneser_ney(discount=0.75) FALSE -400.7109      360
#>          1   10           kneser_ney(discount=0.75) FALSE -400.7109      360
#>          2    3 floor(ymin=0.001, rule=interpolate) FALSE -403.7448      360
#>          2    5 floor(ymin=0.001, rule=interpolate) FALSE -403.7448      360
#>          2   10 floor(ymin=0.001, rule=interpolate) FALSE -403.7448      360
#>          2    3           kneser_ney(discount=0.75) FALSE -404.4129      360
#>  perplexity n_nodes_avg folds_failed
#>    3.043421           4            0
#>    3.043421           4            0
#>    3.043421           4            0
#>    3.043737           4            0
#>    3.043737           4            0
#>    3.043737           4            0
#>    3.069496          13            0
#>    3.069496          13            0
#>    3.069496          13            0
#>    3.075198          13            0
#> 
#> best (min perplexity):
#>  max_depth nmin                           smoothing prune    logLik n_scored
#>          1    3 floor(ymin=0.001, rule=interpolate) FALSE -400.6736      360
#>  perplexity n_nodes_avg folds_failed
#>    3.043421           4            0
# }
```
