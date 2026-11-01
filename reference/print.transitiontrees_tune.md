# Print a context tree Tuning Grid

Print a context tree Tuning Grid

## Usage

``` r
# S3 method for class 'transitiontrees_tune'
print(x, n = 10L, ...)
```

## Arguments

- x:

  A `transitiontrees_tune` object.

- n:

  Integer. Number of configurations to print. Default 10.

- ...:

  Ignored.

## Value

`x` invisibly.

## Examples

``` r
set.seed(1)
m  <- matrix(sample(c("A","B","C"), 30 * 12, replace = TRUE), 30, 12)
tg <- tune_tree(m, max_depth = 1:2, min_count = 5L, prune = FALSE,
                folds = 3L)
print(tg, n = 3L)
#> <transitiontrees_tune>  2 configurations
#>  max_depth nmin                           smoothing prune    logLik n_scored
#>          1    5 floor(ymin=0.001, rule=interpolate) FALSE -402.2564      360
#>          2    5 floor(ymin=0.001, rule=interpolate) FALSE -417.2043      360
#>  perplexity n_nodes_avg folds_failed
#>    3.056831           4            0
#>    3.186429          13            0
#> 
#> best (min perplexity):
#>  max_depth nmin                           smoothing prune    logLik n_scored
#>          1    5 floor(ymin=0.001, rule=interpolate) FALSE -402.2564      360
#>  perplexity n_nodes_avg folds_failed
#>    3.056831           4            0
```
