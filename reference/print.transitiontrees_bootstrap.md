# Print a context tree Bootstrap

Print a context tree Bootstrap

## Usage

``` r
# S3 method for class 'transitiontrees_bootstrap'
print(x, n = 10L, digits = 3L, ...)
```

## Arguments

- x:

  A `transitiontrees_bootstrap` object.

- n:

  Integer. Number of top pathways to print. Default 10.

- digits:

  Integer. Numeric digits for the printed table. Default 3.

- ...:

  Ignored.

## Value

`x` invisibly.

## Examples

``` r
set.seed(1)
seqs <- replicate(40, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
boot <- bootstrap_pathways(context_tree(seqs, max_depth = 1L),
                           iter = 50L)
print(boot, n = 5L)
#> <transitiontrees_bootstrap>  50 resamples
#>   stability  : count in [0.50, 1.50] x observed, p < 0.05
#>   informative: G^2 > qchisq(0.95, df=k-1) = 5.99, threshold 0.80
#>   pathways   : 4 total, 4 stable, 0 informative, 0 both
#> 
#> top pathways (stable + informative first):
#>  pathway depth count p_stability stability_rate stable informative_rate
#>  (start)     0   400        0.02              1   TRUE               NA
#>        A     1   133        0.02              1   TRUE             0.02
#>        B     1   115        0.02              1   TRUE             0.16
#>        C     1   112        0.02              1   TRUE             0.12
#>  informative mean_G2 ci_G2_lo ci_G2_hi
#>        FALSE     NaN       NA       NA
#>        FALSE   1.808    0.064    5.095
#>        FALSE   3.768    0.325   13.784
#>        FALSE   3.124    0.264    7.594
```
