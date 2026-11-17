# Print a pathway bootstrap

[`print()`](https://rdrr.io/r/base/print.html) prints the number of
resamples, the stability and informativeness criteria, the number of
stable and informative pathways, and the first `n` pathways (stable and
informative first).

## Usage

``` r
# S3 method for class 'transitiontrees_bootstrap'
print(x, n = 10L, digits = 3L, ...)
```

## Arguments

- x:

  A `transitiontrees_bootstrap` object.

- n:

  Number of pathways to print. Default 10.

- digits:

  Rounding digits for the printed table. Default 3.

- ...:

  Ignored.

## Value

`x`, invisibly.

## See also

[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md),
[`as.data.frame.transitiontrees_bootstrap`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees_bootstrap.md)
for every pathway.

## Examples

``` r
boot <- bootstrap_pathways(context_tree(trajectories, max_depth = 2),
                           iter = 50)
print(boot, n = 5)
#> <transitiontrees_bootstrap>  50 resamples
#>   stability  : count in [0.50, 1.50] x observed, p < 0.05
#>   informative: G^2 > qchisq(0.95, df=k-1) = 5.99, threshold 0.80
#>   pathways   : 13 total, 13 stable, 9 informative, 9 both
#> 
#> top pathways (stable + informative first):
#>            pathway depth count p_stability stability_rate stable
#>            Average     1   751        0.02              1   TRUE
#>             Active     1   658        0.02              1   TRUE
#>   Active -> Active     2   433        0.02              1   TRUE
#>         Disengaged     1   325        0.02              1   TRUE
#>  Active -> Average     2   160        0.02              1   TRUE
#>  informative_rate informative mean_G2 ci_G2_lo ci_G2_hi
#>              1.00        TRUE 123.075   77.316  171.008
#>              1.00        TRUE 321.870  262.744  389.429
#>              0.98        TRUE  16.095    7.373   27.094
#>              1.00        TRUE 182.304  116.511  274.970
#>              1.00        TRUE  28.374   12.988   43.920
#> # ... 8 more pathways (use summary(x) for full table)
```
