# Coerce a context tree Bootstrap to a Tidy Data Frame

Uniform tidy-extract: returns the per-pathway summary table
(`object$summary`), so `as.data.frame(boot)` and `summary(boot)` are
interchangeable extractors.

## Usage

``` r
# S3 method for class 'transitiontrees_bootstrap'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `transitiontrees_bootstrap`.

- row.names, optional:

  Ignored.

- ...:

  Ignored.

## Value

A data.frame with one row per pathway; see
[`bootstrap_pathways`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md)
for the full column vocabulary.

## Examples

``` r
set.seed(1)
seqs <- replicate(40, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
boot <- bootstrap_pathways(context_tree(seqs, max_depth = 1L),
                           iter = 50L)
as.data.frame(boot)
#>   pathway depth count likely_next next_probability  divergence
#> 1 (start)     0   400           A        0.3675000          NA
#> 2       A     1   133           A        0.3759398 0.002588974
#> 3       B     1   115           B        0.3739130 0.014486242
#> 4       C     1   112           C        0.3839286 0.013756666
#>   changes_prediction        G2 p_stability stability_rate stable
#> 1                 NA        NA  0.01960784              1   TRUE
#> 2              FALSE 0.4773477  0.01960784              1   TRUE
#> 3               TRUE 2.3094525  0.01960784              1   TRUE
#> 4               TRUE 2.1359284  0.01960784              1   TRUE
#>   informative_rate informative flip_consistency mean_count  sd_count
#> 1               NA       FALSE               NA     400.00  0.000000
#> 2             0.02       FALSE             0.76     133.72  9.311042
#> 3             0.16       FALSE             0.56     115.84  8.903932
#> 4             0.12       FALSE             0.74     110.44 11.234350
#>   ci_count_lo ci_count_hi mean_next_probability sd_next_probability
#> 1     400.000     400.000             0.3734000          0.02025198
#> 2     116.125     149.000             0.3948461          0.03014323
#> 3     102.675     130.775             0.4013751          0.03362609
#> 4      88.450     127.775             0.3921896          0.02825135
#>   ci_next_probability_lo ci_next_probability_hi mean_divergence sd_divergence
#> 1              0.3425000              0.4100000             NaN            NA
#> 2              0.3490834              0.4605229     0.009813034   0.008027262
#> 3              0.3587512              0.4829695     0.023290566   0.019873529
#> 4              0.3526574              0.4508381     0.020619118   0.014650138
#>   ci_divergence_lo ci_divergence_hi  mean_G2    sd_G2   ci_G2_lo  ci_G2_hi
#> 1               NA               NA      NaN       NA         NA        NA
#> 2     0.0003512243       0.02845943 1.807901 1.452510 0.06425388  5.095159
#> 3     0.0022515494       0.08461442 3.768096 3.191071 0.32472904 13.784101
#> 4     0.0016409745       0.05128583 3.124366 2.223177 0.26427453  7.594148
```
