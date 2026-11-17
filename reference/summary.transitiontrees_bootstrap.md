# Summarize a pathway bootstrap

[`summary()`](https://rdrr.io/r/base/summary.html) returns the
per-pathway table of a
[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md)
result; identical to `as.data.frame(object)`.

## Usage

``` r
# S3 method for class 'transitiontrees_bootstrap'
summary(object, ...)
```

## Arguments

- object:

  A `transitiontrees_bootstrap` object.

- ...:

  Ignored.

## Value

A data frame with one row per pathway of the original tree; see
[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md)
for the columns.

## See also

[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md),
[`as.data.frame.transitiontrees_bootstrap`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees_bootstrap.md).

## Examples

``` r
boot <- bootstrap_pathways(context_tree(trajectories, max_depth = 1),
                           iter = 50)
summary(boot)
#>      pathway depth count likely_next next_probability divergence
#> 1    Average     1   751     Average        0.6098535  0.1135625
#> 2     Active     1   658      Active        0.6975684  0.3494872
#> 3 Disengaged     1   325  Disengaged        0.4830769  0.4030656
#> 4    (start)     0  1870     Average        0.4347594         NA
#>   changes_prediction       G2 p_stability stability_rate stable
#> 1              FALSE 118.2307  0.01960784              1   TRUE
#> 2               TRUE 318.7958  0.01960784              1   TRUE
#> 3               TRUE 181.5994  0.01960784              1   TRUE
#> 4                 NA       NA  0.01960784              1   TRUE
#>   informative_rate informative flip_consistency mean_count sd_count ci_count_lo
#> 1                1        TRUE             0.92     753.96 45.62445     688.675
#> 2                1        TRUE             0.92     653.68 51.90668     553.025
#> 3                1        TRUE             0.76     328.12 36.56509     267.025
#> 4               NA       FALSE               NA    1871.76 44.06515    1805.050
#>   ci_count_hi mean_next_probability sd_next_probability ci_next_probability_lo
#> 1     833.975             0.6122477          0.02418576              0.5645986
#> 2     744.525             0.6978621          0.02730191              0.6424099
#> 3     409.625             0.4924669          0.03847987              0.4413658
#> 4    1955.775             0.4369780          0.02131449              0.4022123
#>   ci_next_probability_hi mean_divergence sd_divergence ci_divergence_lo
#> 1              0.6557725       0.1184564    0.02819815       0.07340206
#> 2              0.7446661       0.3569934    0.04404941       0.28268814
#> 3              0.5690410       0.4004506    0.08032705       0.26755922
#> 4              0.4725535             NaN            NA               NA
#>   ci_divergence_hi  mean_G2    sd_G2  ci_G2_lo ci_G2_hi
#> 1        0.1635825 123.0750 26.92879  77.31592 171.0082
#> 2        0.4343801 321.8705 34.51486 262.74446 389.4290
#> 3        0.5618206 182.3036 43.55024 116.51096 274.9701
#> 4               NA      NaN       NA        NA       NA
```
