# Convert a pathway bootstrap to a data frame

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the per-pathway summary of a
[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md)
result, or the values of one statistic in every resample.

## Usage

``` r
# S3 method for class 'transitiontrees_bootstrap'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("summary", "resamples"),
  stat = c("count", "next_probability", "divergence", "G2", "changes_prediction"),
  ...
)
```

## Arguments

- x:

  A `transitiontrees_bootstrap` object.

- row.names, optional:

  Ignored.

- what:

  `"summary"` (default) for one row per pathway, or `"resamples"` for
  one row per resample and pathway.

- stat:

  For `what = "resamples"`, the statistic: `"count"` (default),
  `"next_probability"`, `"divergence"`, `"G2"` or
  `"changes_prediction"`.

- ...:

  Ignored.

## Value

A data frame. For `what = "summary"`, one row per pathway with the
columns described in
[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md).
For `what = "resamples"`, one row per resample and pathway with columns
`resample`, `pathway`, `stat` and `value` (`NA` where the statistic is
undefined in that resample).

## Conditions

- `transitiontrees_no_data`:

  (error) `what = "resamples"` for a bootstrap run with
  `keep_resamples = FALSE`.

## See also

[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md),
[`plot_pathway_resamples()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathway_resamples.md).

## Examples

``` r
boot <- bootstrap_pathways(context_tree(trajectories, max_depth = 2),
                           iter = 50)
head(as.data.frame(boot))
#>             pathway depth count likely_next next_probability divergence
#> 1           Average     1   751     Average        0.6098535 0.11356246
#> 2            Active     1   658      Active        0.6975684 0.34948716
#> 3  Active -> Active     2   433      Active        0.7852194 0.02860157
#> 4        Disengaged     1   325  Disengaged        0.4830769 0.40306556
#> 5 Active -> Average     2   160     Average        0.5187500 0.12282588
#> 6 Average -> Active     2   144      Active        0.5000000 0.14727560
#>   changes_prediction        G2 p_stability stability_rate stable
#> 1              FALSE 118.23068  0.01960784              1   TRUE
#> 2               TRUE 318.79579  0.01960784              1   TRUE
#> 3              FALSE  17.16853  0.01960784              1   TRUE
#> 4               TRUE 181.59944  0.01960784              1   TRUE
#> 5              FALSE  27.24365  0.01960784              1   TRUE
#> 6              FALSE  29.40010  0.01960784              1   TRUE
#>   informative_rate informative flip_consistency mean_count sd_count ci_count_lo
#> 1             1.00        TRUE             0.92     753.96 45.62445     688.675
#> 2             1.00        TRUE             0.92     653.68 51.90668     553.025
#> 3             0.98        TRUE             1.00     430.62 46.53543     350.450
#> 4             1.00        TRUE             0.76     328.12 36.56509     267.025
#> 5             1.00        TRUE             0.98     157.62 14.35397     131.800
#> 6             1.00        TRUE             0.76     142.14 13.65913     116.450
#>   ci_count_hi mean_next_probability sd_next_probability ci_next_probability_lo
#> 1     833.975             0.6122477          0.02418576              0.5645986
#> 2     744.525             0.6978621          0.02730191              0.6424099
#> 3     502.875             0.7796172          0.02806291              0.7255254
#> 4     409.625             0.4924669          0.03847987              0.4413658
#> 5     183.875             0.5290255          0.04000390              0.4689122
#> 6     164.650             0.5319784          0.02680429              0.4937135
#>   ci_next_probability_hi mean_divergence sd_divergence ci_divergence_lo
#> 1              0.6557725      0.11845638   0.028198150       0.07340206
#> 2              0.7446661      0.35699344   0.044049407       0.28268814
#> 3              0.8298631      0.02708422   0.009997882       0.01285652
#> 4              0.5690410      0.40045064   0.080327048       0.26755922
#> 5              0.5886575      0.13062843   0.044186845       0.05549939
#> 6              0.5758983      0.13601462   0.057345683       0.05451544
#>   ci_divergence_hi   mean_G2     sd_G2   ci_G2_lo  ci_G2_hi
#> 1       0.16358250 123.07497 26.928788  77.315916 171.00818
#> 2       0.43438012 321.87048 34.514865 262.744459 389.42898
#> 3       0.04450157  16.09519  5.917881   7.372756  27.09372
#> 4       0.56182055 182.30360 43.550238 116.510956 274.97014
#> 5       0.21192655  28.37409  9.640051  12.987756  43.91998
#> 6       0.24576993  26.37677 10.415025  11.114705  48.07441
head(as.data.frame(boot, what = "resamples", stat = "G2"))
#>   resample pathway stat value
#> 1        1 (start)   G2    NA
#> 2        2 (start)   G2    NA
#> 3        3 (start)   G2    NA
#> 4        4 (start)   G2    NA
#> 5        5 (start)   G2    NA
#> 6        6 (start)   G2    NA
```
