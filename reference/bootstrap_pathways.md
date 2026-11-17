# Bootstrap the stability and informativeness of pathways

`bootstrap_pathways()` assesses, by resampling sequences, whether each
context of a fitted tree is stable and informative: whether a sample of
the same size would reproduce it, and whether its departure from the
shorter history is reproducible. The procedure extends the edge-level
bootstrap of transition network analysis (Saqr et al., 2025) from single
transitions to variable-length pathways, and separates pathways that
reflect a regularity of the process from those that reflect only the
particular sequences that were sampled.

Each iteration resamples whole sequences with replacement (Efron &
Tibshirani, 1993) and recounts every pathway of the original tree from
raw counts, without smoothing or a minimum count. A pathway is *stable*
when the chosen statistic (its count by default, or the probability of
its most likely next state, or its divergence from its parent) stays
within `consistency_range` times the observed value in enough resamples:
`p_stability` is the share of resamples outside that band, counting the
observed sample once (\\(b + 1) / (iter + 1)\\), and the pathway is
stable when it is below `1 - stability_threshold`. A pathway is
*informative* when its likelihood-ratio statistic \\G^2\\ against its
parent context exceeds the chi-square critical value at `alpha`
(\\\|A\| - 1\\ degrees of freedom for alphabet \\A\\) in at least
`informative_threshold` of the resamples. Percentile intervals are
reported for the count, the next-state probability, the divergence and
\\G^2\\.

The two flags are read together. A stable and informative pathway is
reproducible and predicts differently from its shorter history; a stable
pathway that is not informative recurs but adds nothing to its parent's
prediction; an informative pathway that is not stable is distinctive but
carried by a subset of sequences; a pathway that is neither is fragile.
The result is a `transitiontrees_bootstrap` object, read with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) (one row
per pathway) and drawn with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

The bootstrap evaluates the pathways of one tree; to choose the depth
use
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md),
to remove contexts that do not differ from their parent use
[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md),
and to compare groups use
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md).

## Usage

``` r
bootstrap_pathways(
  tree,
  iter = 1000L,
  stat = c("count", "next_probability", "divergence"),
  consistency_range = c(0.5, 1.5),
  stability_threshold = 0.95,
  informative_threshold = 0.8,
  alpha = 0.05,
  ci_level = 0.05,
  seed = 1L,
  keep_resamples = TRUE,
  progress = FALSE
)
```

## Arguments

- tree:

  An unweighted `transitiontrees` object that carries its sequences (any
  tree from
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)).

- iter:

  Number of bootstrap resamples. Default 1000.

- stat:

  Statistic used for stability: `"count"` (default),
  `"next_probability"` (probability of the most likely next state) or
  `"divergence"` (divergence from the parent context, bits).

- consistency_range:

  Lower and upper multiples of the observed value that bound a
  consistent resample. Default `c(0.5, 1.5)`.

- stability_threshold:

  A pathway is stable when `p_stability < 1 - stability_threshold`.
  Default 0.95.

- informative_threshold:

  A pathway is informative when the share of resamples with a
  significant \\G^2\\ is at least this. Default 0.80.

- alpha:

  Significance level of the \\G^2\\ test against the parent context.
  Default 0.05.

- ci_level:

  Total tail probability of the percentile intervals. Default 0.05 (95%
  intervals).

- seed:

  Seed for the resampling, or `NULL` to use the current random-number
  state. Default 1. A non-`NULL` seed sets the session's random-number
  state.

- keep_resamples:

  Keep the per-resample values (default `TRUE`); needed by
  `as.data.frame(x, what = "resamples")` and
  [`plot_pathway_resamples()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathway_resamples.md).
  `FALSE` saves memory and leaves the summary table unchanged.

- progress:

  Show a progress bar. Default `FALSE`.

## Value

An object of class `transitiontrees_bootstrap`. `as.data.frame(x)` (or
`summary(x)`) returns one row per pathway of the original tree, stable
and informative pathways first, with columns `pathway`, `depth`,
`count`, `likely_next`, `next_probability`, `divergence` (bits),
`changes_prediction`, `G2` (`NA` for the root), `p_stability`,
`stability_rate`, `stable`, `informative_rate`, `informative`,
`flip_consistency` (share of resamples reproducing
`changes_prediction`), and for each of `count`, `next_probability`,
`divergence` and `G2` the bootstrap `mean_*`, `sd_*`, `ci_*_lo` and
`ci_*_hi`. `as.data.frame(x, what = "resamples")` returns the
per-resample values;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
\\G^2\\ intervals.

## Conditions

- `transitiontrees_bad_input`:

  (error) `tree` is not a `transitiontrees` object, or has no pathways.

- `transitiontrees_weighted_unsupported`:

  (error) the tree was fitted with weights, which resampling would drop.

- `transitiontrees_no_data`:

  (error) the tree does not carry its sequences.

## References

Saqr, M., López-Pernas, S., Törmänen, T., Kaliisa, R., Misiejuk, K., &
Tikka, S. (2025). Transition network analysis: A novel framework for
modeling, visualizing, and identifying the temporal patterns of learners
and learning processes. In *Proceedings of the 15th International
Learning Analytics and Knowledge Conference (LAK '25)* (pp. 351-361).
[doi:10.1145/3706468.3706513](https://doi.org/10.1145/3706468.3706513)

Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman & Hall/CRC.
[doi:10.1201/9780429246593](https://doi.org/10.1201/9780429246593)

## See also

[`as.data.frame.transitiontrees_bootstrap`](https://pak.dynasite.org/transitiontrees/reference/as.data.frame.transitiontrees_bootstrap.md),
[`plot.transitiontrees_bootstrap`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees_bootstrap.md)
and
[`plot_pathway_resamples()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathway_resamples.md)
for the result;
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
for the observed pathway table;
[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md)
for significance-based pruning.

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 2)
boot <- bootstrap_pathways(tree, iter = 50)
boot
#> <transitiontrees_bootstrap>  50 resamples
#>   stability  : count in [0.50, 1.50] x observed, p < 0.05
#>   informative: G^2 > qchisq(0.95, df=k-1) = 5.99, threshold 0.80
#>   pathways   : 13 total, 13 stable, 9 informative, 9 both
#> 
#> top pathways (stable + informative first):
#>                   pathway depth count p_stability stability_rate stable
#>                   Average     1   751       0.020           1.00   TRUE
#>                    Active     1   658       0.020           1.00   TRUE
#>          Active -> Active     2   433       0.020           1.00   TRUE
#>                Disengaged     1   325       0.020           1.00   TRUE
#>         Active -> Average     2   160       0.020           1.00   TRUE
#>         Average -> Active     2   144       0.020           1.00   TRUE
#>     Disengaged -> Average     2   122       0.020           1.00   TRUE
#>      Active -> Disengaged     2    23       0.020           1.00   TRUE
#>  Disengaged -> Disengaged     2   139       0.039           0.98   TRUE
#>                   (start)     0  1870       0.020           1.00   TRUE
#>  informative_rate informative mean_G2 ci_G2_lo ci_G2_hi
#>              1.00        TRUE 123.075   77.316  171.008
#>              1.00        TRUE 321.870  262.744  389.429
#>              0.98        TRUE  16.095    7.373   27.094
#>              1.00        TRUE 182.304  116.511  274.970
#>              1.00        TRUE  28.374   12.988   43.920
#>              1.00        TRUE  26.377   11.115   48.074
#>              1.00        TRUE  33.461   14.364   50.922
#>              0.86        TRUE  12.942    1.817   30.322
#>              0.98        TRUE  20.802    7.285   35.161
#>                NA       FALSE     NaN       NA       NA
#> # ... 3 more pathways (use summary(x) for full table)
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
```
