# Bootstrap Pathway Stability and Informativeness

Non-parametric sequence bootstrap for a fitted `transitiontrees`.
Methodologically built on Saqr, Tikka & López-Pernas (2025), extending
an edge-level bootstrap framework to variable-depth pathways.

The bootstrap tracks every pathway in the original tree. Each iteration
resamples whole sequences with replacement, aggregates raw counts per
depth, and reads each pathway's count vector directly from the resample.
**No smoothing, no `nmin` filter, no extra parameters** inside the loop
– the bootstrap operates on counts analogously to an edge-weight
bootstrap.

Two complementary measures are reported per pathway:

- `p_stability`:

  Bootstrap-estimated probability that the chosen `stat` (default
  `count`) falls outside `[cr[1] * observed, cr[2] * observed]`, with a
  +1 correction. This is a stability p-value: small values mean the
  pathway rarely fails the chosen reproducibility criterion under
  sequence-level resampling.

- `stability_rate`:

  Uncorrected descriptive companion: the fraction of resamples where the
  chosen `stat` lies inside the consistency band. A pathway whose count
  drops to zero in a resample fails the band test automatically.

- `informative_rate`:

  Fraction of resamples where the pathway's empirical \\G^2\\
  likelihood-ratio statistic against its parent context exceeds the
  chi-square critical value at level `alpha` (df = `|alphabet| - 1`).
  Tests *reproducibly significant divergence from the shorter-history
  baseline*.

Read `stable` and `informative` together:

- `stable && informative`: reproducible and predictively distinctive
  pathway.

- `stable && !informative`: reproducible pathway count / statistic, but
  not predictively distinctive from its parent.

- `!stable && informative`: sharp or divergent pathway carried by an
  unstable subset of sequences.

- `!stable && !informative`: weak or sample-fragile pathway.

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

  A fitted `transitiontrees` carrying `tree$data`.

- iter:

  Integer. Number of bootstrap iterations. Default `1000`.

- stat:

  Character. Pathway statistic on which `p_stability` is measured. One
  of `"count"` (default; the edge-weight analogue), `"next_probability"`
  (the most-likely next-state probability), or `"divergence"`.

- consistency_range:

  Numeric vector of length 2. Multiplicative tolerance band around the
  observed value. Default `c(0.5, 1.5)`: a resample counts as consistent
  when its statistic stays within half-to-one-and-a-half times the
  observed value.

- stability_threshold:

  Numeric in \\(0, 1)\\. Backward- compatible stability-rate threshold.
  A pathway is `stable = TRUE` when
  `p_stability < 1 - stability_threshold`. Default `0.95`: at least 95%
  of resamples must fall inside the (wide) `consistency_range` band.

- informative_threshold:

  Numeric in \\(0, 1)\\. A pathway is `informative = TRUE` when
  `informative_rate >= informative_threshold`. Default `0.80`: a pathway
  must clear the \\G^2\\ critical value in at least 80% of resamples to
  be called informative.

- alpha:

  Numeric in \\(0, 1)\\. Significance level for the \\G^2\\ test against
  parent. Default `0.05`.

- ci_level:

  Numeric in \\(0, 1)\\. Tail probability for the bootstrap CIs on
  `count`, `next_probability`, `divergence`, `G2`. Default `0.05` (95%
  CI).

- seed:

  Integer or `NULL`. RNG seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html); `NULL` uses the
  current RNG state. Default `1L`.

- keep_resamples:

  Logical. If `TRUE` (default), the per-iteration resample matrices
  `M_count`, `M_next_probability`, `M_divergence`, `M_G2`,
  `M_changes_prediction` are retained on the returned object. Set to
  `FALSE` to drop them (each is `iter x n_pathways`) when memory
  matters; the summary table is unaffected.

- progress:

  Logical. Show a progress bar. Default `FALSE`.

## Value

A `transitiontrees_bootstrap` object: a list with

- summary:

  Per-pathway data.frame (one row per pathway of the original tree; the
  root row has `NA` `G2`), sorted so that `stable & informative`
  pathways come first then by `stability_rate` descending. Columns: the
  canonical `pathway`, `depth`, `count`, `likely_next`,
  `next_probability`, `divergence`, `changes_prediction` (computed from
  raw counts, no smoothing); `G2` (likelihood-ratio statistic against
  the parent); `p_stability`, `stability_rate`, `stable`;
  `informative_rate`, `informative`; `flip_consistency` (share of
  resamples reproducing the original `changes_prediction` flag); and,
  for each of `count`, `next_probability`, `divergence` and `G2`, the
  bootstrap `mean_*`, `sd_*`, `ci_*_lo` and `ci_*_hi` (percentile
  interval at `ci_level`). Returned by
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

- pathways_orig:

  Empirical original-pathway statistics (no smoothing) as a tidy
  data.frame.

- M_count, M_next_probability, M_divergence, M_G2, M_changes_prediction:

  Raw resample matrices: `iter x n_pathways`, columns named by pathway.

- iter, stat, consistency_range, stability_threshold,
  informative_threshold, alpha_g2, ci_level, seed, g2_critical_value:

  Configuration (`alpha_g2` stores `alpha`).

## References

Saqr, M., Tikka, S., & López-Pernas, S. (2025). Transition Network
Analysis. In *Proceedings of the 15th International Learning Analytics
and Knowledge Conference (LAK '25)*.
[doi:10.1145/3706468.3706513](https://doi.org/10.1145/3706468.3706513)

Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman & Hall.

## Examples

``` r
set.seed(1)
seqs <- replicate(40, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 1L)
boot <- bootstrap_pathways(tree, iter = 50L)
summary(boot)
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
