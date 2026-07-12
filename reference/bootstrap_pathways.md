# Bootstrap Pathway Stability and Informativeness

Non-parametric sequence bootstrap for a fitted `transitiontrees`.
Methodologically built on Saqr, Tikka & López-Pernas (2025), extending
an edge-level bootstrap framework to variable-depth pathways.

The bootstrap tracks every pathway in the original tree. Each iteration
resamples whole sequences with replacement, aggregates raw counts per
depth, and reads each pathway's count vector directly from the resample.
**No smoothing, no `nmin` filter, no extra parameters** inside the loop
— the bootstrap operates on counts analogously to an edge-weight
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
  chi-square critical value at level `alpha_g2` (df = `|alphabet| - 1`).
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
  be called informative, relaxing the earlier 0.95 gate that suppressed
  many genuinely informative deeper pathways.

- alpha:

  Numeric in \\(0, 1)\\. Significance level for the \\G^2\\ test against
  parent. Default `0.05`.

- ci_level:

  Numeric in \\(0, 1)\\. Tail probability for the bootstrap CIs on
  `count`, `next_probability`, `divergence`, `G2`. Default `0.05` (95%
  CI).

- seed:

  Integer or `NULL`. RNG seed. Default `1L`.

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

  Per-pathway data.frame, sorted so that `stable & informative` pathways
  come first then by `stability_rate` descending.

- pathways_orig:

  Empirical original-pathway statistics (no smoothing) as a tidy
  data.frame.

- M_count, M_next_probability, M_divergence, M_G2, M_changes_prediction:

  Raw resample matrices: `iter x n_pathways`, columns named by pathway.

- iter, stat, consistency_range, stability_threshold,
  informative_threshold, alpha_g2, ci_level, seed, g2_critical_value:

  Configuration.

## References

Saqr, M., Tikka, S., & López-Pernas, S. (2025). Transition Network
Analysis. *LAK '25*, doi:10.1145/3706468.3706513.

## Examples

``` r
seqs <- replicate(40, sample(c("A", "B", "C"), 10, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 1L)
boot <- bootstrap_pathways(tree, iter = 50L)
summary(boot)
#>   pathway depth count likely_next next_probability  divergence
#> 1 (start)     0   400           B        0.3475000          NA
#> 2       B     1   125           C        0.3840000 0.011169100
#> 3       C     1   118           B        0.3813559 0.005931254
#> 4       A     1   117           B        0.3846154 0.007124792
#>   changes_prediction        G2 p_stability stability_rate stable
#> 1                 NA        NA  0.01960784              1   TRUE
#> 2               TRUE 1.9354576  0.01960784              1   TRUE
#> 3              FALSE 0.9702508  0.01960784              1   TRUE
#> 4              FALSE 1.1556159  0.01960784              1   TRUE
#>   informative_rate informative flip_consistency mean_count sd_count ci_count_lo
#> 1               NA       FALSE               NA     400.00 0.000000     400.000
#> 2             0.26       FALSE             0.56     125.90 9.754643     106.125
#> 3             0.02       FALSE             0.62     118.46 8.473753     101.125
#> 4             0.06       FALSE             0.68     115.64 8.183508     101.000
#>   ci_count_hi mean_next_probability sd_next_probability ci_next_probability_lo
#> 1     400.000             0.3613500          0.01461417              0.3386250
#> 2     141.775             0.3993568          0.03621698              0.3423915
#> 3     133.000             0.3970899          0.03796610              0.3474359
#> 4     131.325             0.3993929          0.03319253              0.3507759
#>   ci_next_probability_hi mean_divergence sd_divergence ci_divergence_lo
#> 1              0.3819375             NaN            NA               NA
#> 2              0.4625969      0.02204286    0.01898974     0.0008431866
#> 3              0.4641532      0.01206216    0.01271343     0.0002884114
#> 4              0.4580042      0.01479350    0.01402766     0.0003557582
#>   ci_divergence_hi  mean_G2    sd_G2   ci_G2_lo ci_G2_hi
#> 1               NA      NaN       NA         NA       NA
#> 2       0.05599427 3.791584 3.208980 0.15637108 9.892732
#> 3       0.03455939 1.945474 1.988654 0.04500692 5.413773
#> 4       0.04422685 2.387460 2.335918 0.05716678 7.151517
```
