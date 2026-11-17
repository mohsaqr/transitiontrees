# Re-smooth the next-state probabilities of a fitted tree

`smooth_tree()` re-estimates the next-state probabilities of every
context of a fitted tree under a new smoothing scheme, without
refitting. Smoothing moves part of the probability toward unobserved
states or toward the distribution of the shorter context, the parent in
the tree (Chen & Goodman, 1999). The next-state distribution of a
context is estimated from the continuations observed after it; raw
proportions give probability zero to every continuation that was not
observed, which makes any new sequence containing such a continuation
impossible under the model, and they are unstable in contexts with few
occurrences.

Counts and the set of contexts stay as they are, and only the
probabilities change. Parents are smoothed before their children, so the
back-off schemes interpolate with an already smoothed parent. The five
schemes are given under Details. It returns a tree of the same class
with re-smoothed probabilities.

Smoothing changes predictions, perplexity and the `divergence` column of
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md),
but never which contexts the tree contains; that is set by `min_count`
and by pruning. To compare schemes on one dataset use
[`compare_smoothing()`](https://pak.dynasite.org/transitiontrees/reference/compare_smoothing.md)
(in-sample) or
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
(cross-validated).

## Usage

``` r
smooth_tree(tree, smoothing = "floor")
```

## Arguments

- tree:

  A `transitiontrees` object, or a `transitiontrees_group`, whose trees
  are each re-smoothed.

- smoothing:

  A method name, or a list of the method name and its hyperparameters,
  such as `list("kneser_ney", discount = 0.5)`. Methods and default
  hyperparameters: `"floor"` (default; `ymin = 0.001`,
  `rule = "interpolate"`), `"laplace"` (`alpha = 1`), `"kneser_ney"`
  (`discount = 0.75`), `"witten_bell"` (none) and `"jelinek_mercer"`
  (`lambda = 0.5`).

## Value

An object of the same class as `tree` with re-smoothed probabilities;
counts and contexts are unchanged.

## Details

For a context with counts \\n_s\\ over \\k\\ states, total \\n\\ and
parent distribution \\q_s\\ (uniform for the root):

- `"floor"`:

  With `rule = "interpolate"` (default), a distribution with an
  unobserved state becomes \\(1 - k\\y\_{min}) p_s + y\_{min}\\, so that
  each unobserved state has probability exactly `ymin`; distributions
  without zeros are left unchanged. `ymin` must be below \\1/k\\. With
  `rule = "cap"`, every probability is raised to at least `ymin` and the
  vector is renormalized.

- `"laplace"`:

  Additive smoothing: \\(n_s + \alpha) / (n + k\alpha)\\.

- `"kneser_ney"`:

  Absolute discounting with back-off to the parent: \\\max(n_s - d,
  0)/n\\ plus the removed mass times \\q_s\\. The parent's smoothed
  distribution is used for the back-off, not the continuation counts of
  the original method (Kneser & Ney, 1995).

- `"witten_bell"`:

  Interpolation with the parent, \\(1 - \lambda) n_s / n + \lambda
  q_s\\, with \\\lambda = r / (r + n)\\ and \\r\\ the number of distinct
  states observed after the context (Witten & Bell, 1991).

- `"jelinek_mercer"`:

  Interpolation with the parent at a fixed weight: \\(1 - \lambda) n_s /
  n + \lambda q_s\\.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree, the method
is unknown, a hyperparameter is outside its range, or `ymin` is not
below \\1/k\\ under `rule = "interpolate"`.

## References

Begleiter, R., El-Yaniv, R., & Yona, G. (2004). On prediction using
variable order Markov models. *Journal of Artificial Intelligence
Research*, 22, 385-421.
[doi:10.1613/jair.1491](https://doi.org/10.1613/jair.1491)

Chen, S. F., & Goodman, J. (1999). An empirical study of smoothing
techniques for language modeling. *Computer Speech and Language*, 13(4),
359-393.
[doi:10.1006/csla.1999.0128](https://doi.org/10.1006/csla.1999.0128)

Kneser, R., & Ney, H. (1995). Improved backing-off for M-gram language
modeling. In *1995 International Conference on Acoustics, Speech, and
Signal Processing* (Vol. 1, pp. 181-184).
[doi:10.1109/ICASSP.1995.479394](https://doi.org/10.1109/ICASSP.1995.479394)

Witten, I. H., & Bell, T. C. (1991). The zero-frequency problem:
Estimating the probabilities of novel events in adaptive text
compression. *IEEE Transactions on Information Theory*, 37(4),
1085-1094. [doi:10.1109/18.87000](https://doi.org/10.1109/18.87000)

## See also

[`compare_smoothing()`](https://pak.dynasite.org/transitiontrees/reference/compare_smoothing.md)
to compare schemes;
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md),
whose `smoothing` argument applies a scheme at fit time;
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
smooth_tree(tree, "kneser_ney")
#> <transitiontrees>  38 nodes, depth <= 3, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 136 sequences, 1870 observations
#>   smoothing: kneser_ney(discount=0.75)   min_count = 5
#> (start)   n=1870   -> Average (0.43)
#> |-- Active    n=658    -> Active (0.70)
#> |   |-- Active    n=433    -> Active (0.79)
#> |   |   |-- Active    n=316    -> Active (0.84)
#> |   |   |-- Average   n=70     -> Active (0.54)
#> |   |   `-- Disengaged  n=15     -> Active (0.73)
#> |   |-- Average   n=144    -> Active (0.51)
#> |   |   |-- Active    n=53     -> Average (0.55)
#> |   |   |-- Average   n=66     -> Active (0.58)
#> |   |   `-- Disengaged  n=12     -> Average (0.83)
#> |   `-- Disengaged  n=29     -> Active (0.55)
#> |       |-- Active    n=6      -> Active (0.84)
#> |       |-- Average   n=13     -> Active (0.58)
#> |       `-- Disengaged  n=9      -> Average (0.57)
#> |-- Average   n=751    -> Average (0.61)
#> |   |-- Active    n=160    -> Average (0.52)
#> |   |   |-- Active    n=74     -> Average (0.56)
#> |   |   |-- Average   n=64     -> Active (0.45)
#> |   |   `-- Disengaged  n=10     -> Average (0.54)
#> |   |-- Average   n=419    -> Average (0.68)
#> |   |   |-- Active    n=80     -> Average (0.58)
#> |   |   |-- Average   n=248    -> Average (0.72)
#> |   |   `-- Disengaged  n=63     -> Average (0.66)
#> |   `-- Disengaged  n=122    -> Average (0.52)
#> |       |-- Active    n=7      -> Disengaged (0.55)
#> |       |-- Average   n=64     -> Disengaged (0.49) 
#> ... 12 more nodes (use as.data.frame(x) or summary(x))
smooth_tree(tree, list("laplace", alpha = 0.5))
#> <transitiontrees>  38 nodes, depth <= 3, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 136 sequences, 1870 observations
#>   smoothing: laplace(alpha=0.5)   min_count = 5
#> (start)   n=1870   -> Average (0.43)
#> |-- Active    n=658    -> Active (0.70)
#> |   |-- Active    n=433    -> Active (0.78)
#> |   |   |-- Active    n=316    -> Active (0.83)
#> |   |   |-- Average   n=70     -> Active (0.52)
#> |   |   `-- Disengaged  n=15     -> Active (0.64)
#> |   |-- Average   n=144    -> Active (0.50)
#> |   |   |-- Active    n=53     -> Average (0.54)
#> |   |   |-- Average   n=66     -> Active (0.57)
#> |   |   `-- Disengaged  n=12     -> Average (0.78)
#> |   `-- Disengaged  n=29     -> Active (0.51)
#> |       |-- Active    n=6      -> Active (0.73)
#> |       |-- Average   n=13     -> Active (0.52)
#> |       `-- Disengaged  n=9      -> Average (0.52)
#> |-- Average   n=751    -> Average (0.61)
#> |   |-- Active    n=160    -> Average (0.52)
#> |   |   |-- Active    n=74     -> Average (0.55)
#> |   |   |-- Average   n=64     -> Active (0.45)
#> |   |   `-- Disengaged  n=10     -> Average (0.48)
#> |   |-- Average   n=419    -> Average (0.68)
#> |   |   |-- Active    n=80     -> Average (0.57)
#> |   |   |-- Average   n=248    -> Average (0.72)
#> |   |   `-- Disengaged  n=63     -> Average (0.64)
#> |   `-- Disengaged  n=122    -> Average (0.51)
#> |       |-- Active    n=7      -> Disengaged (0.53)
#> |       |-- Average   n=64     -> Disengaged (0.48) 
#> ... 12 more nodes (use as.data.frame(x) or summary(x))
```
