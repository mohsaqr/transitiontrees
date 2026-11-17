# Compare smoothing schemes on one dataset

`compare_smoothing()` fits a tree under each of several smoothing
schemes, with every other setting held fixed, and reports one row per
scheme with the size of the tree and its in-sample perplexity. Given an
already fitted tree, it re-smooths that tree under each scheme instead,
which keeps its contexts, for example after pruning.

The number of contexts is the same for every scheme, because smoothing
changes probabilities, not which contexts exist. In-sample perplexity
favors schemes that stay close to the observed proportions and therefore
is not a criterion for choosing one;
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
compares schemes on held-out data. See
[`smooth_tree()`](https://pak.dynasite.org/transitiontrees/reference/smooth_tree.md)
for the definitions of the schemes.

## Usage

``` r
compare_smoothing(
  data,
  smoothing = c("floor", "laplace", "kneser_ney", "witten_bell", "jelinek_mercer"),
  ...
)
```

## Arguments

- data:

  Sequences in any form accepted by
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md),
  fitted afresh under each scheme, or a fitted `transitiontrees`,
  re-smoothed under each scheme.

- smoothing:

  Names of the schemes to compare. Default all five: `"floor"`,
  `"laplace"`, `"kneser_ney"`, `"witten_bell"`, `"jelinek_mercer"`, each
  with its default hyperparameters.

- ...:

  Further arguments to
  [`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md),
  such as `max_depth` and `min_count`, held fixed across schemes. Not
  used when `data` is a fitted tree.

## Value

A `data.frame` with one row per scheme, in the order of `smoothing`, and
columns `smoothing`, `n_nodes` (number of contexts) and `perplexity`
(in-sample).

## Conditions

`transitiontrees_bad_input` when `smoothing` is not a non-empty
character vector of known methods, or when `data` cannot be fitted (see
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)).

## References

Chen, S. F., & Goodman, J. (1999). An empirical study of smoothing
techniques for language modeling. *Computer Speech and Language*, 13(4),
359-393.
[doi:10.1006/csla.1999.0128](https://doi.org/10.1006/csla.1999.0128)

## See also

[`smooth_tree()`](https://pak.dynasite.org/transitiontrees/reference/smooth_tree.md)
for the schemes;
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md)
for cross-validated selection;
[`compare_pruning()`](https://pak.dynasite.org/transitiontrees/reference/compare_pruning.md)
for the analogous comparison of pruning criteria.

## Examples

``` r
compare_smoothing(trajectories, max_depth = 3)
#>        smoothing n_nodes perplexity
#> 1          floor      38   2.244394
#> 2        laplace      38   2.250477
#> 3     kneser_ney      38   2.249175
#> 4    witten_bell      38   2.245722
#> 5 jelinek_mercer      38   2.281600
compare_smoothing(trajectories, smoothing = c("floor", "kneser_ney"),
                  max_depth = 2, min_count = 10)
#>    smoothing n_nodes perplexity
#> 1      floor      13   2.310147
#> 2 kneser_ney      13   2.310592
```
