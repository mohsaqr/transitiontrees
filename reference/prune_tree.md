# Prune a transition tree to the contexts that change the prediction

`prune_tree()` removes the contexts of a fitted tree that predict the
next state no differently from their shorter suffix, so that the long
contexts that remain are those in which the additional history changes
the prediction. This is the principle behind the learning of prediction
suffix trees (Ron, Singer & Tishby, 1996) and the context algorithm for
variable-length Markov chains (Bühlmann & Wyner, 1999). A fitted tree
holds every context that occurs often enough, and many such contexts add
parameters without adding information about what comes next.

It compares each leaf with its parent, the context one state shorter,
and removes the leaf when the comparison does not pass the chosen
criterion; leaves exposed by a removal are tested in turn until no leaf
can be removed. The root is never removed. With `"G2"` a leaf is kept
when the likelihood-ratio statistic of its counts against the parent's
distribution exceeds the chi-square critical value at `alpha` with \\k -
1\\ degrees of freedom (\\k\\ states); each leaf is tested separately,
without correction for multiple tests. With `"KL"` a leaf is kept when
its Kullback-Leibler divergence from the parent, in nats, exceeds
`threshold`. With `"AIC"` and `"BIC"` a leaf is kept when its
log-likelihood, penalized by \\2(k - 1)\\ or \\\log(n)(k - 1)\\ (\\n\\
the leaf's count), is better than that of the parent's prediction
(Akaike, 1974; Schwarz, 1978). The probabilities of the contexts that
remain are not changed. It returns a tree of the same class with fewer
contexts.

To see how strongly each criterion prunes a tree, use
[`compare_pruning()`](https://pak.dynasite.org/transitiontrees/reference/compare_pruning.md);
to choose between pruned and unpruned trees on held-out data, use
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md);
to draw the `"G2"` decision along the suffix chain of one pathway, use
[`plot_pruning()`](https://pak.dynasite.org/transitiontrees/reference/plot_pruning.md).
The `divergence` column of
[`tree_pathways()`](https://pak.dynasite.org/transitiontrees/reference/tree_pathways.md)
is in bits, whereas the `"KL"` threshold is in nats; divide a threshold
in nats by `log(2)` to read it in bits.

## Usage

``` r
prune_tree(
  tree,
  criterion = c("G2", "KL", "AIC", "BIC"),
  alpha = 0.05,
  threshold = 0.005
)
```

## Arguments

- tree:

  A `transitiontrees` object, or a `transitiontrees_group`, whose trees
  are each pruned.

- criterion:

  `"G2"` (default), `"KL"`, `"AIC"` or `"BIC"`.

- alpha:

  Significance level for `"G2"`, in (0, 1). Default 0.05.

- threshold:

  Minimum divergence in nats for `"KL"`. Default 0.005.

## Value

An object of the same class as `tree`, with fewer contexts. A pruned
tree prints `[pruned]` and the criterion in its header, and
[`n_nodes()`](https://pak.dynasite.org/transitiontrees/reference/n_nodes.md)
gives its size.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree, `alpha` is
not a single number in (0, 1), or `threshold` is not a single
non-negative number.

## References

Akaike, H. (1974). A new look at the statistical model identification.
*IEEE Transactions on Automatic Control*, 19(6), 716-723.
[doi:10.1109/TAC.1974.1100705](https://doi.org/10.1109/TAC.1974.1100705)

Bühlmann, P., & Wyner, A. J. (1999). Variable length Markov chains. *The
Annals of Statistics*, 27(2), 480-513.
[doi:10.1214/aos/1018031204](https://doi.org/10.1214/aos/1018031204)

Ron, D., Singer, Y., & Tishby, N. (1996). The power of amnesia: Learning
probabilistic automata with variable memory length. *Machine Learning*,
25(2-3), 117-149.
[doi:10.1023/A:1026490906255](https://doi.org/10.1023/A%3A1026490906255)

Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of
Statistics*, 6(2), 461-464.
[doi:10.1214/aos/1176344136](https://doi.org/10.1214/aos/1176344136)

## See also

[`compare_pruning()`](https://pak.dynasite.org/transitiontrees/reference/compare_pruning.md),
[`plot_pruning()`](https://pak.dynasite.org/transitiontrees/reference/plot_pruning.md),
[`tune_tree()`](https://pak.dynasite.org/transitiontrees/reference/tune_tree.md),
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 4)
prune_tree(tree, criterion = "G2", alpha = 0.05)
#> <transitiontrees>  25 nodes, depth <= 4, 3 states  [pruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 136 sequences, 1870 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#>   pruned by: G2   alpha = 0.05
#> (start)   n=1870   -> Average (0.43)
#> |-- Active    n=658    -> Active (0.70)
#> |   |-- Active    n=433    -> Active (0.79)
#> |   |   |-- Active    n=316    -> Active (0.84)
#> |   |   |   `-- Average   n=37     -> Active (0.59)
#> |   |   `-- Average   n=70     -> Active (0.53)
#> |   `-- Average   n=144    -> Active (0.50)
#> |       `-- Disengaged  n=12     -> Average (0.83)
#> |-- Average   n=751    -> Average (0.61)
#> |   |-- Active    n=160    -> Average (0.52)
#> |   |   `-- Disengaged  n=10     -> Average (0.50)
#> |   |-- Average   n=419    -> Average (0.68)
#> |   |   |-- Active    n=80     -> Average (0.57)
#> |   |   `-- Average   n=248    -> Average (0.72)
#> |   |       `-- Disengaged  n=37     -> Average (0.70)
#> |   `-- Disengaged  n=122    -> Average (0.52)
#> |       `-- Disengaged  n=31     -> Average (0.71)
#> `-- Disengaged  n=325    -> Disengaged (0.48)
#>     |-- Active    n=23     -> Active (0.39)
#>     |-- Average   n=134    -> Average (0.50)
#>     |   |-- Active    n=17     -> Active (0.41)
#>     |   `-- Average   n=62     -> Average (0.52)
#>     |       |-- Average   n=32     -> Average (0.75)
#>     |       `-- Disengaged  n=13     -> Disengaged (0.85)
#>     `-- Disengaged  n=139    -> Disengaged (0.68) 
n_nodes(prune_tree(tree, criterion = "BIC"))
#> [1] 36
```
