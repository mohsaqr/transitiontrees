# Coerce a context tree Comparison to a Tidy Data Frame

Uniform tidy-extract: returns the per-pathway divergence breakdown
(`object$pathways`) – the consumer-facing detail behind the scalar
`pdist` and the permutation `p_value`.

## Usage

``` r
# S3 method for class 'transitiontrees_comparison'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `transitiontrees_comparison`.

- row.names, optional:

  Ignored.

- ...:

  Ignored.

## Value

A data.frame with one row per context in either tree, sorted by
`divergence_sym` descending, and columns `pathway`, `count_a`, `count_b`
(the context's count in each tree, 0 when absent), `divergence_ab`,
`divergence_ba` (the two directed Kullback-Leibler divergences, nats)
and `divergence_sym` (their mean).

## Examples

``` r
grp <- context_tree(trajectories, group = rep(c("early", "late"), each = 69),
                    max_depth = 2L)
as.data.frame(compare_trees(grp, iter = 50L))
#>                     pathway count_a count_b divergence_ab divergence_ba
#> 1      Disengaged -> Active      21       8  0.5989051815  0.1733317287
#> 2      Active -> Disengaged      15       8  0.2391915314  0.1851861729
#> 3     Disengaged -> Average      71      51  0.0497042199  0.0666472792
#> 4                   Average     412     339  0.0273928234  0.0306543011
#> 5                   (start)     969     901  0.0281546463  0.0286219940
#> 6        Average -> Average     240     179  0.0209653873  0.0236481623
#> 7  Disengaged -> Disengaged      91      48  0.0212288392  0.0170358506
#> 8         Average -> Active      61      83  0.0197559778  0.0165026367
#> 9                    Active     290     368  0.0156234716  0.0119266850
#> 10               Disengaged     198     127  0.0107696516  0.0109754459
#> 11        Active -> Average      70      90  0.0068759640  0.0076289833
#> 12         Active -> Active     183     250  0.0033837180  0.0028464661
#> 13    Average -> Disengaged      80      54  0.0005517024  0.0005599975
#>    divergence_sym
#> 1    0.3861184551
#> 2    0.2121888522
#> 3    0.0581757495
#> 4    0.0290235622
#> 5    0.0283883201
#> 6    0.0223067748
#> 7    0.0191323449
#> 8    0.0181293072
#> 9    0.0137750783
#> 10   0.0108725487
#> 11   0.0072524737
#> 12   0.0031150921
#> 13   0.0005558499
```
