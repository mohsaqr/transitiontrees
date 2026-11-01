# 3. Ecosystem compatibility: TraMineR, tna, and Nestimate

`transitiontrees` is built to slot into the wider `Dynalytics`  
ecosystem (`tna`, `Nestimate`, `cograph`) as well as `TraMineR`. You do
**not** have to re-export or re-format your data:
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
accepts a fitted state-sequence or transition/network object
**directly**, auto-detected by its S3 class. It reads the sequences the
object already carries and fits the suffix tree on exactly those, over
the same state set.

It also speaks `TraMineR`’s missing-data convention with no dependency
on the package: in wide input, the void code `%` and the missing code
`*` (alongside `NA` and `""`) are treated as gaps, never as states. So a
matrix exported from a `TraMineR` state-sequence via
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) drops straight in.

This vignette shows the hand-offs on one shared dataset and confirms
they agree. The `tna` and `Nestimate` sections only run if those
packages are installed; the wide-data section needs no extra package.

## The shared dataset

The bundled `engagement` object is a wide `data.frame`: 1000 weekly
engagement sequences over three states, with `NA` marking dropout.

``` r

data(engagement)
class(engagement)
#> [1] "data.frame"
dim(engagement)
#> [1] 1000   25
```

`tna` ships a dataset with the same name. Because `data(engagement)` has
already put this package’s copy in the workspace, the `tna` and
`Nestimate` routes below read exactly the same object as the wide route.

## Route A – wide sequence data, directly

Hand the wide frame straight to
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

``` r

tree_wide <- context_tree(engagement, max_depth = 2L, min_count = 5L)
tree_wide
#> <transitiontrees>  13 nodes, depth <= 2, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 1000 sequences, 24555 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=24555  -> Active (0.51)
#> |-- Active    n=11894  -> Active (0.86)
#> |   |-- Active    n=9766   -> Active (0.86)
#> |   |-- Average   n=1473   -> Active (0.86)
#> |   `-- Disengaged  n=312    -> Active (0.85)
#> |-- Average   n=5063   -> Average (0.55)
#> |   |-- Active    n=1031   -> Average (0.54)
#> |   |-- Average   n=2670   -> Average (0.54)
#> |   `-- Disengaged  n=1043   -> Average (0.57)
#> `-- Disengaged  n=6598   -> Disengaged (0.79)
#>     |-- Active    n=581    -> Disengaged (0.82)
#>     |-- Average   n=719    -> Disengaged (0.73)
#>     `-- Disengaged  n=4960   -> Disengaged (0.79)
```

The banner lists the alphabet the data defines and the node count.

A `TraMineR` user reaches this same tree with
`context_tree(as.matrix(seq))` – the void/missing codes are dropped
automatically, so there is no need to declare the alphabet or strip the
void by hand.

## Route B – a `tna` transition-network object

A `tna` model carries its underlying sequences in its `$data` slot.
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
reads them and decodes through the model’s label set.

``` r

library(tna)
#> 'tna' package version 1.3.1
#> ------------------------------------------------------
#>   Tikka, S., López-Pernas, S., and Saqr, M. (2025). 
#>   tna: An R Package for Transition Network Analysis.
#>   Applied Psychological Measurement.
#>   https://doi.org/10.1177/01466216251348840
#> ------------------------------------------------------
#> Please type 'citation("tna")' for more citation information.
#> See the package website at https://sonsoles.me/tna/
#> 
#> Attaching package: 'tna'
#> The following object is masked _by_ '.GlobalEnv':
#> 
#>     engagement
#> The following objects are masked from 'package:transitiontrees':
#> 
#>     engagement, group_regulation_long
model_tna <- tna(engagement)
class(model_tna)
#> [1] "tna"

tree_tna <- context_tree(model_tna, max_depth = 2L, min_count = 5L)
tree_tna
#> <transitiontrees>  13 nodes, depth <= 2, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 1000 sequences, 24555 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=24555  -> Active (0.51)
#> |-- Active    n=11894  -> Active (0.86)
#> |   |-- Active    n=9766   -> Active (0.86)
#> |   |-- Average   n=1473   -> Active (0.86)
#> |   `-- Disengaged  n=312    -> Active (0.85)
#> |-- Average   n=5063   -> Average (0.55)
#> |   |-- Active    n=1031   -> Average (0.54)
#> |   |-- Average   n=2670   -> Average (0.54)
#> |   `-- Disengaged  n=1043   -> Average (0.57)
#> `-- Disengaged  n=6598   -> Disengaged (0.79)
#>     |-- Active    n=581    -> Disengaged (0.82)
#>     |-- Average   n=719    -> Disengaged (0.73)
#>     `-- Disengaged  n=4960   -> Disengaged (0.79)
```

## Route C – a `Nestimate` network object

[`Nestimate::build_tna()`](https://pak.dynasite.org/Nestimate/reference/build_tna.html)
(and the other `build_*()` constructors) return a `netobject` that
likewise carries the sequence frame and a `$nodes` label table. Same
hand-off:

``` r

library(Nestimate)
#> 
#> Attaching package: 'Nestimate'
#> The following objects are masked from 'package:tna':
#> 
#>     cluster_data, group_regulation_long, plot_mosaic
model_nest <- build_tna(engagement)
class(model_nest)
#> [1] "netobject"       "cograph_network"

tree_nest <- context_tree(model_nest, max_depth = 2L, min_count = 5L)
tree_nest
#> <transitiontrees>  13 nodes, depth <= 2, 3 states  [unpruned]
#>   alphabet : Active, Average, Disengaged
#>   fit on   : 1000 sequences, 24555 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=24555  -> Active (0.51)
#> |-- Active    n=11894  -> Active (0.86)
#> |   |-- Active    n=9766   -> Active (0.86)
#> |   |-- Average   n=1473   -> Active (0.86)
#> |   `-- Disengaged  n=312    -> Active (0.85)
#> |-- Average   n=5063   -> Average (0.55)
#> |   |-- Active    n=1031   -> Average (0.54)
#> |   |-- Average   n=2670   -> Average (0.54)
#> |   `-- Disengaged  n=1043   -> Average (0.57)
#> `-- Disengaged  n=6598   -> Disengaged (0.79)
#>     |-- Active    n=581    -> Disengaged (0.82)
#>     |-- Average   n=719    -> Disengaged (0.73)
#>     `-- Disengaged  n=4960   -> Disengaged (0.79)
```

## They agree

Because all three objects wrap the *same* sequences, the fitted trees
are identical – same alphabet, same pathways, same counts and
probabilities, same fit.

``` r

identical(tree_pathways(tree_wide), tree_pathways(tree_tna))
#> [1] TRUE
identical(tree_pathways(tree_tna), tree_pathways(tree_nest))
#> [1] TRUE

data.frame(
  route   = c("wide", "tna", "Nestimate"),
  n_nodes = c(n_nodes(tree_wide), n_nodes(tree_tna), n_nodes(tree_nest)),
  rbind(model_fit(tree_wide), model_fit(tree_tna), model_fit(tree_nest)))
#>       route n_nodes    logLik df  nobs      AIC      BIC perplexity
#> 1      wide      13 -16184.05 26 24555 32420.09 32630.92    1.93304
#> 2       tna      13 -16184.05 26 24555 32420.09 32630.92    1.93304
#> 3 Nestimate      13 -16184.05 26 24555 32420.09 32630.92    1.93304
```

The tree shares one symbol space with the source model, so the whole
pathway API –
[`common_pathways()`](https://pak.dynasite.org/transitiontrees/reference/common_pathways.md),
[`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md),
[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md),
the plots – composes onto an already-estimated network with no
conversion step.

## The boundary: sequences, never aggregated transitions

The hand-off works only when the object actually **carries its
sequences**. A *pure graph* projection – nodes, edges and weights with
the raw sequences nulled out (an aggregated transition network) – is
**rejected**, on purpose: the original sequences cannot be recovered
from edge weights, so fabricating them would be silently wrong. The same
invariant rejects a bare numeric transition matrix. A `tna` model built
from a transition matrix alone is exactly such a pure graph:

``` r

states <- c("Active", "Average", "Disengaged")
transitions <- matrix(c(0.6, 0.3, 0.1,
                        0.2, 0.5, 0.3,
                        0.1, 0.3, 0.6), nrow = 3, byrow = TRUE,
                      dimnames = list(states, states))
graph_only <- tna(transitions)   # weights only, no sequences
context_tree(graph_only)         # -> informative error, not a fabricated fit
#> Error:
#> ! This tna object carries no sequence data anywhere (no usable $data / $sequences / $seqdata / embedded netobject). It looks like a pure graph - an aggregated transition network - and transitiontrees fits on raw sequences, not on aggregated transitions: the original sequences cannot be recovered from edge weights (the same reason numeric transition matrices are rejected). Pass a sequence-bearing object instead, e.g. a wide-format transition/network object that still carries its raw sequences, or the original wide sequence data.frame.
```

`transitiontrees` fits on raw sequences. If you have only an aggregated
network, go back to the sequence data it was built from.

## Group objects

The grouped constructors compose too.
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
recognises a `group_tna` / `netobject_group` (a named list of per-group
models) and fits one tree per group, returning a `transitiontrees_group`
that
[`prune_tree()`](https://pak.dynasite.org/transitiontrees/reference/prune_tree.md),
[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md),
and
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
consume directly – see the *Advanced analysis* vignette for the group
workflow.
