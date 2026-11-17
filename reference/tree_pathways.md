# List the pathways of a fitted transition tree

`tree_pathways()` lists the pathways of a fitted tree. A pathway is a
context of the tree: the run of states that precedes a point at which
the next state is predicted, written from oldest to newest
(`"A -> B -> C"`). For each pathway the tree records how often it
occurred and the distribution of the state that followed. The pathway
table therefore shows which histories occur in the data, what tends to
follow each of them, and where a longer history changes what is expected
next.

It returns one row per pathway with its depth and count, its most likely
next state and that state's probability, the Kullback-Leibler divergence
of its next-state distribution from that of its parent (the pathway
without its oldest state), and whether the most likely next state
differs from the parent's. The divergence, in bits, measures how much
information the additional state of history carries about the next
state. The root, written `"(start)"`, conditions on no history and holds
the overall distribution of states.

The ranked views
[`common_pathways()`](https://pak.dynasite.org/transitiontrees/reference/common_pathways.md),
[`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md)
and
[`sharp_pathways()`](https://pak.dynasite.org/transitiontrees/reference/sharp_pathways.md)
return the top rows of this table by count, divergence and next-state
probability.
[`tree_dependence()`](https://pak.dynasite.org/transitiontrees/reference/tree_dependence.md)
adds the entropy of each pathway's prediction,
[`query_pathway()`](https://pak.dynasite.org/transitiontrees/reference/query_pathway.md)
looks up one pathway, and
[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md)
assesses how stable each row is under resampling of sequences.

## Usage

``` r
tree_pathways(
  tree,
  min_count = 1L,
  sort_by = c("count", "divergence", "depth"),
  decreasing = TRUE,
  ...
)
```

## Arguments

- tree:

  A `transitiontrees` object.

- min_count:

  Keep pathways that occurred at least this many times. Default 1.

- sort_by:

  Column to sort by: `"count"` (default), `"divergence"` or `"depth"`.

- decreasing:

  Sort in decreasing order. Default `TRUE`.

- ...:

  Not used.

## Value

A `data.frame` with one row per pathway and columns

- pathway:

  The context, oldest state first; `"(start)"` for the root.

- depth:

  Number of states in the pathway.

- count:

  Number of times the pathway was followed by a state (weighted for a
  weighted fit).

- likely_next:

  Most likely next state.

- next_probability:

  Its smoothed probability.

- divergence:

  Kullback-Leibler divergence, in bits, of the pathway's next-state
  distribution from its parent's; `NA` for the root, and `Inf` when the
  pathway gives probability to a state that the parent gives probability
  zero.

- changes_prediction:

  Whether `likely_next` differs from the parent's most likely next
  state; `NA` for the root.

A tree with no pathway returns a zero-row table with these columns.

## Conditions

`transitiontrees_bad_input` when `tree` is not a fitted tree.

## See also

[`common_pathways()`](https://pak.dynasite.org/transitiontrees/reference/common_pathways.md),
[`divergent_pathways()`](https://pak.dynasite.org/transitiontrees/reference/divergent_pathways.md),
[`sharp_pathways()`](https://pak.dynasite.org/transitiontrees/reference/sharp_pathways.md),
[`tree_dependence()`](https://pak.dynasite.org/transitiontrees/reference/tree_dependence.md),
[`plot_pathways()`](https://pak.dynasite.org/transitiontrees/reference/plot_pathways.md),
[`bootstrap_pathways()`](https://pak.dynasite.org/transitiontrees/reference/bootstrap_pathways.md).

## Examples

``` r
tree <- context_tree(trajectories, max_depth = 3)
tree_pathways(tree)
#>                                   pathway depth count likely_next
#> 1                                 (start)     0  1870     Average
#> 2                                 Average     1   751     Average
#> 3                                  Active     1   658      Active
#> 4                        Active -> Active     2   433      Active
#> 5                      Average -> Average     2   419     Average
#> 6                              Disengaged     1   325  Disengaged
#> 7              Active -> Active -> Active     3   316      Active
#> 8           Average -> Average -> Average     3   248     Average
#> 9                       Active -> Average     2   160     Average
#> 10                      Average -> Active     2   144      Active
#> 11               Disengaged -> Disengaged     2   139  Disengaged
#> 12                  Average -> Disengaged     2   134     Average
#> 13                  Disengaged -> Average     2   122     Average
#> 14 Disengaged -> Disengaged -> Disengaged     3    81  Disengaged
#> 15           Active -> Average -> Average     3    80     Average
#> 16            Active -> Active -> Average     3    74     Average
#> 17            Average -> Active -> Active     3    70      Active
#> 18           Average -> Average -> Active     3    66      Active
#> 19           Average -> Active -> Average     3    64      Active
#> 20       Average -> Disengaged -> Average     3    64  Disengaged
#> 21       Disengaged -> Average -> Average     3    63     Average
#> 22       Average -> Average -> Disengaged     3    62     Average
#> 23            Active -> Average -> Active     3    53     Average
#> 24    Disengaged -> Average -> Disengaged     3    46     Average
#> 25    Average -> Disengaged -> Disengaged     3    45  Disengaged
#> 26    Disengaged -> Disengaged -> Average     3    31     Average
#> 27                   Disengaged -> Active     2    29      Active
#> 28                   Active -> Disengaged     2    23      Active
#> 29        Active -> Average -> Disengaged     3    17      Active
#> 30         Disengaged -> Active -> Active     3    15      Active
#> 31        Average -> Disengaged -> Active     3    13      Active
#> 32        Disengaged -> Average -> Active     3    12     Average
#> 33        Disengaged -> Active -> Average     3    10     Average
#> 34     Disengaged -> Disengaged -> Active     3     9     Average
#> 35         Active -> Active -> Disengaged     3     9      Active
#> 36     Active -> Disengaged -> Disengaged     3     7  Disengaged
#> 37        Active -> Disengaged -> Average     3     7  Disengaged
#> 38         Active -> Disengaged -> Active     3     6      Active
#>    next_probability  divergence changes_prediction
#> 1         0.4347594          NA                 NA
#> 2         0.6098535 0.113562458              FALSE
#> 3         0.6975684 0.349487162               TRUE
#> 4         0.7852194 0.028601567              FALSE
#> 5         0.6778043 0.014667961              FALSE
#> 6         0.4830769 0.403065556               TRUE
#> 7         0.8354430 0.011491875              FALSE
#> 8         0.7177419 0.006242261              FALSE
#> 9         0.5187500 0.122825881              FALSE
#> 10        0.5000000 0.147275598              FALSE
#> 11        0.6762590 0.109773627              FALSE
#> 12        0.5000000 0.046118168               TRUE
#> 13        0.5163934 0.176186405              FALSE
#> 14        0.7037037 0.008106021              FALSE
#> 15        0.5750000 0.101876320              FALSE
#> 16        0.5540541 0.024836699              FALSE
#> 17        0.5285714 0.233006163              FALSE
#> 18        0.5757576 0.022257659              FALSE
#> 19        0.4531250 0.020767156               TRUE
#> 20        0.4843750 0.029405449               TRUE
#> 21        0.6507937 0.011442647              FALSE
#> 22        0.5161290 0.026852413              FALSE
#> 23        0.5471698 0.013380529               TRUE
#> 24        0.6086957 0.048844466              FALSE
#> 25        0.6000000 0.079870032              FALSE
#> 26        0.7085484 0.198404690              FALSE
#> 27        0.5172414 0.129904345              FALSE
#> 28        0.3913043 0.347809514               TRUE
#> 29        0.4117647 0.313295862               TRUE
#> 30        0.6666667 0.077077745              FALSE
#> 31        0.5384615 0.006062059              FALSE
#> 32        0.8318333 0.393334181               TRUE
#> 33        0.5000000 0.547770136              FALSE
#> 34        0.5555556 0.106026615               TRUE
#> 35        0.4444444 0.024571720              FALSE
#> 36        0.5714286 0.312801206              FALSE
#> 37        0.5707143 0.201375668               TRUE
#> 38        0.8318333 0.677361037              FALSE
tree_pathways(tree, min_count = 20, sort_by = "divergence")
#>                                   pathway depth count likely_next
#> 1                              Disengaged     1   325  Disengaged
#> 2                                  Active     1   658      Active
#> 3                    Active -> Disengaged     2    23      Active
#> 4             Average -> Active -> Active     3    70      Active
#> 5     Disengaged -> Disengaged -> Average     3    31     Average
#> 6                   Disengaged -> Average     2   122     Average
#> 7                       Average -> Active     2   144      Active
#> 8                    Disengaged -> Active     2    29      Active
#> 9                       Active -> Average     2   160     Average
#> 10                                Average     1   751     Average
#> 11               Disengaged -> Disengaged     2   139  Disengaged
#> 12           Active -> Average -> Average     3    80     Average
#> 13    Average -> Disengaged -> Disengaged     3    45  Disengaged
#> 14    Disengaged -> Average -> Disengaged     3    46     Average
#> 15                  Average -> Disengaged     2   134     Average
#> 16       Average -> Disengaged -> Average     3    64  Disengaged
#> 17                       Active -> Active     2   433      Active
#> 18       Average -> Average -> Disengaged     3    62     Average
#> 19            Active -> Active -> Average     3    74     Average
#> 20           Average -> Average -> Active     3    66      Active
#> 21           Average -> Active -> Average     3    64      Active
#> 22                     Average -> Average     2   419     Average
#> 23            Active -> Average -> Active     3    53     Average
#> 24             Active -> Active -> Active     3   316      Active
#> 25       Disengaged -> Average -> Average     3    63     Average
#> 26 Disengaged -> Disengaged -> Disengaged     3    81  Disengaged
#> 27          Average -> Average -> Average     3   248     Average
#> 28                                (start)     0  1870     Average
#>    next_probability  divergence changes_prediction
#> 1         0.4830769 0.403065556               TRUE
#> 2         0.6975684 0.349487162               TRUE
#> 3         0.3913043 0.347809514               TRUE
#> 4         0.5285714 0.233006163              FALSE
#> 5         0.7085484 0.198404690              FALSE
#> 6         0.5163934 0.176186405              FALSE
#> 7         0.5000000 0.147275598              FALSE
#> 8         0.5172414 0.129904345              FALSE
#> 9         0.5187500 0.122825881              FALSE
#> 10        0.6098535 0.113562458              FALSE
#> 11        0.6762590 0.109773627              FALSE
#> 12        0.5750000 0.101876320              FALSE
#> 13        0.6000000 0.079870032              FALSE
#> 14        0.6086957 0.048844466              FALSE
#> 15        0.5000000 0.046118168               TRUE
#> 16        0.4843750 0.029405449               TRUE
#> 17        0.7852194 0.028601567              FALSE
#> 18        0.5161290 0.026852413              FALSE
#> 19        0.5540541 0.024836699              FALSE
#> 20        0.5757576 0.022257659              FALSE
#> 21        0.4531250 0.020767156               TRUE
#> 22        0.6778043 0.014667961              FALSE
#> 23        0.5471698 0.013380529               TRUE
#> 24        0.8354430 0.011491875              FALSE
#> 25        0.6507937 0.011442647              FALSE
#> 26        0.7037037 0.008106021              FALSE
#> 27        0.7177419 0.006242261              FALSE
#> 28        0.4347594          NA                 NA
```
