# Pathways from a Fitted Tree

Returns a tidy data.frame with one row per pathway (= context) in the
tree. The pathway is the sequence of states the tree conditions on; each
row reports the count, depth, modal next state, and how surprising the
next-state distribution is relative to a shorter history. This is the
substantive consumer-facing API of a fitted tree – a ranked list of
trajectories that the data actually supports, with a consistent
interpretive frame.

Renamed from `pathways()` to avoid a naming collision with a sibling
package.

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

  Integer. Drop pathways with fewer than this many occurrences. Default
  1.

- sort_by:

  Character. One of `"count"` (default), `"divergence"`, or `"depth"`.
  Sorts the returned data.frame.

- decreasing:

  Logical. Default `TRUE`.

- ...:

  Ignored.

## Value

A data.frame with columns `pathway` (arrow notation, e.g.
`"A -> B -> C"`; the root is reported as `"(start)"`), `depth` (history
length), `count`, `likely_next` (the most likely next state),
`next_probability` (its probability), `divergence` (Kullback-Leibler
divergence from the parent context's prediction, in bits; `NA` for the
root and when the parent context is absent, and `Inf` when the pathway
places probability on a state the parent predicts with probability 0),
and `changes_prediction` (logical, did the most likely next state change
vs the parent context?). The empty case returns a 0-row data.frame with
the same schema.

## Details

Each row is a pathway – a (possibly empty) sequence of states ending at
the point where a prediction is made. The `divergence` column quantifies
how much more information the pathway carries than its (k-1)-suffix in
bits. `changes_prediction = TRUE` marks pathways where the longer
history changes which next state is most likely.

## Examples

``` r
set.seed(1)
seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 3)
tree_pathways(tree)
#>        pathway depth count likely_next next_probability  divergence
#> 1      (start)     0   600           A        0.3666667          NA
#> 2            A     1   200           A        0.3800000 0.003825681
#> 3            B     1   181           A        0.3977901 0.013246404
#> 4            C     1   169           C        0.3491124 0.008133498
#> 5       B -> A     2    67           A        0.4776119 0.030761672
#> 6       A -> A     2    66           A        0.3787879 0.003552395
#> 7       A -> C     2    60           B        0.3500000 0.001283401
#> 8       B -> B     2    58           B        0.4310345 0.021097272
#> 9       A -> B     2    56           A        0.4285714 0.011114708
#> 10      C -> C     2    51           B        0.3725490 0.018405476
#> 11      C -> B     2    51           A        0.4313725 0.021492643
#> 12      C -> A     2    49           C        0.4081633 0.059682891
#> 13      B -> C     2    42           C        0.4285714 0.049377543
#> 14 B -> A -> A     3    31           C        0.3870968 0.023604265
#> 15 A -> B -> A     3    22           A        0.6363636 0.076529600
#> 16 B -> B -> B     3    22           A        0.4090909 0.136253810
#> 17 C -> B -> A     3    21           A        0.3809524 0.030903327
#> 18 A -> A -> B     3    20           B        0.6000000 0.154084905
#> 19 A -> B -> B     3    20           B        0.7000000 0.278634412
#> 20 B -> A -> C     3    19           C        0.4736842 0.061596157
#> 21 C -> C -> B     3    18           A        0.3888889 0.012294568
#> 22 A -> C -> B     3    18           A        0.3888889 0.030794436
#> 23 A -> A -> C     3    18           B        0.4444444 0.039629872
#> 24 A -> A -> A     3    18           A        0.5000000 0.046411356
#> 25 A -> C -> C     3    17           B        0.4117647 0.007127040
#> 26 C -> A -> C     3    17           A        0.4705882 0.111283847
#> 27 B -> B -> A     3    17           A        0.4117647 0.014045673
#> 28 C -> C -> A     3    17           B        0.4117647 0.040427740
#> 29 C -> A -> B     3    16           A        0.5625000 0.061854013
#> 30 A -> C -> A     3    15           C        0.4000000 0.001864762
#> 31 B -> C -> C     3    15           A        0.4666667 0.043190419
#> 32 B -> C -> A     3    15           C        0.5333333 0.046134174
#> 33 B -> A -> B     3    15           A        0.3333333 0.076821907
#> 34 B -> B -> C     3    14           A        0.4285714 0.064179031
#> 35 C -> B -> C     3    14           C        0.5000000 0.027630139
#> 36 C -> B -> B     3    12           B        0.4166667 0.037355561
#> 37 C -> C -> C     3    12           B        0.4166667 0.036270818
#> 38 C -> A -> A     3    10           A        0.5000000 0.054909881
#> 39 A -> B -> C     3     9           C        0.4444444 0.083684788
#> 40 B -> C -> B     3     8           A        0.6241250 0.492922140
#>    changes_prediction
#> 1                  NA
#> 2               FALSE
#> 3               FALSE
#> 4                TRUE
#> 5               FALSE
#> 6               FALSE
#> 7                TRUE
#> 8                TRUE
#> 9               FALSE
#> 10               TRUE
#> 11              FALSE
#> 12               TRUE
#> 13              FALSE
#> 14               TRUE
#> 15              FALSE
#> 16               TRUE
#> 17              FALSE
#> 18               TRUE
#> 19              FALSE
#> 20               TRUE
#> 21              FALSE
#> 22              FALSE
#> 23              FALSE
#> 24              FALSE
#> 25              FALSE
#> 26               TRUE
#> 27              FALSE
#> 28               TRUE
#> 29              FALSE
#> 30              FALSE
#> 31               TRUE
#> 32              FALSE
#> 33              FALSE
#> 34               TRUE
#> 35              FALSE
#> 36              FALSE
#> 37              FALSE
#> 38              FALSE
#> 39              FALSE
#> 40              FALSE
```
