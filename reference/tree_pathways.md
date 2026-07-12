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
# \donttest{
seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
                  simplify = FALSE)
tree <- context_tree(seqs, max_depth = 3)
tree_pathways(tree)
#>        pathway depth count likely_next next_probability   divergence
#> 1      (start)     0   600           B        0.3566667           NA
#> 2            A     1   192           A        0.3489583 0.0039721730
#> 3            B     1   189           B        0.3703704 0.0005886290
#> 4            C     1   169           B        0.4023669 0.0064778128
#> 5       A -> A     2    64           B        0.3750000 0.0164432891
#> 6       C -> B     2    62           B        0.3709677 0.0001503802
#> 7       B -> B     2    60           B        0.3500000 0.0099442693
#> 8       B -> A     2    59           A        0.4745763 0.0661658747
#> 9       A -> C     2    59           B        0.4237288 0.0159206898
#> 10      A -> B     2    53           B        0.3773585 0.0040428332
#> 11      B -> C     2    50           B        0.4000000 0.0010788187
#> 12      C -> A     2    48           B        0.4166667 0.0321948465
#> 13      C -> C     2    45           B        0.3777778 0.0020431387
#> 14 B -> A -> A     3    27           B        0.4074074 0.0130513888
#> 15 A -> C -> B     3    24           B        0.4166667 0.0267393032
#> 16 A -> A -> C     3    22           A        0.4545455 0.0390314185
#> 17 A -> A -> B     3    21           A        0.5238095 0.1516225820
#> 18 C -> B -> A     3    21           A        0.4761905 0.0218642752
#> 19 C -> B -> B     3    20           B        0.4500000 0.0461640356
#> 20 B -> B -> C     3    20           B        0.4000000 0.0000000000
#> 21 A -> C -> A     3    19           B        0.4736842 0.1252438272
#> 22 B -> B -> B     3    18           A        0.4444444 0.0667826336
#> 23 A -> A -> A     3    18           B        0.3888889 0.0167113500
#> 24 B -> C -> B     3    18           C        0.3888889 0.0330749797
#> 25 A -> B -> A     3    17           A        0.4117647 0.1815244382
#> 26 C -> A -> B     3    17           C        0.4117647 0.0954280571
#> 27 B -> A -> C     3    17           A        0.3529412 0.0251772797
#> 28 B -> B -> A     3    16           C        0.5000000 0.1021618220
#> 29 A -> B -> B     3    16           C        0.4375000 0.0239915862
#> 30 C -> B -> C     3    16           A        0.4375000 0.0610866837
#> 31 C -> C -> B     3    14           A        0.5000000 0.1152159016
#> 32 B -> C -> A     3    14           B        0.5000000 0.1127253274
#> 33 B -> C -> C     3    13           B        0.4615385 0.0230290645
#> 34 A -> C -> C     3    13           A        0.3846154 0.0163039366
#> 35 C -> C -> A     3    13           A        0.4615385 0.2490203212
#> 36 C -> A -> C     3    13           B        0.5384615 0.1123652365
#> 37 C -> C -> C     3    12           B        0.4166667 0.0067516151
#> 38 A -> B -> C     3    11           B        0.4545455 0.0533933372
#> 39 C -> A -> A     3    11           A        0.3636364 0.0389841559
#> 40 B -> A -> B     3     9           C        0.4444444 0.1206430645
#>    changes_prediction
#> 1                  NA
#> 2                TRUE
#> 3               FALSE
#> 4               FALSE
#> 5                TRUE
#> 6               FALSE
#> 7               FALSE
#> 8               FALSE
#> 9               FALSE
#> 10              FALSE
#> 11              FALSE
#> 12               TRUE
#> 13              FALSE
#> 14              FALSE
#> 15              FALSE
#> 16               TRUE
#> 17               TRUE
#> 18              FALSE
#> 19              FALSE
#> 20              FALSE
#> 21              FALSE
#> 22               TRUE
#> 23              FALSE
#> 24               TRUE
#> 25              FALSE
#> 26               TRUE
#> 27               TRUE
#> 28               TRUE
#> 29               TRUE
#> 30               TRUE
#> 31               TRUE
#> 32              FALSE
#> 33              FALSE
#> 34               TRUE
#> 35               TRUE
#> 36              FALSE
#> 37              FALSE
#> 38              FALSE
#> 39               TRUE
#> 40               TRUE
# }
```
