# Compute the divergence between two context trees

`tree_distance()` computes the divergence between two context trees
fitted to different sets of sequences, a measure of how differently they
predict the next state. The distance is the Kullback-Leibler divergence
between the two trees' next-state distributions (symmetric by default),
computed for every context present in either tree and averaged with
weights equal to the context's combined count in the two trees. A
context missing from one tree is evaluated there at its longest fitted
suffix, which is the distribution that tree would use for prediction.
The value is 0 when the two trees predict identically and grows as their
predictions separate.

It returns the distance only, in nats. It has no reference distribution,
and it is positive even for two samples of one process; use
[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md)
for a permutation test and a per-context breakdown, or
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
for more than two groups.

## Usage

``` r
tree_distance(tree_a, tree_b, symmetric = TRUE)
```

## Arguments

- tree_a, tree_b:

  `transitiontrees` objects with the same alphabet.

- symmetric:

  `TRUE` (default) for the symmetric divergence \\\frac{1}{2}(D\_{KL}(A
  \\ B) + D\_{KL}(B \\ A))\\; `FALSE` for \\D\_{KL}(A \\ B)\\ alone.

## Value

A single non-negative number, in nats. It is `Inf` when, at a context
with a positive count, one tree gives probability zero to a next state
that the other tree gives positive probability.

## Conditions

- `transitiontrees_bad_input`:

  (error) either argument is not a `transitiontrees` object.

- `transitiontrees_alphabet_mismatch`:

  (error) the two trees have different alphabets.

## See also

[`compare_trees()`](https://pak.dynasite.org/transitiontrees/reference/compare_trees.md)
for the permutation test;
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md)
for several groups.

## Examples

``` r
high <- context_tree(subset(group_regulation_long, Achiever == "High"),
                     actor = "Actor", time = "Time", action = "Action",
                     max_depth = 2)
low  <- context_tree(subset(group_regulation_long, Achiever == "Low"),
                     actor = "Actor", time = "Time", action = "Action",
                     max_depth = 2)
tree_distance(high, low)
#> [1] 0.0705289
tree_distance(high, low, symmetric = FALSE)
#> [1] 0.06927832
```
