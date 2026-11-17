# Print a group of fitted transition trees

[`print()`](https://rdrr.io/r/base/print.html) prints one line per group
of a `transitiontrees_group`: the number of contexts, the depth of the
tree, and the numbers of sequences and observations it was fitted on.

## Usage

``` r
# S3 method for class 'transitiontrees_group'
print(x, ...)
```

## Arguments

- x:

  A `transitiontrees_group` object.

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`summary.transitiontrees_group()`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_group.md)
for the fit of each group;
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md).

## Examples

``` r
groups <- context_tree(group_regulation_long, actor = "Actor",
                       time = "Time", action = "Action",
                       group = "Achiever", max_depth = 2)
print(groups)
#> <transitiontrees_group>  2 groups
#>   High            80 nodes, depth <= 2, 1000 seq, 13721 obs
#>   Low             82 nodes, depth <= 2, 1000 seq, 13812 obs
```
