# Plot Each Tree in a context tree Group

Draw every member of a `transitiontrees_group` in turn via
[`plot.transitiontrees`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees.md).
Each member's plot is printed (so the call produces one figure per
group, e.g. in an R Markdown chunk), captioned with its group name; the
named list of plot objects is returned invisibly for further use.

## Usage

``` r
# S3 method for class 'transitiontrees_group'
plot(x, ...)
```

## Arguments

- x:

  A `transitiontrees_group`.

- ...:

  Passed to
  [`plot.transitiontrees`](https://pak.dynasite.org/transitiontrees/reference/plot.transitiontrees.md)
  (e.g. `style`).

## Value

Invisibly, a named list (one entry per group, in group order) of the
per-member plot objects.

## Examples

``` r
# \donttest{
m   <- matrix(sample(c("A","B","C"), 200, replace = TRUE), 40, 5)
grp <- context_tree(m, group = rep(c("x","y"), each = 20),
                    max_depth = 2L)
plot(grp)                         # one tree per group


# }
```
