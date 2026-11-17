# Print the summary of a group of transition trees

[`print()`](https://rdrr.io/r/base/print.html) prints the table returned
by
[`summary.transitiontrees_group()`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_group.md)
under a one-line header, with numeric columns rounded.

## Usage

``` r
# S3 method for class 'summary.transitiontrees_group'
print(x, digits = 3L, ...)
```

## Arguments

- x:

  A `summary.transitiontrees_group` object.

- digits:

  Number of significant digits for numeric columns. Default 3.

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`summary.transitiontrees_group()`](https://pak.dynasite.org/transitiontrees/reference/summary.transitiontrees_group.md).

## Examples

``` r
groups <- context_tree(group_regulation_long, actor = "Actor",
                       time = "Time", action = "Action",
                       group = "Achiever", max_depth = 2)
print(summary(groups), digits = 2)
#> <transitiontrees_group summary>  2 groups
#> 
#>  group n_seq n_obs n_nodes max_depth pruned logLik  df  nobs   AIC   BIC
#>   High  1000 13721      80         2  FALSE -22773 640 13721 46826 51643
#>    Low  1000 13812      82         2  FALSE -22968 656 13812 47248 52190
#>  perplexity
#>         5.3
#>         5.3
```
