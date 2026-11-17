# Summarize a group of fitted transition trees

[`summary()`](https://rdrr.io/r/base/summary.html) summarizes a
`transitiontrees_group` with one row per group: the size of the data
each tree was fitted on, the size and depth of the tree, and its
in-sample fit (log-likelihood, number of parameters, AIC, BIC and
perplexity, as returned by
[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md)).
The in-sample fit describes each tree; to test whether the groups
differ, use
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md).

## Usage

``` r
# S3 method for class 'transitiontrees_group'
summary(object, ...)
```

## Arguments

- object:

  A `transitiontrees_group` object.

- ...:

  Not used.

## Value

An object of class `summary.transitiontrees_group`, a `data.frame` with
one row per group and columns `group`, `n_seq`, `n_obs`, `n_nodes`,
`max_depth`, `pruned`, `logLik`, `df`, `nobs`, `AIC`, `BIC` and
`perplexity`.

## See also

[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md),
[`model_fit()`](https://pak.dynasite.org/transitiontrees/reference/model_fit.md),
[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md).

## Examples

``` r
groups <- context_tree(group_regulation_long, actor = "Actor",
                       time = "Time", action = "Action",
                       group = "Achiever", max_depth = 2)
summary(groups)
#> <transitiontrees_group summary>  2 groups
#> 
#>  group n_seq n_obs n_nodes max_depth pruned logLik  df  nobs   AIC   BIC
#>   High  1000 13721      80         2  FALSE -22773 640 13721 46826 51643
#>    Low  1000 13812      82         2  FALSE -22968 656 13812 47248 52190
#>  perplexity
#>        5.26
#>        5.27
```
