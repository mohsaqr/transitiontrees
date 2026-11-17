# Regulation of collaborative learning (long event log)

`group_regulation_long` is a coded dataset of the regulation actions of
2,000 students working in 200 small groups across three courses, one row
per action with a timestamp. Each action is one of nine regulation
moves: `adapt`, `cohesion`, `consensus`, `coregulate`, `discuss`,
`emotion`, `monitor`, `plan` and `synthesis`. Each student is classified
as a high or low achiever (1,000 each), which makes the data suited to
comparing the regulation pathways of two groups.

The data are in long format. Name the columns in
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md)
(`actor = "Actor"`, `time = "Time"`, `action = "Action"`) or reshape
them with
[`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md);
with the default 15-minute session gap, each student's actions form one
sequence of 2 to 26 actions.

## Usage

``` r
group_regulation_long
```

## Format

A `data.frame` with 27,533 rows (actions) and 6 columns:

- Actor:

  Integer student identifier, 1 to 2,000.

- Achiever:

  Achievement level, `"High"` or `"Low"`, constant within a student.

- Group:

  Numeric collaboration group, 1 to 200; each student belongs to one
  group.

- Course:

  Course, `"A"`, `"B"` or `"C"`.

- Time:

  `POSIXct` timestamp of the action.

- Action:

  The regulation move, the state of the sequence.

## Source

Identical to `group_regulation_long` in the tna package (there a
tibble), where it is described as generated synthetically from that
package's `group_regulation` data.

## See also

[`compare_groups()`](https://pak.dynasite.org/transitiontrees/reference/compare_groups.md),
[`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md).

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
