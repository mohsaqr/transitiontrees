# Collaborative-regulation events (long format)

A long, one-row-per-event log of collaborative regulation moves, with
timestamps. Used to demonstrate long-format loading: reshape with
[`prepare_input()`](https://pak.dynasite.org/transitiontrees/reference/prepare_input.md)
(or name `action =` in
[`context_tree()`](https://pak.dynasite.org/transitiontrees/reference/context_tree.md))
to split each actor's events into time-gap sessions. Bundled example
dataset.

## Usage

``` r
group_regulation_long
```

## Format

A `data.frame` with 27533 rows and 6 columns:

- Actor:

  integer; the learner.

- Achiever:

  character; an achievement-level covariate.

- Group:

  numeric; the collaboration group.

- Course:

  character; the course.

- Time:

  POSIXct; the event timestamp.

- Action:

  character; the regulation move (the state).

## Source

Bundled example dataset.

## Examples

``` r
data(group_regulation_long)
context_tree(group_regulation_long, actor = "Actor", time = "Time",
             action = "Action", max_depth = 3L)
#> <transitiontrees>  471 nodes, depth <= 3, 9 states  [unpruned]
#>   alphabet : adapt, cohesion, consensus, coregulate, discuss, emotion, monitor, plan, synthesis
#>   fit on   : 2000 sequences, 27533 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=27533  -> consensus (0.25)
#> |-- adapt     n=509    -> consensus (0.47)
#> |   |-- consensus  n=27     -> cohesion (0.48)
#> |   |   |-- discuss   n=6      -> consensus (0.66)
#> |   |   `-- plan      n=7      -> consensus (0.71)
#> |   |-- coregulate  n=28     -> consensus (0.50)
#> |   |   `-- consensus  n=21     -> consensus (0.43)
#> |   |-- discuss   n=259    -> consensus (0.47)
#> |   |   |-- adapt     n=6      -> cohesion (0.50)
#> |   |   |-- cohesion  n=8      -> consensus (0.37)
#> |   |   |-- consensus  n=60     -> consensus (0.53)
#> |   |   |-- coregulate  n=37     -> consensus (0.43)
#> |   |   |-- discuss   n=48     -> consensus (0.52)
#> |   |   |-- emotion   n=14     -> consensus (0.50)
#> |   |   |-- monitor   n=25     -> consensus (0.40)
#> |   |   `-- plan      n=26     -> consensus (0.42)
#> |   |-- emotion   n=6      -> consensus (0.83)
#> |   |-- monitor   n=16     -> consensus (0.50)
#> |   |   `-- plan      n=5      -> cohesion (0.40)
#> |   |-- plan      n=6      -> consensus (0.66)
#> |   `-- synthesis  n=140    -> consensus (0.48)
#> |       |-- consensus  n=8      -> cohesion (0.37)
#> |       |-- coregulate  n=5      -> consensus (0.60)
#> |       |-- discuss   n=107    -> consensus (0.45)
#> |       `-- monitor   n=6      -> consensus (0.66)
#> |-- cohesion  n=1695   -> consensus (0.50) 
#> ... 445 more nodes (use as.data.frame(x) or summary(x))
```
