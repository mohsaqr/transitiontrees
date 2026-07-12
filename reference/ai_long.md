# AI-collaboration messages (long format)

A long, one-row-per-message log from an AI-assisted collaboration study,
with Unix timestamps and an explicit session id. Used to demonstrate
long-format loading with Unix time and sessions. Bundled example
dataset.

## Usage

``` r
ai_long
```

## Format

A `data.frame` with 8551 rows and 9 columns, including `project`,
`session_id`, `timestamp` (Unix seconds), `code` / `cluster` (the state
at two granularities), and `code_order` / `order_in_session`
(within-sequence order).

## Source

Bundled example dataset.

## Examples

``` r
data(ai_long)
context_tree(ai_long, actor = "project", time = "timestamp",
             action = "code", max_depth = 2L)
#> <transitiontrees>  62 nodes, depth <= 2, 8 states  [unpruned]
#>   alphabet : Ask, Delegate, Execute, Explain, Investigate, Plan, Repair, Report
#>   fit on   : 160 sequences, 8551 observations
#>   smoothing: floor(ymin=0.001, rule=interpolate)   min_count = 5
#> (start)   n=8551   -> Execute (0.38)
#> |-- Ask       n=99     -> Explain (0.44)
#> |   |-- Execute   n=32     -> Execute (0.37)
#> |   |-- Explain   n=9      -> Explain (0.66)
#> |   |-- Investigate  n=21     -> Explain (0.43)
#> |   |-- Plan      n=15     -> Explain (0.40)
#> |   `-- Repair    n=13     -> Explain (0.84)
#> |-- Delegate  n=290    -> Plan (0.60)
#> |   |-- Ask       n=5      -> Plan (0.60)
#> |   |-- Delegate  n=9      -> Plan (0.55)
#> |   |-- Execute   n=88     -> Plan (0.59)
#> |   |-- Explain   n=7      -> Plan (0.43)
#> |   |-- Investigate  n=43     -> Plan (0.58)
#> |   |-- Plan      n=88     -> Plan (0.65)
#> |   |-- Repair    n=11     -> Execute (0.54)
#> |   `-- Report    n=12     -> Plan (0.41)
#> |-- Execute   n=3188   -> Execute (0.49)
#> |   |-- Ask       n=23     -> Execute (0.52)
#> |   |-- Delegate  n=51     -> Execute (0.37)
#> |   |-- Execute   n=1537   -> Execute (0.53)
#> |   |-- Explain   n=129    -> Investigate (0.37)
#> |   |-- Investigate  n=545    -> Execute (0.49)
#> |   |-- Plan      n=695    -> Execute (0.44)
#> |   |-- Repair    n=118    -> Execute (0.45)
#> |   `-- Report    n=62     -> Execute (0.53)
#> |-- Explain   n=512    -> Investigate (0.32) 
#> ... 36 more nodes (use as.data.frame(x) or summary(x))
```
