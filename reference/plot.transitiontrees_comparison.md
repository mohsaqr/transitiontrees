# Plot a context tree Comparison

Visualises the permutation-test result. Histogram of the null
distribution of transitiontrees distances under shuffled group labels; a
vertical line marks the observed distance; the panel header carries the
p-value.

## Usage

``` r
# S3 method for class 'transitiontrees_comparison'
plot(x, bins = 30L, ...)
```

## Arguments

- x:

  A `transitiontrees_comparison` object.

- bins:

  Integer. Histogram bins. Default 30.

- ...:

  Ignored.

## Value

A ggplot object.
