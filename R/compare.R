# ---- compare_trees(): two-tree divergence + permutation test ----

#' @noRd
.pt_lookup_prob <- function(tree, ctx) {
  ## Lookup prob at ctx; fall back to longest-matching-suffix.
  p <- if (ctx %in% names(tree$nodes)) {
    tree$nodes[[ctx]]$prob
  } else {
    parts <- if (identical(ctx, .ROOT)) character(0) else
      strsplit(ctx, " -> ", fixed = TRUE)[[1L]]
    matched <- .ct_match_context(tree, parts)
    tree$nodes[[matched]]$prob
  }
  stats::setNames(as.numeric(p), tree$alphabet)
}

#' @noRd
.pt_lookup_count <- function(tree, ctx) {
  if (ctx %in% names(tree$nodes)) return(tree$nodes[[ctx]]$n)
  0L
}

#' @noRd
.pt_distance_breakdown <- function(tree_a, tree_b) {
  ctx_union <- unique(c(names(tree_a$nodes), names(tree_b$nodes)))
  alphabet <- tree_a$alphabet
  p_a <- lapply(ctx_union, function(c) .pt_lookup_prob(tree_a, c)[alphabet])
  p_b <- lapply(ctx_union, function(c) .pt_lookup_prob(tree_b, c)[alphabet])
  n_a <- vapply(ctx_union, function(c) .pt_lookup_count(tree_a, c),
                numeric(1))
  n_b <- vapply(ctx_union, function(c) .pt_lookup_count(tree_b, c),
                numeric(1))
  KL_ab <- mapply(.ct_kl, p_a, p_b, USE.NAMES = FALSE)
  KL_ba <- mapply(.ct_kl, p_b, p_a, USE.NAMES = FALSE)
  data.frame(
    pathway        = ifelse(ctx_union == .ROOT, .ROOT_LABEL, ctx_union),
    count_a        = n_a, count_b = n_b,
    divergence_ab  = KL_ab,
    divergence_ba  = KL_ba,
    divergence_sym = 0.5 * (KL_ab + KL_ba),
    stringsAsFactors = FALSE
  )
}

#' @noRd
.pt_check_alphabet_compatible <- function(tree_a, tree_b) {
  if (!setequal(tree_a$alphabet, tree_b$alphabet))
    .tt_abort("Trees have incompatible alphabets.",
              class = "transitiontrees_alphabet_mismatch")
}

#' @noRd
.pt_reapply_pruning <- function(tree, template) {
  if (!isTRUE(template$pruned)) return(tree)
  pr <- template$pruning
  if (is.null(pr) || is.null(pr$criterion))
    .tt_abort("Pruned comparison tree is missing pruning metadata.")
  prune_tree(tree, criterion = pr$criterion,
                 alpha = pr$alpha %||% 0.05,
                 threshold = pr$threshold %||% 0.005)
}

#' Compute the divergence between two context trees
#'
#' @description
#' \code{tree_distance()} computes the divergence between two context trees
#' fitted to different sets of sequences, a measure of how differently they
#' predict the next state. The distance is the Kullback-Leibler divergence
#' between the two trees' next-state distributions (symmetric by default),
#' computed for every context present in either tree and averaged with weights
#' equal to the context's combined count in the two trees. A context missing
#' from one tree is evaluated there at its longest fitted suffix, which is the
#' distribution that tree would use for prediction. The value is 0 when the two
#' trees predict identically and grows as their predictions separate.
#'
#' It returns the distance only, in nats. It has no reference distribution, and
#' it is positive even for two samples of one process; use
#' \code{\link{compare_trees}()} for a permutation test and a per-context
#' breakdown, or \code{\link{compare_groups}()} for more than two groups.
#'
#' @param tree_a,tree_b \code{transitiontrees} objects with the same
#'   alphabet.
#' @param symmetric \code{TRUE} (default) for the symmetric divergence
#'   \eqn{\frac{1}{2}(D_{KL}(A \| B) + D_{KL}(B \| A))}{(KL(A || B) +
#'   KL(B || A)) / 2}; \code{FALSE} for \eqn{D_{KL}(A \| B)}{KL(A || B)}
#'   alone.
#' @return A single non-negative number, in nats. It is \code{Inf} when,
#'   at a context with a positive count, one tree gives probability zero
#'   to a next state that the other tree gives positive probability.
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_bad_input}}{(error) either argument is
#'     not a \code{transitiontrees} object.}
#'   \item{\code{transitiontrees_alphabet_mismatch}}{(error) the two trees
#'     have different alphabets.}
#' }
#' @seealso \code{\link{compare_trees}()} for the permutation test;
#'   \code{\link{compare_groups}()} for several groups.
#' @examples
#' high <- context_tree(subset(group_regulation_long, Achiever == "High"),
#'                      actor = "Actor", time = "Time", action = "Action",
#'                      max_depth = 2)
#' low  <- context_tree(subset(group_regulation_long, Achiever == "Low"),
#'                      actor = "Actor", time = "Time", action = "Action",
#'                      max_depth = 2)
#' tree_distance(high, low)
#' tree_distance(high, low, symmetric = FALSE)
#' @export
tree_distance <- function(tree_a, tree_b, symmetric = TRUE) {
  .tt_check_class(tree_a, "transitiontrees", "tree_a")
  .tt_check_class(tree_b, "transitiontrees", "tree_b")
  .pt_check_alphabet_compatible(tree_a, tree_b)

  bd <- .pt_distance_breakdown(tree_a, tree_b)
  per <- if (isTRUE(symmetric)) bd$divergence_sym else bd$divergence_ab
  w   <- bd$count_a + bd$count_b
  total <- sum(w)
  if (total == 0) return(0)
  ## Drop any infinite contributions from zero weights to avoid NaN.
  msk <- w > 0
  sum(per[msk] * w[msk]) / total
}

#' Test whether two context trees differ, with a permutation test
#'
#' @description
#' \code{compare_trees()} tests whether two fitted context trees, one per group,
#' predict differently overall, with a permutation test. The test statistic is
#' the count-weighted symmetric Kullback-Leibler divergence between the trees'
#' next-state distributions, the quantity returned by
#' \code{\link{tree_distance}()}; two groups of sequences can pass through the
#' same states and still differ in what tends to follow a given history.
#'
#' The divergence has no known sampling distribution and is positive even when
#' both groups come from one process, so the reference distribution is built by
#' permutation. The sequences of the two trees are pooled, reassigned at random
#' to two groups of the original sizes, and two trees are refitted with each
#' original tree's depth, minimum count, smoothing and pruning; this is repeated
#' \code{iter} times. The p-value is the share of permuted distances at least as
#' large as the observed one, counting the observed value once (\eqn{(1 + b) /
#' (iter + 1)}). A per-context table shows which histories contribute most to
#' the distance. The result is a \code{transitiontrees_comparison} object, read
#' with \code{as.data.frame()} and \code{summary()} and drawn with
#' \code{plot()}.
#'
#' The comparison is pairwise and treats sequences as exchangeable. For more
#' than two groups, for separate tests of what follows a context and how often
#' it is reached, for per-context p-values with false discovery rate control, or
#' for a permutation stratified by person when one person contributes several
#' sequences, use \code{\link{compare_groups}()}.
#'
#' @param tree_a,tree_b Unweighted \code{transitiontrees} objects fitted to
#'   the two groups. Alternatively a \code{transitiontrees_group} with two
#'   groups as \code{tree_a} and \code{tree_b = NULL} (default); its two
#'   trees are compared in group order.
#' @param iter Number of permutations. Default 200.
#' @param seed Seed for the permutations. Default 1. The call sets the
#'   session's random-number state.
#' @param symmetric \code{TRUE} (default) for the symmetric divergence;
#'   \code{FALSE} for \eqn{D_{KL}(A \| B)}{KL(A || B)} alone.
#'
#' @return An object of class \code{transitiontrees_comparison}. Its
#'   tables are returned by \code{as.data.frame(x, what = "pathways")}
#'   (one row per context in either tree: counts in each tree and the two
#'   directed and the symmetric divergence) and
#'   \code{as.data.frame(x, what = "null")} (one row per permutation);
#'   \code{summary()} gives one row with the observed distance, the null
#'   mean, standard deviation and 95th percentile, and the p-value;
#'   \code{plot()} draws the null distribution.
#'
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_bad_input}}{(error) an argument is not a
#'     \code{transitiontrees} object; \code{tree_b} missing; a
#'     \code{transitiontrees_group} without exactly two groups, or given
#'     together with \code{tree_b}; \code{iter} not a positive number.}
#'   \item{\code{transitiontrees_weighted_unsupported}}{(error) a tree was
#'     fitted with weights; the permutation refits from the raw sequences
#'     and would drop them.}
#'   \item{\code{transitiontrees_alphabet_mismatch}}{(error) the trees have
#'     different alphabets.}
#' }
#'
#' @seealso \code{\link{summary.transitiontrees_comparison}},
#'   \code{\link{as.data.frame.transitiontrees_comparison}} and
#'   \code{\link{plot.transitiontrees_comparison}} for the result;
#'   \code{\link{tree_distance}()} for the distance alone;
#'   \code{\link{compare_groups}()} for several groups and per-context
#'   tests.
#'
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' cmp <- compare_trees(grp, iter = 20)
#' cmp
#' summary(cmp)
#' @export
compare_trees <- function(tree_a, tree_b = NULL, iter = 200L,
                              seed = 1L, symmetric = TRUE) {
  ## A 2-group transitiontrees_group as the sole argument -> compare its pair.
  if (inherits(tree_a, "transitiontrees_group")) {
    if (!is.null(tree_b))
      .tt_abort("Pass a 2-group 'transitiontrees_group' as the only tree argument, ",
           "not a group plus 'tree_b'.")
    if (length(tree_a) != 2L)
      .tt_abort("compare_trees() is pairwise; the transitiontrees_group must ",
           "have exactly 2 groups (got ", length(tree_a), ").")
    grp    <- tree_a
    tree_a <- grp[[1L]]
    tree_b <- grp[[2L]]
  }
  if (is.null(tree_b))
    .tt_abort("'tree_b' is required (or pass a 2-group 'transitiontrees_group' as ",
         "the first argument).")
  .tt_check_class(tree_a, "transitiontrees", "tree_a")
  .tt_check_class(tree_b, "transitiontrees", "tree_b")
  .pt_assert_unweighted(tree_a, "compare_trees()")
  .pt_assert_unweighted(tree_b, "compare_trees()")
  if (!is.numeric(iter) || length(iter) != 1L ||
      is.na(iter) || iter < 1)
    .tt_abort("'iter' must be a single positive integer.")
  .pt_check_alphabet_compatible(tree_a, tree_b)

  pdist_obs <- tree_distance(tree_a, tree_b, symmetric = symmetric)
  bd        <- .pt_distance_breakdown(tree_a, tree_b)

  pooled  <- c(tree_a$data, tree_b$data)
  n_a     <- length(tree_a$data)
  n_total <- length(pooled)
  alphabet <- tree_a$alphabet

  ## Refit using each observed tree's hyperparameters and pruning state.
  set.seed(seed)
  null_dist <- vapply(seq_len(iter), function(b) {
    perm  <- sample.int(n_total)
    a_idx <- perm[seq_len(n_a)]
    b_idx <- perm[(n_a + 1L):n_total]
    tr_a  <- context_tree(pooled[a_idx],
                          max_depth = tree_a$max_depth,
                          min_count = tree_a$nmin,
                          smoothing = tree_a$smoothing,
                          alphabet  = alphabet)
    tr_b  <- context_tree(pooled[b_idx],
                          max_depth = tree_b$max_depth,
                          min_count = tree_b$nmin,
                          smoothing = tree_b$smoothing,
                          alphabet  = alphabet)
    tr_a <- .pt_reapply_pruning(tr_a, tree_a)
    tr_b <- .pt_reapply_pruning(tr_b, tree_b)
    tree_distance(tr_a, tr_b, symmetric = symmetric)
  }, numeric(1))

  p_value <- (1 + sum(null_dist >= pdist_obs)) / (iter + 1)

  pathways <- bd[order(-bd$divergence_sym), , drop = FALSE]
  rownames(pathways) <- NULL
  structure(
    list(pdist     = pdist_obs,
         null_dist = null_dist,
         p_value   = p_value,
         pathways  = pathways),
    class = "transitiontrees_comparison"
  )
}

#' Plot the permutation distribution of a two-tree comparison
#'
#' @description
#' \code{plot()} draws the permutation null distribution of the distance between
#' two context trees as a histogram, with the observed distance as a vertical
#' line and the permutation p-value in the subtitle. An observed line far in the
#' right tail indicates that the two groups predict differently beyond what
#' random reassignment of sequences produces.
#'
#' @param x A \code{transitiontrees_comparison} object from
#'   \code{\link{compare_trees}()}.
#' @param bins Number of histogram bins. Default 30.
#' @param ... Ignored.
#' @return A ggplot object: a histogram of the permuted distances with the
#'   observed distance marked.
#' @seealso \code{\link{compare_trees}()},
#'   \code{\link{as.data.frame.transitiontrees_comparison}} for the
#'   permuted values.
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' plot(compare_trees(grp, iter = 20), bins = 10)
#' @export
plot.transitiontrees_comparison <- function(x, bins = 30L, ...) {
  null_df <- data.frame(d = x$null_dist)
  ggplot2::ggplot(null_df, ggplot2::aes(x = .data$d)) +
    ggplot2::geom_histogram(bins = bins, fill = "#0072B2",
                             colour = "white", alpha = 0.85) +
    ggplot2::geom_vline(xintercept = x$pdist, colour = "#D55E00",
                        linewidth = 1.1) +
    ggplot2::annotate("text", x = x$pdist, y = Inf,
                       label = sprintf("observed = %.3f",
                                        x$pdist),
                       hjust = -0.05, vjust = 1.4, colour = "#D55E00",
                       size = 3.6, fontface = "bold") +
    ggplot2::labs(
      x = "symmetric KL between trees", y = "permutations",
      title = "context tree comparison: observed vs. null distribution",
      subtitle = sprintf("permutation p-value = %.3f  (n = %d)",
                          x$p_value, length(x$null_dist))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40",
                                             size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Convert a two-tree comparison to a data frame
#'
#' @description
#' \code{as.data.frame()} returns the per-context divergence table or the
#' permutation null distribution of a \code{\link{compare_trees}()} result.
#'
#' @param x A \code{transitiontrees_comparison} object.
#' @param row.names,optional Ignored.
#' @param what \code{"pathways"} (default) for the per-context table, or
#'   \code{"null"} for the permuted distances.
#' @param ... Ignored.
#' @return A data frame. For \code{what = "pathways"}, one row per context
#'   present in either tree, sorted by \code{divergence_sym} (largest
#'   first), with columns \code{pathway}, \code{count_a} and
#'   \code{count_b} (the context's count in each tree, 0 when absent),
#'   \code{divergence_ab} and \code{divergence_ba} (the two directed
#'   Kullback-Leibler divergences, nats) and \code{divergence_sym} (their
#'   mean). For \code{what = "null"}, one row per permutation with columns
#'   \code{permutation} and \code{distance}.
#' @seealso \code{\link{compare_trees}()},
#'   \code{\link{summary.transitiontrees_comparison}}.
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' cmp <- compare_trees(grp, iter = 20)
#' as.data.frame(cmp)
#' head(as.data.frame(cmp, what = "null"))
#' @export
as.data.frame.transitiontrees_comparison <- function(x, row.names = NULL,
                                               optional = FALSE,
                                               what = c("pathways", "null"),
                                               ...) {
  what <- match.arg(what)
  switch(what,
    pathways = x$pathways,
    null     = data.frame(permutation = seq_along(x$null_dist),
                          distance    = as.numeric(x$null_dist)))
}

#' Summarize a two-tree comparison
#'
#' @description
#' \code{summary()} reduces a \code{\link{compare_trees}()} result to one row:
#' the observed distance between the two trees, the mean, standard deviation and
#' 95th percentile of the permutation null distribution, and the permutation
#' p-value.
#'
#' @param object A \code{transitiontrees_comparison} object.
#' @param ... Ignored.
#' @return An object of class \code{summary.transitiontrees_comparison}:
#'   a one-row data frame with columns \code{distance}, \code{null_mean},
#'   \code{null_sd}, \code{null_q95}, \code{p_value}, \code{iter} and
#'   \code{n_pathways} (contexts in the per-context table).
#' @seealso \code{\link{compare_trees}()},
#'   \code{\link{as.data.frame.transitiontrees_comparison}}.
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' summary(compare_trees(grp, iter = 20))
#' @export
summary.transitiontrees_comparison <- function(object, ...) {
  .tt_check_class(object, "transitiontrees_comparison", "object")
  nd <- as.numeric(object$null_dist)
  out <- data.frame(
    distance   = object$pdist,
    null_mean  = mean(nd),
    null_sd    = if (length(nd) > 1L) stats::sd(nd) else NA_real_,
    null_q95   = unname(stats::quantile(nd, 0.95, names = FALSE)),
    p_value    = object$p_value,
    iter       = length(nd),
    n_pathways = nrow(object$pathways),
    row.names  = NULL)
  class(out) <- c("summary.transitiontrees_comparison", "data.frame")
  out
}

#' Print a summary of a two-tree comparison
#'
#' @description
#' \code{print()} prints the one-row summary of a \code{\link{compare_trees}()}
#' result as a table.
#'
#' @param x A \code{summary.transitiontrees_comparison} object.
#' @param digits Significant digits for numeric columns. Default 3.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{summary.transitiontrees_comparison}}.
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' print(summary(compare_trees(grp, iter = 20)), digits = 2)
#' @export
print.summary.transitiontrees_comparison <- function(x, digits = 3L, ...) {
  cat("<transitiontrees_comparison summary>\n")
  print.data.frame(as.data.frame(unclass(x)), row.names = FALSE,
                   digits = digits)
  invisible(x)
}

#' Print a two-tree comparison
#'
#' @description
#' \code{print()} prints the number of permutations, the observed distance, the
#' mean of the null distribution, the p-value, and the \code{n} contexts that
#' contribute most to the distance.
#'
#' @param x A \code{transitiontrees_comparison} object.
#' @param digits Significant digits for the distances and rounding digits
#'   for the per-context table. Default 3.
#' @param n Number of contexts to print. Default 6.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{compare_trees}()},
#'   \code{\link{as.data.frame.transitiontrees_comparison}} for every
#'   context.
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' print(compare_trees(grp, iter = 20), n = 3)
#' @export
print.transitiontrees_comparison <- function(x, digits = 3L, n = 6L, ...) {
  cat(sprintf("<transitiontrees_comparison>  iter = %d\n",
              length(x$null_dist)))
  cat(sprintf("  observed distance : %s\n",
              format(x$pdist, digits = digits)))
  cat(sprintf("  null mean         : %s\n",
              format(mean(x$null_dist), digits = digits)))
  cat(sprintf("  p-value           : %s\n",
              format(x$p_value, digits = digits)))
  cat("\ntop divergent pathways:\n")
  bd <- utils::head(x$pathways, n)
  bd$divergence_ab  <- round(bd$divergence_ab, digits)
  bd$divergence_ba  <- round(bd$divergence_ba, digits)
  bd$divergence_sym <- round(bd$divergence_sym, digits)
  print.data.frame(bd, row.names = FALSE)
  invisible(x)
}
