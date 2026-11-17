# ---- Pruning ----

#' Prune a transition tree to the contexts that change the prediction
#'
#' @description
#' \code{prune_tree()} removes the contexts of a fitted tree that predict the
#' next state no differently from their shorter suffix, so that the long
#' contexts that remain are those in which the additional history changes the
#' prediction. This is the principle behind the learning of prediction suffix
#' trees (Ron, Singer & Tishby, 1996) and the context algorithm for
#' variable-length Markov chains (\enc{Bühlmann}{Buhlmann} & Wyner, 1999). A
#' fitted tree holds every context that occurs often enough, and many such
#' contexts add parameters without adding information about what comes next.
#'
#' It compares each leaf with its parent, the context one state shorter, and
#' removes the leaf when the comparison does not pass the chosen criterion;
#' leaves exposed by a removal are tested in turn until no leaf can be removed.
#' The root is never removed. With \code{"G2"} a leaf is kept when the
#' likelihood-ratio statistic of its counts against the parent's distribution
#' exceeds the chi-square critical value at \code{alpha} with \eqn{k - 1}
#' degrees of freedom (\eqn{k} states); each leaf is tested separately, without
#' correction for multiple tests. With \code{"KL"} a leaf is kept when its
#' Kullback-Leibler divergence from the parent, in nats, exceeds
#' \code{threshold}. With \code{"AIC"} and \code{"BIC"} a leaf is kept when its
#' log-likelihood, penalized by \eqn{2(k - 1)} or \eqn{\log(n)(k - 1)} (\eqn{n}
#' the leaf's count), is better than that of the parent's prediction (Akaike,
#' 1974; Schwarz, 1978). The probabilities of the contexts that remain are not
#' changed. It returns a tree of the same class with fewer contexts.
#'
#' To see how strongly each criterion prunes a tree, use
#' \code{\link{compare_pruning}()}; to choose between pruned and unpruned
#' trees on held-out data, use \code{\link{tune_tree}()}; to draw the
#' \code{"G2"} decision along the suffix chain of one pathway, use
#' \code{\link{plot_pruning}()}.
#' The \code{divergence} column of \code{\link{tree_pathways}()} is in bits,
#' whereas the \code{"KL"} threshold is in nats; divide a threshold in nats
#' by \code{log(2)} to read it in bits.
#'
#' @param tree A \code{transitiontrees} object, or a
#'   \code{transitiontrees_group}, whose trees are each pruned.
#' @param criterion \code{"G2"} (default), \code{"KL"}, \code{"AIC"} or
#'   \code{"BIC"}.
#' @param alpha Significance level for \code{"G2"}, in (0, 1). Default
#'   0.05.
#' @param threshold Minimum divergence in nats for \code{"KL"}. Default
#'   0.005.
#'
#' @return An object of the same class as \code{tree}, with fewer contexts.
#'   A pruned tree prints \code{[pruned]} and the criterion in its header,
#'   and \code{\link{n_nodes}()} gives its size.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree,
#' \code{alpha} is not a single number in (0, 1), or \code{threshold} is
#' not a single non-negative number.
#'
#' @references
#' Akaike, H. (1974). A new look at the statistical model identification.
#' \emph{IEEE Transactions on Automatic Control}, 19(6), 716-723.
#' \doi{10.1109/TAC.1974.1100705}
#'
#' \enc{Bühlmann}{Buhlmann}, P., & Wyner, A. J. (1999). Variable length
#' Markov chains. \emph{The Annals of Statistics}, 27(2), 480-513.
#' \doi{10.1214/aos/1018031204}
#'
#' Ron, D., Singer, Y., & Tishby, N. (1996). The power of amnesia: Learning
#' probabilistic automata with variable memory length. \emph{Machine
#' Learning}, 25(2-3), 117-149. \doi{10.1023/A:1026490906255}
#'
#' Schwarz, G. (1978). Estimating the dimension of a model. \emph{The
#' Annals of Statistics}, 6(2), 461-464. \doi{10.1214/aos/1176344136}
#'
#' @seealso \code{\link{compare_pruning}()}, \code{\link{plot_pruning}()},
#'   \code{\link{tune_tree}()}, \code{\link{context_tree}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 4)
#' prune_tree(tree, criterion = "G2", alpha = 0.05)
#' n_nodes(prune_tree(tree, criterion = "BIC"))
#' @export
prune_tree <- function(tree,
                           criterion = c("G2", "KL", "AIC", "BIC"),
                           alpha = 0.05, threshold = 0.005) {
  ## A transitiontrees_group prunes each member, preserving the group wrapper.
  if (inherits(tree, "transitiontrees_group")) {
    out <- lapply(tree, prune_tree, criterion = criterion,
                  alpha = alpha, threshold = threshold)
    return(structure(out, class = class(tree),
                     group = attr(tree, "group")))
  }
  .tt_check_class(tree, "transitiontrees", "tree")
  criterion <- match.arg(criterion)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1)
    .tt_abort("'alpha' must be a single number in (0, 1).")
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
      threshold < 0)
    .tt_abort("'threshold' must be a single non-negative number.")

  ## Critical value for G2: chi-square at (k-1) df
  k_states <- length(tree$alphabet)
  g2_cv    <- stats::qchisq(1 - alpha, df = k_states - 1L)

  nodes <- tree$nodes
  edges <- tree$edges
  alpha_size <- k_states

  ## Repeat until no further removals
  repeat {
    leaves <- setdiff(names(nodes), edges$parent)
    leaves <- leaves[leaves != .ROOT]  # never the root
    drop   <- character(0)

    for (leaf in leaves) {
      parts  <- strsplit(leaf, " -> ", fixed = TRUE)[[1L]]
      parent <- if (length(parts) == 1L) .ROOT else
        paste(parts[-1L], collapse = " -> ")
      if (!parent %in% names(nodes)) next

      child  <- nodes[[leaf]]
      par    <- nodes[[parent]]
      keep <- switch(
        criterion,
        G2  = .ct_g2(child$counts, par$prob) > g2_cv,
        KL  = .ct_kl(child$prob, par$prob) > threshold,
        AIC = {
          ll_child  <- sum(child$counts * log(pmax(child$prob,
                                                   .Machine$double.eps)))
          ll_parent <- sum(child$counts * log(pmax(par$prob,
                                                   .Machine$double.eps)))
          aic_child  <- 2 * (alpha_size - 1L) - 2 * ll_child
          aic_parent <- -2 * ll_parent
          aic_child < aic_parent
        },
        BIC = {
          ll_child  <- sum(child$counts * log(pmax(child$prob,
                                                   .Machine$double.eps)))
          ll_parent <- sum(child$counts * log(pmax(par$prob,
                                                   .Machine$double.eps)))
          bic_child  <- log(child$n) * (alpha_size - 1L) - 2 * ll_child
          bic_parent <- -2 * ll_parent
          bic_child < bic_parent
        }
      )
      if (!keep) drop <- c(drop, leaf)
    }
    if (length(drop) == 0L) break
    nodes <- nodes[!names(nodes) %in% drop]
    edges <- edges[!edges$child %in% drop, , drop = FALSE]
  }

  tree$nodes  <- nodes
  tree$edges  <- edges
  tree$pruned <- TRUE
  tree$pruning <- list(criterion = criterion, alpha = alpha,
                       threshold = threshold)
  tree
}

#' Compare pruning criteria on one tree
#'
#' @description
#' \code{compare_pruning()} prunes one fitted tree under each of several
#' criteria and reports how many contexts each keeps. The table shows how
#' strongly the criteria differ on the data at hand before one of them is
#' applied with \code{\link{prune_tree}()}; the criteria themselves are
#' described there. A smaller tree is not by itself a better one: to compare
#' pruned trees by how well they predict new sequences, use
#' \code{\link{tune_tree}()}.
#'
#' @param tree A \code{transitiontrees} object, usually unpruned.
#' @param criterion Criteria to compare. Default all four: \code{"G2"},
#'   \code{"KL"}, \code{"AIC"}, \code{"BIC"}.
#' @param alpha Significance level for \code{"G2"}. Default 0.05.
#' @param threshold Minimum divergence in nats for \code{"KL"}. Default
#'   0.005.
#'
#' @return A \code{data.frame} with one row per criterion, in the order of
#'   \code{criterion}, and columns \code{criterion}, \code{n_nodes} (contexts
#'   kept, root included) and \code{reduction_pct} (percentage of the
#'   original contexts removed).
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree,
#' \code{criterion} is not a non-empty character vector, or \code{alpha}
#' or \code{threshold} is invalid.
#'
#' @seealso \code{\link{prune_tree}()}, \code{\link{plot_pruning}()},
#'   \code{\link{tune_tree}()}, \code{\link{compare_smoothing}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 4)
#' compare_pruning(tree)
#' compare_pruning(tree, criterion = c("G2", "BIC"), alpha = 0.01)
#' @export
compare_pruning <- function(tree,
                            criterion = c("G2", "KL", "AIC", "BIC"),
                            alpha = 0.05, threshold = 0.005) {
  .tt_check_class(tree, "transitiontrees", "tree")
  if (!is.character(criterion) || length(criterion) < 1L)
    .tt_abort("'criterion' must be a non-empty character vector of ",
         "criterion names.")
  n0 <- length(tree$nodes)
  n  <- vapply(criterion, function(cr)
    length(prune_tree(tree, criterion = cr, alpha = alpha,
                          threshold = threshold)$nodes), integer(1))
  data.frame(
    criterion     = criterion,
    n_nodes       = n,
    reduction_pct = round(100 * (1 - n / n0), 1),
    stringsAsFactors = FALSE,
    row.names     = NULL)
}
