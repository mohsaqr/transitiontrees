# ---- tree_dependence(): per-context KL diagnostic ----

#' @noRd
.pt_empty_dependence_df <- function() {
  data.frame(
    pathway            = character(0),
    depth              = integer(0),
    count              = numeric(0),
    divergence         = numeric(0),
    entropy            = numeric(0),
    entropy_before     = numeric(0),
    entropy_drop       = numeric(0),
    likely_next        = character(0),
    likely_before      = character(0),
    changes_prediction = logical(0),
    stringsAsFactors   = FALSE
  )
}

#' Measure how much each step of history changes the prediction
#'
#' @description
#' \code{tree_dependence()} measures how much each step of history changes the
#' prediction of a fitted tree, that is, where the process is path dependent. A
#' process is path dependent where what happens next depends on more than the
#' current state. In a transition tree this shows as a context whose next-state
#' distribution differs from that of its parent, the same context without its
#' oldest state. The difference has two sides: how far the distribution moves
#' (its divergence) and whether it becomes more or less certain (its change in
#' entropy).
#'
#' It returns one row for every context other than the root, with the
#' Kullback-Leibler divergence of its next-state distribution from its parent's,
#' the Shannon entropy of both distributions and their difference, and the most
#' likely next state under each (Cover & Thomas, 2006). A positive
#' \code{entropy_drop} means the additional state of history makes the next
#' state more predictable; a negative one means it makes it less predictable.
#' The divergence is the quantity that \code{\link{prune_tree}()} compares with
#' its threshold under \code{criterion = "KL"} (there in nats).
#'
#' For the pathway table with counts and ranked views, use
#' \code{\link{tree_pathways}()} and \code{\link{divergent_pathways}()}; to
#' see the full next-state distributions along one pathway, use
#' \code{\link{plot_pruning}()}.
#'
#' @param tree A \code{transitiontrees} object.
#' @param base Base of the logarithm: 2 (default) for bits, \code{exp(1)}
#'   for nats.
#' @param sort_by Column to sort by, in decreasing order:
#'   \code{"divergence"} (default), \code{"entropy_drop"}, \code{"entropy"},
#'   \code{"count"} or \code{"depth"}.
#' @param top Number of rows to keep after sorting. \code{NULL} (default)
#'   keeps all.
#'
#' @return A \code{data.frame} with one row per non-root context whose
#'   parent is in the tree, and columns \code{pathway}, \code{depth},
#'   \code{count}, \code{divergence} (from the parent's distribution),
#'   \code{entropy} (of this context's distribution), \code{entropy_before}
#'   (of the parent's), \code{entropy_drop} (\code{entropy_before} minus
#'   \code{entropy}), \code{likely_next}, \code{likely_before} (the parent's
#'   most likely next state) and \code{changes_prediction}. A tree with
#'   only a root returns a zero-row table with these columns.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree.
#'
#' @references
#' Cover, T. M., & Thomas, J. A. (2006). \emph{Elements of information
#' theory} (2nd ed.). Wiley.
#'
#' @seealso \code{\link{tree_pathways}()}, \code{\link{divergent_pathways}()},
#'   \code{\link{prune_tree}()}, \code{\link{plot_divergence}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' tree_dependence(tree, top = 5)
#' tree_dependence(prune_tree(tree), sort_by = "entropy_drop")
#' @export
tree_dependence <- function(tree, base = 2,
                                sort_by = c("divergence", "entropy_drop",
                                            "entropy", "count", "depth"),
                                top = NULL) {
  .tt_check_class(tree, "transitiontrees", "tree")
  stopifnot(
    "`base` must be a single positive number other than 1" =
      is.numeric(base) && length(base) == 1L && !is.na(base) &&
      base > 0 && abs(base - 1) > sqrt(.Machine$double.eps)
  )
  sort_by <- match.arg(sort_by)

  H <- function(p) {
    msk <- p > 0
    -sum(p[msk] * log(p[msk], base = base))
  }
  KL <- function(p, q) {
    msk <- p > 0
    if (any(p[msk] > 0 & q[msk] == 0)) return(Inf)
    sum(p[msk] * log(p[msk] / q[msk], base = base))
  }

  rows <- lapply(setdiff(names(tree$nodes), .ROOT), function(ctx) {
    parts  <- strsplit(ctx, " -> ", fixed = TRUE)[[1L]]
    parent <- if (length(parts) == 1L) .ROOT else
      paste(parts[-1L], collapse = " -> ")
    if (!parent %in% names(tree$nodes)) return(NULL)

    node <- tree$nodes[[ctx]]
    par  <- tree$nodes[[parent]]
    p    <- node$prob
    q    <- par$prob

    modal_next   <- tree$alphabet[which.max(p)]
    modal_parent <- tree$alphabet[which.max(q)]

    data.frame(
      pathway            = ctx,
      depth              = node$depth,
      count              = node$n,
      divergence         = KL(p, q),
      entropy            = H(p),
      entropy_before     = H(q),
      entropy_drop       = H(q) - H(p),
      likely_next        = modal_next,
      likely_before      = modal_parent,
      changes_prediction = modal_next != modal_parent,
      stringsAsFactors   = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0L) return(.pt_empty_dependence_df())
  out <- out[order(-out[[sort_by]]), , drop = FALSE]
  rownames(out) <- NULL
  if (!is.null(top)) out <- utils::head(out, n = as.integer(top))
  out
}
