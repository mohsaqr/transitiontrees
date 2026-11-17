# ---- Pathway-centric API ----

#' @noRd
.pt_empty_pathways_df <- function() {
  data.frame(
    pathway            = character(0),
    depth              = integer(0),
    count              = numeric(0),
    likely_next        = character(0),
    next_probability   = numeric(0),
    divergence         = numeric(0),
    changes_prediction = logical(0),
    stringsAsFactors   = FALSE
  )
}

#' List the pathways of a fitted transition tree
#'
#' @description
#' \code{tree_pathways()} lists the pathways of a fitted tree. A pathway is a
#' context of the tree: the run of states that precedes a point at which the
#' next state is predicted, written from oldest to newest (\code{"A -> B ->
#' C"}). For each pathway the tree records how often it occurred and the
#' distribution of the state that followed. The pathway table therefore shows
#' which histories occur in the data, what tends to follow each of them, and
#' where a longer history changes what is expected next.
#'
#' It returns one row per pathway with its depth and count, its most likely next
#' state and that state's probability, the Kullback-Leibler divergence of its
#' next-state distribution from that of its parent (the pathway without its
#' oldest state), and whether the most likely next state differs from the
#' parent's. The divergence, in bits, measures how much information the
#' additional state of history carries about the next state. The root, written
#' \code{"(start)"}, conditions on no history and holds the overall distribution
#' of states.
#'
#' The ranked views \code{\link{common_pathways}()},
#' \code{\link{divergent_pathways}()} and \code{\link{sharp_pathways}()}
#' return the top rows of this table by count, divergence and next-state
#' probability. \code{\link{tree_dependence}()} adds the entropy of each
#' pathway's prediction, \code{\link{query_pathway}()} looks up one pathway,
#' and \code{\link{bootstrap_pathways}()} assesses how stable each row is
#' under resampling of sequences.
#'
#' @param tree A \code{transitiontrees} object.
#' @param min_count Keep pathways that occurred at least this many times.
#'   Default 1.
#' @param sort_by Column to sort by: \code{"count"} (default),
#'   \code{"divergence"} or \code{"depth"}.
#' @param decreasing Sort in decreasing order. Default \code{TRUE}.
#' @param ... Not used.
#'
#' @return A \code{data.frame} with one row per pathway and columns
#'   \describe{
#'     \item{pathway}{The context, oldest state first; \code{"(start)"} for
#'       the root.}
#'     \item{depth}{Number of states in the pathway.}
#'     \item{count}{Number of times the pathway was followed by a state
#'       (weighted for a weighted fit).}
#'     \item{likely_next}{Most likely next state.}
#'     \item{next_probability}{Its smoothed probability.}
#'     \item{divergence}{Kullback-Leibler divergence, in bits, of the
#'       pathway's next-state distribution from its parent's; \code{NA} for
#'       the root, and \code{Inf} when the pathway gives probability to a
#'       state that the parent gives probability zero.}
#'     \item{changes_prediction}{Whether \code{likely_next} differs from the
#'       parent's most likely next state; \code{NA} for the root.}
#'   }
#'   A tree with no pathway returns a zero-row table with these columns.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree.
#'
#' @seealso \code{\link{common_pathways}()},
#'   \code{\link{divergent_pathways}()}, \code{\link{sharp_pathways}()},
#'   \code{\link{tree_dependence}()}, \code{\link{plot_pathways}()},
#'   \code{\link{bootstrap_pathways}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' tree_pathways(tree)
#' tree_pathways(tree, min_count = 20, sort_by = "divergence")
#' @export
tree_pathways <- function(tree, min_count = 1L,
                              sort_by = c("count", "divergence", "depth"),
                              decreasing = TRUE, ...) {
  .tt_check_class(tree, "transitiontrees", "tree")
  sort_by   <- match.arg(sort_by)
  min_count <- as.integer(min_count)

  rows <- lapply(names(tree$nodes), function(ctx) {
    info <- tree$nodes[[ctx]]
    if (info$n < min_count) return(NULL)

    pathway <- if (identical(ctx, .ROOT)) .ROOT_LABEL else ctx
    parts <- if (identical(ctx, .ROOT)) character(0) else
      strsplit(ctx, " -> ", fixed = TRUE)[[1L]]

    modal_idx  <- which.max(info$prob)
    modal_next <- tree$alphabet[modal_idx]
    prob_next  <- info$prob[modal_idx]

    ## KL vs parent (NA for root)
    if (identical(ctx, .ROOT)) {
      KL    <- NA_real_
      flips <- NA
    } else {
      parent <- if (length(parts) == 1L) .ROOT else
        paste(parts[-1L], collapse = " -> ")
      par_info <- tree$nodes[[parent]]
      if (is.null(par_info)) {
        KL    <- NA_real_
        flips <- NA
      } else {
        p <- info$prob; q <- par_info$prob
        msk <- p > 0
        ## Reported in bits (log base 2) for the pathway table; the
        ## internal .ct_kl() helper used by pruning works in nats. The
        ## two unit conventions are intentional.
        KL <- if (any(p[msk] > 0 & q[msk] == 0)) Inf else
          sum(p[msk] * log(p[msk] / q[msk], base = 2))
        modal_parent <- tree$alphabet[which.max(par_info$prob)]
        flips <- modal_next != modal_parent
      }
    }

    data.frame(
      pathway            = pathway,
      depth              = info$depth,
      count              = info$n,
      likely_next        = modal_next,
      next_probability   = prob_next,
      divergence         = KL,
      changes_prediction = flips,
      stringsAsFactors   = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0L) return(.pt_empty_pathways_df())

  key <- switch(sort_by,
                count      = out$count,
                divergence = ifelse(is.na(out$divergence), -Inf,
                                    out$divergence),
                depth      = out$depth)
  ord <- order(key, decreasing = decreasing)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


# ---- High-level convenience functions -----------------------------------

#' List the most frequent pathways of a tree
#'
#' @description
#' \code{common_pathways()} returns the pathways that occurred most often,
#' optionally at one depth. Frequent pathways are the typical routes through the
#' process: the histories most learners or sessions actually pass through, and
#' the state each of them most often leads to. Frequency says nothing about
#' whether a pathway's prediction differs from that of its shorter suffix; for
#' that, see \code{\link{divergent_pathways}()}.
#'
#' @param tree A \code{transitiontrees} object.
#' @param top Number of pathways to return. Default 10.
#' @param depth Return only pathways of this depth. \code{NULL} (default)
#'   returns all depths.
#' @param min_count Keep pathways that occurred at least this many times.
#'   Default 1.
#'
#' @return A \code{data.frame} with the columns of
#'   \code{\link{tree_pathways}()}, one row per pathway, sorted by
#'   \code{count} in decreasing order, at most \code{top} rows.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree.
#'
#' @seealso \code{\link{tree_pathways}()},
#'   \code{\link{divergent_pathways}()}, \code{\link{sharp_pathways}()},
#'   \code{\link{plot_pathways}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' common_pathways(tree, top = 5)
#' common_pathways(tree, top = 5, depth = 3)
#' @export
common_pathways <- function(tree, top = 10L, depth = NULL, min_count = 1L) {
  .tt_check_class(tree, "transitiontrees", "tree")
  out <- tree_pathways(tree, min_count = min_count, sort_by = "count",
                           decreasing = TRUE)
  if (!is.null(depth))
    out <- out[out$depth == as.integer(depth), , drop = FALSE]
  out <- utils::head(out, n = as.integer(top))
  rownames(out) <- NULL
  out
}

#' List the pathways whose history most changes the prediction
#'
#' @description
#' \code{divergent_pathways()} returns the pathways whose next-state
#' distribution departs most from that of their parent, the same pathway without
#' its oldest state, measured by the Kullback-Leibler divergence in bits. These
#' are the places where remembering one more state of history changes what is
#' expected next, and therefore where the process depends on more than its most
#' recent state. With \code{flips_only = TRUE}, only pathways whose most likely
#' next state differs from the parent's are returned.
#'
#' Divergence is computed from smoothed probabilities and can be large for
#' rare pathways; set \code{min_count} to exclude them, and use
#' \code{\link{bootstrap_pathways}()} to check whether a divergence is stable
#' under resampling. The root, which has no parent, is never returned.
#'
#' @param tree A \code{transitiontrees} object.
#' @param top Number of pathways to return. Default 10.
#' @param min_count Keep pathways that occurred at least this many times.
#'   Default 1.
#' @param flips_only Return only pathways whose most likely next state
#'   differs from the parent's. Default \code{FALSE}.
#'
#' @return A \code{data.frame} with the columns of
#'   \code{\link{tree_pathways}()}, one row per pathway, sorted by
#'   \code{divergence} in decreasing order, at most \code{top} rows.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree.
#'
#' @seealso \code{\link{tree_pathways}()}, \code{\link{tree_dependence}()},
#'   \code{\link{plot_divergence}()}, \code{\link{bootstrap_pathways}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' divergent_pathways(tree, top = 5, min_count = 20)
#' divergent_pathways(tree, flips_only = TRUE)
#' @export
divergent_pathways <- function(tree, top = 10L, min_count = 1L,
                               flips_only = FALSE) {
  .tt_check_class(tree, "transitiontrees", "tree")
  out <- tree_pathways(tree, min_count = min_count,
                           sort_by = "divergence", decreasing = TRUE)
  out <- out[!is.na(out$divergence), , drop = FALSE]
  if (isTRUE(flips_only))
    out <- out[!is.na(out$changes_prediction) & out$changes_prediction,
               , drop = FALSE]
  utils::head(out, n = as.integer(top))
}

#' List the pathways with the most certain next state
#'
#' @description
#' \code{sharp_pathways()} returns the pathways whose most likely next state has
#' the highest probability. A pathway with a next-state probability close to 1
#' is followed almost always by the same state, a near-deterministic step in the
#' process; a probability close to \eqn{1/k} (\eqn{k} states) means the next
#' state is hard to predict from that history. Sharp pathways with small counts
#' are often artefacts of few observations; set \code{min_count} accordingly.
#'
#' @param tree A \code{transitiontrees} object.
#' @param top Number of pathways to return. Default 10.
#' @param min_count Keep pathways that occurred at least this many times.
#'   Default 1.
#'
#' @return A \code{data.frame} with the columns of
#'   \code{\link{tree_pathways}()}, one row per pathway, sorted by
#'   \code{next_probability} in decreasing order, at most \code{top} rows.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree.
#'
#' @seealso \code{\link{tree_pathways}()}, \code{\link{common_pathways}()},
#'   \code{\link{tree_dependence}()} for the entropy of each prediction.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' sharp_pathways(tree, top = 5, min_count = 20)
#' @export
sharp_pathways <- function(tree, top = 10L, min_count = 1L) {
  .tt_check_class(tree, "transitiontrees", "tree")
  out <- tree_pathways(tree, min_count = min_count, sort_by = "count")
  out <- out[order(-out$next_probability), , drop = FALSE]
  rownames(out) <- NULL
  utils::head(out, n = as.integer(top))
}
