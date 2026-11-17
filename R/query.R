# ---- query_pathway() / subtree() / pathway_exists(): introspection ----

#' @noRd
.pt_normalise_pathway <- function(pathway) {
  ## Accept either a character vector of states or a single arrow-
  ## notation string. Returns a single canonical " -> " string. The
  ## special token "(start)" (display label), legacy "(root)", or the
  ## sentinel "<root>" all map to .ROOT.
  if (is.null(pathway) || length(pathway) == 0L) return(.ROOT)
  if (length(pathway) == 1L && is.character(pathway)) {
    p <- pathway
    if (identical(p, .ROOT_LABEL) || identical(p, "(root)") ||
        identical(p, .ROOT)) return(.ROOT)
    return(trimws(p))
  }
  paste(as.character(pathway), collapse = " -> ")
}

#' Look up the next-state distribution after a pathway
#'
#' @description
#' \code{query_pathway()} returns the probability distribution that the fitted
#' tree assigns to the next state after a given history, or the probability of
#' one particular next state. For example, it gives the probability of
#' Disengaged after the history Active, Average, read directly from the fitted
#' model.
#'
#' By default the history need not be a context of the tree: as in
#' prediction, the longest suffix of the history that is in the tree is
#' used, falling back to the root when no suffix matches. With
#' \code{exact = TRUE} only a context of the tree is accepted, and other
#' histories return \code{NA}. To find which context matched, use
#' \code{\link{predict.transitiontrees}()} with \code{tidy = TRUE}; to test
#' whether a pathway is a context, use \code{\link{pathway_exists}()}.
#'
#' @param tree A \code{transitiontrees} object.
#' @param pathway The history, oldest state first: a single string in arrow
#'   notation (\code{"A -> B"}) or a character vector (\code{c("A", "B")}).
#'   \code{"(start)"} or an empty vector denotes the root.
#' @param next_state A state whose probability is returned. \code{NULL}
#'   (default) returns the whole distribution.
#' @param exact Require \code{pathway} to be a context of the tree. Default
#'   \code{FALSE} (use the longest matching suffix).
#'
#' @return With \code{next_state}, a single probability. Otherwise a named
#'   numeric vector with one probability per state, in the order of the
#'   tree's states. \code{NA} for a pathway that is not a context when
#'   \code{exact = TRUE}.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree or
#' \code{next_state} is not one of its states.
#'
#' @seealso \code{\link{pathway_exists}()}, \code{\link{subtree}()},
#'   \code{\link{predict.transitiontrees}()}, \code{\link{tree_pathways}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' query_pathway(tree, "Active -> Average")
#' query_pathway(tree, c("Active", "Average"), next_state = "Disengaged")
#' @export
query_pathway <- function(tree, pathway, next_state = NULL,
                          exact = FALSE) {
  .tt_check_class(tree, "transitiontrees", "tree")
  key <- .pt_normalise_pathway(pathway)

  if (key %in% names(tree$nodes)) {
    p <- tree$nodes[[key]]$prob
  } else if (isTRUE(exact)) {
    p <- rep(NA_real_, length(tree$alphabet))
  } else {
    ## Use longest-matching-suffix fall-back.
    parts <- if (identical(key, .ROOT)) character(0) else
      strsplit(key, " -> ", fixed = TRUE)[[1L]]
    ctx <- .ct_match_context(tree, parts)
    p   <- tree$nodes[[ctx]]$prob
  }
  names(p) <- tree$alphabet

  if (is.null(next_state)) return(p)
  idx <- match(as.character(next_state), tree$alphabet)
  if (is.na(idx))
    .tt_abort("'next_state' (", next_state,
         ") is not in the tree's alphabet.")
  p[[idx]]
}

#' Test whether a pathway is a context of the tree
#'
#' @description
#' \code{pathway_exists()} reports whether a history is itself a context (node)
#' of the fitted tree, that is, whether it occurred at least \code{min_count}
#' times and, for a pruned tree, survived pruning. A history that is not a
#' context is still predicted, from its longest suffix in the tree; see
#' \code{\link{query_pathway}()}.
#'
#' @param tree A \code{transitiontrees} object.
#' @param pathway The history, oldest state first: a string in arrow
#'   notation or a character vector.
#'
#' @return A single logical value.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree.
#'
#' @seealso \code{\link{query_pathway}()}, \code{\link{subtree}()},
#'   \code{\link{tree_pathways}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' pathway_exists(tree, "Active -> Disengaged -> Active")
#' pathway_exists(prune_tree(tree), "Active -> Disengaged -> Active")
#' @export
pathway_exists <- function(tree, pathway) {
  .tt_check_class(tree, "transitiontrees", "tree")
  key <- .pt_normalise_pathway(pathway)
  key %in% names(tree$nodes)
}

#' Extract the branch of a tree below one context
#'
#' @description
#' \code{subtree()} returns the part of a fitted tree that lies below a given
#' context: the context itself and every longer context that ends with it, that
#' is, every history that extends it further into the past. The branch shows how
#' the prediction after a given recent history is refined by what happened
#' before it, and can be summarized, converted to a data frame and plotted like
#' a full tree.
#'
#' Contexts keep their full names, and the states, smoothing and other
#' settings are those of the original tree. The printed header names the
#' context at which the branch was cut.
#'
#' @param tree A \code{transitiontrees} object.
#' @param pathway The context at which to cut, oldest state first: a string
#'   in arrow notation or a character vector. It must be a context of the
#'   tree.
#'
#' @return A \code{transitiontrees} object holding \code{pathway} and its
#'   descendants.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree or
#' \code{pathway} is not one of its contexts.
#'
#' @seealso \code{\link{pathway_exists}()}, \code{\link{query_pathway}()},
#'   \code{\link{plot.transitiontrees}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' branch <- subtree(tree, "Disengaged")
#' summary(branch)
#' as.data.frame(branch, sort_by = "depth")
#' @export
subtree <- function(tree, pathway) {
  .tt_check_class(tree, "transitiontrees", "tree")
  key <- .pt_normalise_pathway(pathway)
  if (!key %in% names(tree$nodes))
    .tt_abort("Pathway '", key, "' is not a node in the tree.")

  edges_by_parent <- .pt_children_of(tree)
  collect <- function(ctx) {
    children <- edges_by_parent[[ctx]]
    if (is.null(children) || length(children) == 0L) return(ctx)
    c(ctx, unlist(lapply(children, collect), use.names = FALSE))
  }
  keep <- collect(key)

  new_nodes <- tree$nodes[keep]
  new_edges <- tree$edges[tree$edges$parent %in% keep &
                           tree$edges$child  %in% keep, , drop = FALSE]
  rownames(new_edges) <- NULL

  new_tree <- tree
  new_tree$nodes <- new_nodes
  new_tree$edges <- new_edges
  attr(new_tree, "local_root") <- key
  new_tree
}

#' Count how often each state occurs
#'
#' @description
#' \code{state_frequencies()} tabulates the states of the sequences a tree was
#' fitted on. The counts are those of the root context, which conditions on no
#' history and counts every observed position once; they form the marginal
#' distribution of states, against which every longer context is compared. For a
#' \code{transitiontrees_group}, the table is computed for each group.
#'
#' @param tree A \code{transitiontrees} or a \code{transitiontrees_group}
#'   object.
#'
#' @return A \code{data.frame} with one row per state (per group and state,
#'   with a leading \code{group} column, for a \code{transitiontrees_group})
#'   and columns \code{state}, \code{count} (observed positions, weighted
#'   for a weighted fit) and \code{proportion}, sorted by decreasing count,
#'   with ties in the order of the tree's states.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is neither a fitted
#' tree nor a group of them.
#'
#' @seealso \code{\link{context_tree}()}, \code{\link{tree_pathways}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' state_frequencies(tree)
#' @export
state_frequencies <- function(tree) {
  if (inherits(tree, "transitiontrees_group")) {
    parts <- lapply(names(tree), function(nm)
      cbind(group = nm, state_frequencies(tree[[nm]]),
            stringsAsFactors = FALSE))
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    return(out)
  }
  .tt_check_class(tree, "transitiontrees", "tree")
  counts <- as.numeric(tree$nodes[[.ROOT]]$counts)
  total  <- sum(counts)
  out <- data.frame(
    state      = tree$alphabet,
    count      = counts,
    proportion = if (total > 0) counts / total else rep(NA_real_, length(counts)),
    stringsAsFactors = FALSE)
  out <- out[order(-out$count, seq_len(nrow(out))), , drop = FALSE]
  rownames(out) <- NULL
  out
}
