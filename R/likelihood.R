# ---- Likelihood, perplexity, and stats S3 generics ----
#
# Implements:
#   logLik.transitiontrees   - in-sample if newdata = NULL, held-out otherwise
#   nobs.transitiontrees
#   perplexity()      - exp(- ll / n)
#   score_sequences() - per-sequence log-lik
#   score_positions() - per-position log-lik

#' @noRd
.pt_unique_counts <- function(tree) {
  ## Each token is attributed to exactly one node (its deepest matching
  ## context); summing unique counts across all nodes equals n_obs.
  children_of <- .pt_children_of(tree)
  lapply(setNames(names(tree$nodes), names(tree$nodes)), function(ctx) {
    info     <- tree$nodes[[ctx]]
    children <- children_of[[ctx]]
    if (is.null(children) || length(children) == 0L) return(info$counts)
    child_sum <- Reduce("+", lapply(children,
                                    function(c) tree$nodes[[c]]$counts))
    info$counts - child_sum
  })
}

#' @noRd
.pt_loglik_in_sample <- function(tree) {
  uc <- .pt_unique_counts(tree)
  contribs <- vapply(names(tree$nodes), function(ctx) {
    counts <- uc[[ctx]]
    p      <- tree$nodes[[ctx]]$prob
    msk    <- counts > 0 & p > 0
    if (!any(msk)) return(0)
    sum(counts[msk] * log(p[msk]))
  }, numeric(1))
  list(ll = sum(contribs), n = sum(vapply(uc, sum, numeric(1))))
}

#' @noRd
.pt_score_walk <- function(tree, newdata) {
  trajs <- .ct_traj(.ct_coerce(newdata))
  if (length(trajs) == 0L)
    .tt_abort("No usable held-out sequences after coercion.")
  alpha <- tree$alphabet

  per_seq <- lapply(seq_along(trajs), function(i) {
    traj <- trajs[[i]]
    L    <- length(traj)
    if (L < 1L) return(NULL)
    ctx_per_pos <- vapply(seq_len(L), function(t) {
      hist <- if (t == 1L) character(0) else traj[seq_len(t - 1L)]
      .ct_match_context(tree, hist)
    }, character(1))
    obs    <- traj
    idx    <- match(obs, alpha)
    keep   <- !is.na(idx)
    if (!any(keep)) return(NULL)
    ctx_k  <- ctx_per_pos[keep]
    pos_k  <- seq_len(L)[keep]
    obs_k  <- obs[keep]
    idx_k  <- idx[keep]
    p_k    <- vapply(seq_along(ctx_k),
                     function(j) tree$nodes[[ctx_k[j]]]$prob[idx_k[j]],
                     numeric(1))
    list(sequence_id     = rep.int(i, length(pos_k)),
         position        = pos_k,
         matched_context = ifelse(ctx_k == .ROOT, .ROOT_LABEL, ctx_k),
         observed        = obs_k,
         predicted_prob  = p_k,
         log_lik         = ifelse(p_k > 0, log(p_k), -Inf))
  })
  per_seq <- per_seq[!vapply(per_seq, is.null, logical(1))]

  if (length(per_seq) == 0L)
    return(data.frame(sequence_id = integer(0), position = integer(0),
                      matched_context = character(0),
                      observed = character(0),
                      predicted_prob = numeric(0),
                      log_lik = numeric(0),
                      stringsAsFactors = FALSE))

  bind <- function(field, fn) fn(lapply(per_seq, `[[`, field))
  data.frame(
    sequence_id     = bind("sequence_id",     function(x) unlist(x, use.names = FALSE)),
    position        = bind("position",        function(x) unlist(x, use.names = FALSE)),
    matched_context = bind("matched_context", function(x) unlist(x, use.names = FALSE)),
    observed        = bind("observed",        function(x) unlist(x, use.names = FALSE)),
    predicted_prob  = bind("predicted_prob",  function(x) unlist(x, use.names = FALSE)),
    log_lik         = bind("log_lik",         function(x) unlist(x, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
}

#' Compute the log-likelihood of a transition tree
#'
#' @description
#' \code{logLik()} computes the log-likelihood of a fitted tree: the sum, over
#' every observed position, of the log probability that the tree assigned to the
#' state that occurred, given the longest context of the tree that matches the
#' preceding history. It is the basis of the information criteria and of
#' perplexity, and allows trees with different depths, pruning or smoothing to
#' be compared on the same data.
#'
#' Without \code{newdata}, the log-likelihood is computed in-sample from the
#' fitted counts. With \code{newdata}, the sequences are scored under the fitted
#' tree, which gives the held-out log-likelihood. A state that is not among the
#' tree's states is not scored, and the position after it is predicted from the
#' root; the reported number of observations counts only scored positions. The
#' number of parameters is the number of contexts times \eqn{k - 1}, for \eqn{k}
#' states.
#'
#' The result is a \code{logLik} object and works with
#' \code{\link[stats]{AIC}()} and \code{\link[stats]{BIC}()}.
#' \code{\link{model_fit}()} returns log-likelihood, AIC, BIC and perplexity
#' together as a table.
#'
#' @param object A \code{transitiontrees} object.
#' @param newdata Sequences to score, in any form accepted by
#'   \code{\link{context_tree}()}. \code{NULL} (default) gives the in-sample
#'   log-likelihood.
#' @param ... Not used.
#'
#' @return An object of class \code{logLik}: the log-likelihood (natural
#'   logarithm), with attributes \code{df} (number of parameters) and
#'   \code{nobs} (number of scored positions).
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{newdata} has an unsupported
#' type or no usable sequence.
#'
#' @seealso \code{\link{model_fit}()}, \code{\link{perplexity}()},
#'   \code{\link{score_sequences}()}, \code{\link{nobs.transitiontrees}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' logLik(tree)
#' AIC(tree)
#' logLik(tree, newdata = engagement)
#' @export
logLik.transitiontrees <- function(object, newdata = NULL, ...) {
  alpha_size <- length(object$alphabet)
  df <- length(object$nodes) * (alpha_size - 1L)
  if (is.null(newdata)) {
    res <- .pt_loglik_in_sample(object)
    val <- res$ll; n <- res$n
  } else {
    pos <- .pt_score_walk(object, newdata)
    val <- sum(pos$log_lik)
    n   <- nrow(pos)
  }
  structure(val, nobs = as.integer(n), df = as.integer(df),
            class = "logLik")
}

#' Count the observations a transition tree was fitted on
#'
#' @description
#' \code{nobs()} returns the number of observed positions in the sequences the
#' tree was fitted on, the number used by \code{\link[stats]{BIC}()}. For a
#' weighted fit it is the weighted total. A held-out log-likelihood,
#' \code{logLik(tree, newdata)}, carries its own count of scored positions.
#'
#' @param object A \code{transitiontrees} object.
#' @param ... Not used.
#'
#' @return A single integer.
#'
#' @seealso \code{\link{logLik.transitiontrees}()},
#'   \code{\link{model_fit}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' nobs(tree)
#' @export
nobs.transitiontrees <- function(object, ...) {
  as.integer(object$n_obs)
}

#' Compute the perplexity of a transition tree
#'
#' @description
#' \code{perplexity()} computes the perplexity of a fitted tree: the exponential
#' of the negative mean log-likelihood per observed position, \eqn{\exp(-\ell /
#' n)}. It reads as an effective number of equally likely next states: a
#' perplexity of \eqn{k} for \eqn{k} states means the tree predicts no better
#' than chance, and a perplexity of 1 means every next state is predicted with
#' certainty. Lower values indicate better prediction.
#'
#' In-sample perplexity (without \code{newdata}) favors deep, lightly smoothed
#' trees that reproduce the fitting data; held-out perplexity, computed on
#' sequences not used for fitting, measures how well the tree generalizes. It
#' returns a single number. \code{\link{tune_tree}()} uses cross-validated
#' perplexity to choose a tree's settings.
#'
#' @param tree A \code{transitiontrees} object.
#' @param newdata Sequences to score, in any form accepted by
#'   \code{\link{context_tree}()}. \code{NULL} (default) gives in-sample
#'   perplexity.
#'
#' @return A single number; \code{NA} when no position could be scored.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree or
#' \code{newdata} has no usable sequence.
#'
#' @seealso \code{\link{model_fit}()}, \code{\link{score_sequences}()},
#'   \code{\link{tune_tree}()}, \code{\link{plot_predictive}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' perplexity(tree)
#' perplexity(tree, newdata = engagement)
#' @export
perplexity <- function(tree, newdata = NULL) {
  .tt_check_class(tree, "transitiontrees", "tree")
  ll <- logLik.transitiontrees(tree, newdata = newdata)
  n  <- attr(ll, "nobs")
  if (n == 0L) return(NA_real_)
  exp(-as.numeric(ll) / n)
}

#' Score each sequence under a fitted tree
#'
#' @description
#' \code{score_sequences()} scores each sequence by how well the fitted tree
#' predicts it: its log-likelihood, the number of positions scored, and its
#' perplexity (\eqn{\exp(-\ell / n)}). Sequences with high perplexity follow
#' pathways the tree finds unlikely, which makes the score useful for finding
#' atypical sequences or for checking how well a tree fitted on one group
#' describes another.
#'
#' Each position is predicted from the longest context of the tree that matches
#' the preceding states, as in \code{\link{predict.transitiontrees}()}; states
#' not among the tree's states are not scored. It returns one row per scored
#' sequence. For the score of each position use \code{\link{score_positions}()};
#' for the most or least typical sequences directly, use
#' \code{\link{mine_sequences}()}.
#'
#' @param tree A \code{transitiontrees} object.
#' @param newdata Sequences to score, in any form accepted by
#'   \code{\link{context_tree}()}.
#'
#' @return A \code{data.frame} with one row per scored sequence and columns
#'   \code{sequence_id} (position of the sequence among the non-empty
#'   sequences of \code{newdata}),
#'   \code{n_scored} (positions scored), \code{log_lik} (natural logarithm)
#'   and \code{perplexity}.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree or
#' \code{newdata} has an unsupported type or no usable sequence.
#'
#' @seealso \code{\link{score_positions}()}, \code{\link{mine_sequences}()},
#'   \code{\link{perplexity}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' head(score_sequences(tree, engagement))
#' @export
score_sequences <- function(tree, newdata) {
  .tt_check_class(tree, "transitiontrees", "tree")
  pos <- .pt_score_walk(tree, newdata)
  if (nrow(pos) == 0L)
    return(data.frame(sequence_id = integer(0), n_scored = integer(0),
                      log_lik = numeric(0), perplexity = numeric(0),
                      stringsAsFactors = FALSE))
  seq_ids <- sort(unique(pos$sequence_id))
  ll  <- vapply(seq_ids,
                function(i) sum(pos$log_lik[pos$sequence_id == i]),
                numeric(1))
  n   <- vapply(seq_ids,
                function(i) sum(pos$sequence_id == i), integer(1))
  pp  <- exp(-ll / n)
  out <- data.frame(sequence_id = seq_ids, n_scored = n,
                    log_lik = ll, perplexity = pp,
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

#' Score each position of each sequence under a fitted tree
#'
#' @description
#' \code{score_positions()} returns, for every position of every sequence, the
#' context of the tree that was used to predict it, the state that occurred, the
#' probability the tree gave that state, and its log-likelihood contribution.
#' Low probabilities mark the moves the tree did not expect, which locates where
#' in a sequence it departs from the typical pathways.
#'
#' With \code{worst}, only the positions with the lowest predicted
#' probability are returned, in increasing order of probability. Summed by
#' sequence, the log-likelihood contributions give
#' \code{\link{score_sequences}()}.
#'
#' @param tree A \code{transitiontrees} object.
#' @param newdata Sequences to score, in any form accepted by
#'   \code{\link{context_tree}()}.
#' @param worst Number of least expected positions to return. \code{NULL}
#'   (default) returns all positions, in sequence order.
#'
#' @return A \code{data.frame} with one row per scored position and columns
#'   \code{sequence_id}, \code{position}, \code{matched_context} (the
#'   context used for the prediction; \code{"(start)"} for the root),
#'   \code{observed} (the state that occurred), \code{predicted_prob} (its
#'   probability) and \code{log_lik} (its natural logarithm).
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree,
#' \code{worst} is not a single positive number, or \code{newdata} has no
#' usable sequence.
#'
#' @seealso \code{\link{score_sequences}()},
#'   \code{\link{predict.transitiontrees}()},
#'   \code{\link{plot_predictive}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' score_positions(tree, engagement, worst = 5)
#' @export
score_positions <- function(tree, newdata, worst = NULL) {
  .tt_check_class(tree, "transitiontrees", "tree")
  if (!is.null(worst) &&
      (!is.numeric(worst) || length(worst) != 1L || is.na(worst) ||
       worst < 1))
    .tt_abort("'worst' must be NULL or a single positive integer.")
  out <- .pt_score_walk(tree, newdata)
  if (!is.null(worst) && nrow(out) > 0L) {
    out <- out[order(out$predicted_prob), , drop = FALSE]
    out <- utils::head(out, n = as.integer(worst))
    rownames(out) <- NULL
  }
  out
}

#' Summarize the fit of a transition tree
#'
#' @description
#' \code{model_fit()} returns the fit statistics of a tree as one row:
#' log-likelihood, number of parameters, number of scored positions, AIC, BIC
#' and perplexity. The row puts trees with different depths, pruning or
#' smoothing on a common footing, and for a group of trees gives one row per
#' group.
#'
#' Without \code{newdata} the statistics are in-sample. With
#' \code{newdata} every statistic is computed on those sequences, with AIC
#' and BIC formed from the held-out log-likelihood and the tree's number of
#' parameters. The number of parameters is the number of contexts times
#' \eqn{k - 1}, for \eqn{k} states (see
#' \code{\link{logLik.transitiontrees}()}).
#'
#' @param tree A \code{transitiontrees} or \code{transitiontrees_group}
#'   object.
#' @param newdata Sequences to evaluate the tree on. \code{NULL} (default)
#'   evaluates it in-sample.
#'
#' @return A \code{data.frame} with columns \code{logLik}, \code{df},
#'   \code{nobs}, \code{AIC}, \code{BIC} and \code{perplexity}: one row for
#'   a tree, or one row per group, with a leading \code{group} column, for a
#'   \code{transitiontrees_group}. When no position can be scored, the
#'   statistics are \code{NA}.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree or
#' group, or \code{newdata} has no usable sequence.
#'
#' @seealso \code{\link{logLik.transitiontrees}()}, \code{\link{perplexity}()},
#'   \code{\link{summary.transitiontrees_group}()}, \code{\link{tune_tree}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' model_fit(tree)
#' model_fit(prune_tree(tree))
#' @export
model_fit <- function(tree, newdata = NULL) {
  if (inherits(tree, "transitiontrees_group")) {
    parts <- lapply(names(tree), function(nm)
      cbind(group = nm, model_fit(tree[[nm]], newdata = newdata),
            stringsAsFactors = FALSE))
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    return(out)
  }
  .tt_check_class(tree, "transitiontrees", "tree")
  ll  <- logLik(tree, newdata = newdata)
  df  <- attr(ll, "df")
  n   <- attr(ll, "nobs")
  llv <- as.numeric(ll)
  ## With zero scored positions the fit scalars are undefined; report
  ## them all as NA rather than a row mixing -Inf (BIC), finite (AIC),
  ## and NA (perplexity).
  if (n == 0L) {
    return(data.frame(
      logLik     = NA_real_, df = df, nobs = 0L,
      AIC        = NA_real_, BIC = NA_real_, perplexity = NA_real_,
      row.names  = NULL
    ))
  }
  data.frame(
    logLik     = llv,
    df         = df,
    nobs       = n,
    AIC        = -2 * llv + 2 * df,
    BIC        = -2 * llv + log(n) * df,
    perplexity = perplexity(tree, newdata = newdata),
    row.names  = NULL
  )
}

#' Count the contexts of a transition tree
#'
#' @description
#' \code{n_nodes()} returns the number of contexts (nodes) in a fitted tree, the
#' root included: the number shown in the printed header and a direct measure of
#' the tree's complexity, for example before and after pruning.
#'
#' @param tree A \code{transitiontrees} or \code{transitiontrees_group}
#'   object.
#'
#' @return A single integer; for a \code{transitiontrees_group}, a named
#'   integer vector with one count per group.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{tree} is not a fitted tree or
#' group.
#'
#' @seealso \code{\link{prune_tree}()}, \code{\link{compare_pruning}()},
#'   \code{\link{summary.transitiontrees_group}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' n_nodes(tree)
#' n_nodes(prune_tree(tree))
#' @export
n_nodes <- function(tree) {
  if (inherits(tree, "transitiontrees_group"))
    return(vapply(tree, function(t) length(t$nodes), integer(1)))
  .tt_check_class(tree, "transitiontrees", "tree")
  length(tree$nodes)
}
