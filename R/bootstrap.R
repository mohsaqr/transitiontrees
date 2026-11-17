# ---- Bootstrap pathway uncertainty ----
#
# Methodologically built on Saqr et al. (2025, LAK '25):
# extends the bootstrap-and-stability framework from edges to
# variable-depth pathways.
#
# Design rules:
# - The bootstrap operates on raw counts only. No `nmin` filter, no
#   smoothing inside the resampling loop. Pathways tracked are
#   exactly those in the original tree; their values in each
#   resample are read directly from the per-depth count matrices.
# - Stability: p_stability = bootstrap tail probability that the
#   chosen `stat` (default `count`) falls outside
#   [cr[1] x observed, cr[2] x observed], with a +1 correction.
#   stability_rate = 1 - uncorrected instability fraction is retained
#   as a descriptive companion. count = 0 (pathway absent in a
#   resample) automatically lands outside the band; no special
#   "appearance rate" concept needed.
# - Informativeness: informative_rate = fraction of resamples where
#   the pathway's empirical likelihood-ratio statistic G^2 against
#   its parent context exceeds the chi-square critical value at
#   level `alpha_g2` (df = |alphabet| - 1). Parallels the
#   parametric `prune_tree(criterion = "G2")` decision, but reproducibly
#   across resamples.
# - Performance: per-sequence (context, next_state) pair-index
#   vectors are precomputed once per depth. Each bootstrap
#   iteration is one `unlist`+`tabulate` per depth followed by a
#   per-pathway lookup; no `data.frame` construction inside the
#   loop.

#' @noRd
.pt_pathway_depth <- function(p) {
  is_root <- p == .ROOT_LABEL
  out <- integer(length(p))
  if (any(!is_root)) {
    parts <- strsplit(p[!is_root], " -> ", fixed = TRUE)
    out[!is_root] <- lengths(parts)
  }
  out
}

#' @noRd
.pt_precompute_for_boot <- function(trajs, max_depth, alphabet) {
  ## For each (sequence s, depth d) precompute the integer-encoded
  ## (context, next_state) pair vector. Keys are pair_idx =
  ## (ctx_idx - 1) * na + next_idx, where ctx_idx indexes the
  ## global context dictionary at depth d. Bootstrap aggregation is
  ## then a single `unlist(pair[idx])` + `tabulate` per depth —
  ## a transition fast-path lifted from edges to
  ## pathways.
  na <- length(alphabet)
  pair_per_depth     <- vector("list", max_depth + 1L)
  ctx_dict_per_depth <- vector("list", max_depth + 1L)
  for (d in seq.int(0L, max_depth)) {
    if (d == 0L) {
      ctx_per_seq <- lapply(trajs,
                            function(traj) rep(.ROOT, length(traj)))
      nxt_per_seq <- lapply(trajs,
                            function(traj) match(traj, alphabet))
    } else {
      ctx_per_seq <- lapply(trajs, function(traj) {
        L <- length(traj)
        if (L < d + 1L) return(character(0))
        starts <- seq_len(L - d)
        vapply(starts, function(t)
          paste(traj[t:(t + d - 1L)], collapse = " -> "),
          character(1))
      })
      nxt_per_seq <- lapply(trajs, function(traj) {
        L <- length(traj)
        if (L < d + 1L) return(integer(0))
        starts <- seq_len(L - d)
        match(traj[starts + d], alphabet)
      })
    }
    all_ctx <- unique(unlist(ctx_per_seq, use.names = FALSE))
    ctx_dict_per_depth[[d + 1L]] <- all_ctx

    pair_per_seq <- mapply(function(ctx_v, nxt_v) {
      keep <- !is.na(nxt_v) & !is.na(ctx_v)
      if (!any(keep)) return(integer(0))
      ctx_v <- ctx_v[keep]; nxt_v <- nxt_v[keep]
      ctx_idx <- match(ctx_v, all_ctx)
      (ctx_idx - 1L) * na + nxt_v
    }, ctx_per_seq, nxt_per_seq, SIMPLIFY = FALSE)
    pair_per_depth[[d + 1L]] <- pair_per_seq
  }
  list(pair      = pair_per_depth,
       ctx_dict  = ctx_dict_per_depth,
       alphabet  = alphabet,
       max_depth = max_depth,
       na        = na)
}

#' @noRd
.pt_resolve_pathway_layout <- function(pathways, precomp) {
  ## Resolve each original-tree pathway to (depth, ctx_idx in dict)
  ## and its parent's (depth, ctx_idx). Returns NA for the root's
  ## parent. Done once before the bootstrap loop.
  na <- precomp$na
  P <- length(pathways)
  pw_depth        <- integer(P)
  pw_ctx_idx      <- integer(P)
  parent_depth    <- integer(P)
  parent_ctx_idx  <- integer(P)
  for (j in seq_len(P)) {
    p <- pathways[[j]]
    if (identical(p, .ROOT)) {
      pw_depth[j]       <- 0L
      pw_ctx_idx[j]     <- match(.ROOT, precomp$ctx_dict[[1L]])
      parent_depth[j]   <- NA_integer_
      parent_ctx_idx[j] <- NA_integer_
    } else {
      d <- length(strsplit(p, " -> ", fixed = TRUE)[[1L]])
      pw_depth[j]       <- d
      pw_ctx_idx[j]     <- match(p, precomp$ctx_dict[[d + 1L]])
      parent_depth[j]   <- d - 1L
      parent_ctx_idx[j] <- match(.pt_parent_ctx(p),
                                 precomp$ctx_dict[[d]])
    }
  }
  list(pw_depth = pw_depth, pw_ctx_idx = pw_ctx_idx,
       parent_depth = parent_depth,
       parent_ctx_idx = parent_ctx_idx)
}

#' @noRd
.pt_resample_count_matrices <- function(precomp, idx) {
  ## Return a list of n_ctx_d x na count matrices, one per depth.
  na <- precomp$na
  max_depth <- precomp$max_depth
  out <- vector("list", max_depth + 1L)
  for (d in seq.int(0L, max_depth)) {
    n_ctx_d <- length(precomp$ctx_dict[[d + 1L]])
    if (n_ctx_d == 0L) {
      out[[d + 1L]] <- matrix(0L, 0L, na)
      next
    }
    pooled <- unlist(precomp$pair[[d + 1L]][idx],
                     use.names = FALSE)
    flat <- tabulate(pooled, nbins = n_ctx_d * na)
    out[[d + 1L]] <- matrix(flat, n_ctx_d, na, byrow = TRUE)
  }
  out
}

#' @noRd
.pt_pathway_stats_from_counts <- function(cnt_by_depth, layout) {
  ## Empirical pathway stats from per-depth count matrices. No
  ## smoothing. Returns parallel numeric vectors over the layout's
  ## pathways.
  P <- length(layout$pw_depth)
  count_total <- integer(P)
  prob_next   <- rep(NA_real_,    P)
  modal       <- rep(NA_integer_, P)
  KL          <- rep(NA_real_,    P)
  G2          <- rep(NA_real_,    P)
  flips       <- rep(NA_integer_, P)

  for (j in seq_len(P)) {
    cnt_vec <- cnt_by_depth[[layout$pw_depth[j] + 1L]][
      layout$pw_ctx_idx[j], ]
    tot <- sum(cnt_vec)
    count_total[j] <- as.integer(tot)
    if (tot == 0L) next
    m <- which.max(cnt_vec)[1L]
    modal[j]     <- m
    prob_next[j] <- cnt_vec[m] / tot

    pd  <- layout$parent_depth[j]
    pix <- layout$parent_ctx_idx[j]
    if (is.na(pd) || is.na(pix)) next
    pcnt <- cnt_by_depth[[pd + 1L]][pix, ]
    ptot <- sum(pcnt)
    if (ptot == 0L) next
    pprob <- pcnt / ptot
    flips[j] <- as.integer(m != which.max(pcnt)[1L])

    p <- cnt_vec / tot
    msk <- p > 0
    if (any(p > 0 & pprob == 0)) {
      KL[j] <- Inf
    } else {
      KL[j] <- sum(p[msk] * log(p[msk] / pprob[msk], base = 2))
    }

    G2[j] <- .ct_g2(cnt_vec, pprob)
  }
  list(count = count_total, prob_next = prob_next, modal = modal,
       KL = KL, G2 = G2, flips = flips)
}

#' Bootstrap the stability and informativeness of pathways
#'
#' @description
#' \code{bootstrap_pathways()} assesses, by resampling sequences, whether each
#' context of a fitted tree is stable and informative: whether a sample of the
#' same size would reproduce it, and whether its departure from the shorter
#' history is reproducible. The procedure extends the edge-level bootstrap of
#' transition network analysis (Saqr et al., 2025) from single transitions to
#' variable-length pathways, and separates pathways that reflect a regularity of
#' the process from those that reflect only the particular sequences that were
#' sampled.
#'
#' Each iteration resamples whole sequences with replacement (Efron &
#' Tibshirani, 1993) and recounts every pathway of the original tree from raw
#' counts, without smoothing or a minimum count. A pathway is \emph{stable} when
#' the chosen statistic (its count by default, or the probability of its most
#' likely next state, or its divergence from its parent) stays within
#' \code{consistency_range} times the observed value in enough resamples:
#' \code{p_stability} is the share of resamples outside that band, counting the
#' observed sample once (\eqn{(b + 1) / (iter + 1)}), and the pathway is stable
#' when it is below \code{1 - stability_threshold}. A pathway is
#' \emph{informative} when its likelihood-ratio statistic \eqn{G^2} against its
#' parent context exceeds the chi-square critical value at \code{alpha}
#' (\eqn{|A| - 1} degrees of freedom for alphabet \eqn{A}) in at least
#' \code{informative_threshold} of the resamples. Percentile intervals are
#' reported for the count, the next-state probability, the divergence and
#' \eqn{G^2}.
#'
#' The two flags are read together. A stable and informative pathway is
#' reproducible and predicts differently from its shorter history; a stable
#' pathway that is not informative recurs but adds nothing to its parent's
#' prediction; an informative pathway that is not stable is distinctive but
#' carried by a subset of sequences; a pathway that is neither is fragile. The
#' result is a \code{transitiontrees_bootstrap} object, read with
#' \code{as.data.frame()} (one row per pathway) and drawn with \code{plot()}.
#'
#' The bootstrap evaluates the pathways of one tree; to choose the depth use
#' \code{\link{tune_tree}()}, to remove contexts that do not differ from their
#' parent use \code{\link{prune_tree}()}, and to compare groups use
#' \code{\link{compare_groups}()}.
#'
#' @param tree An unweighted \code{transitiontrees} object that carries
#'   its sequences (any tree from \code{\link{context_tree}()}).
#' @param iter Number of bootstrap resamples. Default 1000.
#' @param stat Statistic used for stability: \code{"count"} (default),
#'   \code{"next_probability"} (probability of the most likely next state)
#'   or \code{"divergence"} (divergence from the parent context, bits).
#' @param consistency_range Lower and upper multiples of the observed
#'   value that bound a consistent resample. Default \code{c(0.5, 1.5)}.
#' @param stability_threshold A pathway is stable when
#'   \code{p_stability < 1 - stability_threshold}. Default 0.95.
#' @param informative_threshold A pathway is informative when the share of
#'   resamples with a significant \eqn{G^2} is at least this. Default 0.80.
#' @param alpha Significance level of the \eqn{G^2} test against the
#'   parent context. Default 0.05.
#' @param ci_level Total tail probability of the percentile intervals.
#'   Default 0.05 (95\% intervals).
#' @param seed Seed for the resampling, or \code{NULL} to use the current
#'   random-number state. Default 1. A non-\code{NULL} seed sets the
#'   session's random-number state.
#' @param keep_resamples Keep the per-resample values (default
#'   \code{TRUE}); needed by \code{as.data.frame(x, what = "resamples")}
#'   and \code{\link{plot_pathway_resamples}()}. \code{FALSE} saves memory
#'   and leaves the summary table unchanged.
#' @param progress Show a progress bar. Default \code{FALSE}.
#'
#' @return An object of class \code{transitiontrees_bootstrap}.
#'   \code{as.data.frame(x)} (or \code{summary(x)}) returns one row per
#'   pathway of the original tree, stable and informative pathways first,
#'   with columns \code{pathway}, \code{depth}, \code{count},
#'   \code{likely_next}, \code{next_probability}, \code{divergence} (bits),
#'   \code{changes_prediction}, \code{G2} (\code{NA} for the root),
#'   \code{p_stability}, \code{stability_rate}, \code{stable},
#'   \code{informative_rate}, \code{informative}, \code{flip_consistency}
#'   (share of resamples reproducing \code{changes_prediction}), and for
#'   each of \code{count}, \code{next_probability}, \code{divergence} and
#'   \code{G2} the bootstrap \code{mean_*}, \code{sd_*}, \code{ci_*_lo}
#'   and \code{ci_*_hi}. \code{as.data.frame(x, what = "resamples")}
#'   returns the per-resample values; \code{plot()} draws the \eqn{G^2}
#'   intervals.
#'
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_bad_input}}{(error) \code{tree} is not a
#'     \code{transitiontrees} object, or has no pathways.}
#'   \item{\code{transitiontrees_weighted_unsupported}}{(error) the tree
#'     was fitted with weights, which resampling would drop.}
#'   \item{\code{transitiontrees_no_data}}{(error) the tree does not carry
#'     its sequences.}
#' }
#'
#' @references
#' Saqr, M., López-Pernas, S., Törmänen, T., Kaliisa, R., Misiejuk, K.,
#' & Tikka, S. (2025). Transition network analysis: A novel framework for
#' modeling, visualizing, and identifying the temporal patterns of
#' learners and learning processes. In \emph{Proceedings of the 15th
#' International Learning Analytics and Knowledge Conference (LAK '25)}
#' (pp. 351-361). \doi{10.1145/3706468.3706513}
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An Introduction to the
#' Bootstrap}. Chapman & Hall/CRC. \doi{10.1201/9780429246593}
#'
#' @seealso \code{\link{as.data.frame.transitiontrees_bootstrap}},
#'   \code{\link{plot.transitiontrees_bootstrap}} and
#'   \code{\link{plot_pathway_resamples}()} for the result;
#'   \code{\link{tree_pathways}()} for the observed pathway table;
#'   \code{\link{prune_tree}()} for significance-based pruning.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' boot <- bootstrap_pathways(tree, iter = 50)
#' boot
#' head(as.data.frame(boot))
#' @export
bootstrap_pathways <- function(tree,
                               iter                  = 1000L,
                               stat                  = c("count",
                                                         "next_probability",
                                                         "divergence"),
                               consistency_range     = c(0.5, 1.5),
                               stability_threshold   = 0.95,
                               informative_threshold = 0.80,
                               alpha                 = 0.05,
                               ci_level              = 0.05,
                               seed                  = 1L,
                               keep_resamples        = TRUE,
                               progress              = FALSE) {
  .tt_check_class(tree, "transitiontrees", "tree")
  .pt_assert_unweighted(tree, "bootstrap_pathways()")
  if (is.null(tree$data))
    .tt_abort("tree$data is missing; bootstrap requires the original ",
         "sequences. Refit with context_tree() (>=0.1.1 keeps data).",
              class = "transitiontrees_no_data")
  iter <- as.integer(iter)
  stopifnot(
    "`iter` must be a single integer >= 2" =
      length(iter) == 1L && !is.na(iter) && iter >= 2L,
    "`ci_level` must be a single number in (0, 1)" =
      length(ci_level) == 1L && ci_level > 0 && ci_level < 1,
    "`alpha` must be a single number in (0, 1)" =
      length(alpha) == 1L && alpha > 0 && alpha < 1,
    "`stability_threshold` must be a single number in (0, 1)" =
      length(stability_threshold) == 1L &&
      stability_threshold > 0 && stability_threshold < 1,
    "`informative_threshold` must be a single number in (0, 1)" =
      length(informative_threshold) == 1L &&
      informative_threshold > 0 && informative_threshold < 1,
    "`consistency_range` must be two increasing positive numbers" =
      length(consistency_range) == 2L && all(consistency_range > 0) &&
      consistency_range[1L] < consistency_range[2L]
  )
  stat <- match.arg(stat)
  if (!is.null(seed)) set.seed(as.integer(seed))

  trajs    <- tree$data
  n        <- length(trajs)
  alphabet <- tree$alphabet
  na       <- length(alphabet)
  max_depth <- tree$max_depth
  g2_crit   <- stats::qchisq(1 - alpha, df = na - 1L)

  precomp <- .pt_precompute_for_boot(trajs, max_depth, alphabet)

  ## Pathways from the original tree (raw — no resample dependency).
  orig_path_keys <- names(tree$nodes)
  if (length(orig_path_keys) == 0L)
    .tt_abort("Original tree has no pathways. Did fitting succeed?")
  layout <- .pt_resolve_pathway_layout(orig_path_keys, precomp)
  P <- length(orig_path_keys)
  pathway_label <- ifelse(orig_path_keys == .ROOT, .ROOT_LABEL,
                          orig_path_keys)

  ## Original pathway stats — recompute empirically (no smoothing)
  ## so original and bootstrap values are on the same scale.
  full_idx  <- seq_len(n)
  cnt_full  <- .pt_resample_count_matrices(precomp, full_idx)
  orig_stat <- .pt_pathway_stats_from_counts(cnt_full, layout)
  orig_modal_next <- ifelse(is.na(orig_stat$modal), NA_character_,
                            alphabet[orig_stat$modal])
  orig_flips      <- ifelse(is.na(orig_stat$flips), NA, as.logical(orig_stat$flips))

  ## Bootstrap loop.
  M_count <- matrix(0L,           iter, P,
                    dimnames = list(NULL, pathway_label))
  M_prob  <- matrix(NA_real_,     iter, P,
                    dimnames = list(NULL, pathway_label))
  M_KL    <- matrix(NA_real_,     iter, P,
                    dimnames = list(NULL, pathway_label))
  M_G2    <- matrix(NA_real_,     iter, P,
                    dimnames = list(NULL, pathway_label))
  M_flips <- matrix(NA_integer_,  iter, P,
                    dimnames = list(NULL, pathway_label))

  if (progress) pb <- utils::txtProgressBar(min = 0, max = iter,
                                            style = 3)
  ## Loop kept on purpose: resamples are drawn in a fixed RNG order (the
  ## parity contract pins the draw order) and the progress bar ticks per
  ## resample; vapply() would hold every resample's tables at once.
  for (b in seq_len(iter)) {
    idx <- sample.int(n, n, replace = TRUE)
    cnt_b <- .pt_resample_count_matrices(precomp, idx)
    s     <- .pt_pathway_stats_from_counts(cnt_b, layout)
    M_count[b, ] <- s$count
    M_prob[b, ]  <- s$prob_next
    M_KL[b, ]    <- s$KL
    M_G2[b, ]    <- s$G2
    M_flips[b, ] <- s$flips
    if (progress) utils::setTxtProgressBar(pb, b)
  }
  if (progress) close(pb)

  ## CIs (ci_level = 0.05 -> 95%).
  q_probs    <- c(ci_level / 2, 1 - ci_level / 2)
  col_quants <- function(M) apply(M, 2L, stats::quantile,
                                  probs = q_probs, na.rm = TRUE)
  ci_count <- col_quants(M_count)
  ci_prob  <- col_quants(M_prob)
  ci_KL    <- col_quants(M_KL)
  ci_G2    <- col_quants(M_G2)

  mean_count <- colMeans(M_count)
  sd_count   <- apply(M_count,    2L, stats::sd)
  mean_prob  <- colMeans(M_prob,  na.rm = TRUE)
  sd_prob    <- apply(M_prob,     2L, stats::sd, na.rm = TRUE)
  mean_KL    <- colMeans(M_KL,    na.rm = TRUE)
  sd_KL      <- apply(M_KL,       2L, stats::sd, na.rm = TRUE)
  mean_G2    <- colMeans(M_G2,    na.rm = TRUE)
  sd_G2      <- apply(M_G2,       2L, stats::sd, na.rm = TRUE)

  ## Stability rate on the chosen stat (band test).
  M_stat_arr  <- switch(stat,
                        count            = M_count,
                        next_probability = M_prob,
                        divergence       = M_KL)
  orig_stat_v <- switch(stat,
                        count            = orig_stat$count,
                        next_probability = orig_stat$prob_next,
                        divergence       = orig_stat$KL)
  cr_lo <- pmin(orig_stat_v * consistency_range[1],
                orig_stat_v * consistency_range[2])
  cr_hi <- pmax(orig_stat_v * consistency_range[1],
                orig_stat_v * consistency_range[2])
  n_in <- vapply(seq_len(P), function(j) {
    if (is.na(orig_stat_v[j])) return(NA_integer_)
    sum(M_stat_arr[, j] >= cr_lo[j] & M_stat_arr[, j] <= cr_hi[j],
        na.rm = TRUE)
  }, integer(1L))
  stability_rate <- n_in / iter
  n_out <- iter - n_in
  p_stability <- (n_out + 1) / (iter + 1)
  p_stability[is.na(n_in)] <- NA_real_
  stable <- !is.na(p_stability) &
            p_stability < (1 - stability_threshold)

  ## Informative rate: fraction of resamples where empirical G^2
  ## exceeds the chi-square critical value (df = na - 1).
  n_inf <- vapply(seq_len(P), function(j) {
    g <- M_G2[, j]
    valid <- !is.na(g) & is.finite(g)
    if (!any(valid)) return(NA_integer_)
    sum(g[valid] > g2_crit)
  }, integer(1L))
  informative_rate <- n_inf / iter
  informative_rate[is.na(orig_stat$G2)] <- NA_real_
  informative <- !is.na(informative_rate) &
                 informative_rate >= informative_threshold

  ## Modal-flip consistency vs the original empirical modal-flip
  ## flag. Descriptive companion to stability/informativeness.
  flip_consistency <- vapply(seq_len(P), function(j) {
    of <- orig_stat$flips[j]
    if (is.na(of)) return(NA_real_)
    fb <- M_flips[, j]
    valid <- !is.na(fb)
    if (!any(valid)) return(NA_real_)
    mean(fb[valid] == as.integer(of))
  }, numeric(1L))

  summary_df <- data.frame(
    pathway              = pathway_label,
    depth                = layout$pw_depth,
    count                = orig_stat$count,
    likely_next          = orig_modal_next,
    next_probability     = orig_stat$prob_next,
    divergence           = orig_stat$KL,
    changes_prediction   = orig_flips,
    G2                   = orig_stat$G2,
    p_stability          = p_stability,
    stability_rate       = stability_rate,
    stable               = stable,
    informative_rate     = informative_rate,
    informative          = informative,
    flip_consistency     = flip_consistency,
    mean_count           = mean_count,
    sd_count             = sd_count,
    ci_count_lo          = ci_count[1L, ],
    ci_count_hi          = ci_count[2L, ],
    mean_next_probability   = mean_prob,
    sd_next_probability     = sd_prob,
    ci_next_probability_lo  = ci_prob[1L, ],
    ci_next_probability_hi  = ci_prob[2L, ],
    mean_divergence      = mean_KL,
    sd_divergence        = sd_KL,
    ci_divergence_lo     = ci_KL[1L, ],
    ci_divergence_hi     = ci_KL[2L, ],
    mean_G2              = mean_G2,
    sd_G2                = sd_G2,
    ci_G2_lo             = ci_G2[1L, ],
    ci_G2_hi             = ci_G2[2L, ],
    stringsAsFactors     = FALSE
  )
  ## Sort key: 2 (both flags), 1 (one flag), 0 (neither).
  ## Pathways with both stable + informative come first; ties
  ## break on stability_rate then count.
  trust_score <- as.integer(summary_df$stable) +
                 as.integer(summary_df$informative)
  trust_score[is.na(trust_score)] <- 0L
  ord <- order(-trust_score,
               -ifelse(is.na(summary_df$stability_rate), 0,
                       summary_df$stability_rate),
               -ifelse(is.na(summary_df$count), 0,
                       summary_df$count))
  summary_df <- summary_df[ord, , drop = FALSE]
  rownames(summary_df) <- NULL

  pathways_orig <- data.frame(
    pathway            = pathway_label,
    depth              = layout$pw_depth,
    count              = orig_stat$count,
    likely_next        = orig_modal_next,
    next_probability   = orig_stat$prob_next,
    divergence         = orig_stat$KL,
    changes_prediction = orig_flips,
    G2                 = orig_stat$G2,
    stringsAsFactors   = FALSE
  )

  if (!isTRUE(keep_resamples)) {
    M_count <- NULL; M_prob <- NULL; M_KL <- NULL
    M_G2    <- NULL; M_flips <- NULL
  }

  structure(
    list(
      pathways_orig         = pathways_orig,
      summary               = summary_df,
      M_count               = M_count,
      M_next_probability    = M_prob,
      M_divergence          = M_KL,
      M_G2                  = M_G2,
      M_changes_prediction  = M_flips,
      iter                  = iter,
      stat                  = stat,
      consistency_range     = consistency_range,
      stability_threshold   = stability_threshold,
      informative_threshold = informative_threshold,
      alpha_g2              = alpha,
      ci_level              = ci_level,
      g2_critical_value     = g2_crit,
      seed                  = seed
    ),
    class = "transitiontrees_bootstrap"
  )
}

#' Print a pathway bootstrap
#'
#' @description
#' \code{print()} prints the number of resamples, the stability and
#' informativeness criteria, the number of stable and informative pathways, and
#' the first \code{n} pathways (stable and informative first).
#'
#' @param x A \code{transitiontrees_bootstrap} object.
#' @param n Number of pathways to print. Default 10.
#' @param digits Rounding digits for the printed table. Default 3.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{bootstrap_pathways}()},
#'   \code{\link{as.data.frame.transitiontrees_bootstrap}} for every
#'   pathway.
#' @examples
#' boot <- bootstrap_pathways(context_tree(trajectories, max_depth = 2),
#'                            iter = 50)
#' print(boot, n = 5)
#' @export
print.transitiontrees_bootstrap <- function(x, n = 10L, digits = 3L, ...) {
  cat(sprintf(
    "<transitiontrees_bootstrap>  %d resamples\n", x$iter))
  cat(sprintf(
    "  stability  : %s in [%.2f, %.2f] x observed, p < %.2f\n",
    x$stat, x$consistency_range[1L], x$consistency_range[2L],
    1 - x$stability_threshold))
  cat(sprintf(
    "  informative: G^2 > qchisq(%.2f, df=k-1) = %.2f, threshold %.2f\n",
    1 - x$alpha_g2, x$g2_critical_value,
    x$informative_threshold))
  s <- x$summary
  n_total <- nrow(s)
  n_stable <- sum(s$stable, na.rm = TRUE)
  n_info   <- sum(s$informative, na.rm = TRUE)
  n_both   <- sum(s$stable & s$informative, na.rm = TRUE)
  cat(sprintf(
    "  pathways   : %d total, %d stable, %d informative, %d both\n",
    n_total, n_stable, n_info, n_both))
  cat("\ntop pathways (stable + informative first):\n")
  show <- utils::head(s, n)
  show_disp <- show[, c("pathway", "depth", "count",
                        "p_stability", "stability_rate", "stable",
                        "informative_rate", "informative",
                        "mean_G2", "ci_G2_lo", "ci_G2_hi")]
  num_cols <- vapply(show_disp, is.numeric, logical(1))
  show_disp[, num_cols] <- lapply(show_disp[, num_cols], round,
                                  digits = digits)
  print.data.frame(show_disp, row.names = FALSE)
  if (n_total > n) {
    cat(sprintf(
      "# ... %d more pathways (use summary(x) for full table)\n",
      n_total - n))
  }
  invisible(x)
}

#' Summarize a pathway bootstrap
#'
#' @description
#' \code{summary()} returns the per-pathway table of a
#' \code{\link{bootstrap_pathways}()} result; identical to
#' \code{as.data.frame(object)}.
#'
#' @param object A \code{transitiontrees_bootstrap} object.
#' @param ... Ignored.
#' @return A data frame with one row per pathway of the original tree; see
#'   \code{\link{bootstrap_pathways}()} for the columns.
#' @seealso \code{\link{bootstrap_pathways}()},
#'   \code{\link{as.data.frame.transitiontrees_bootstrap}}.
#' @examples
#' boot <- bootstrap_pathways(context_tree(trajectories, max_depth = 1),
#'                            iter = 50)
#' summary(boot)
#' @export
summary.transitiontrees_bootstrap <- function(object, ...) {
  object$summary
}

#' Convert a pathway bootstrap to a data frame
#'
#' @description
#' \code{as.data.frame()} returns the per-pathway summary of a
#' \code{\link{bootstrap_pathways}()} result, or the values of one statistic in
#' every resample.
#'
#' @param x A \code{transitiontrees_bootstrap} object.
#' @param row.names,optional Ignored.
#' @param what \code{"summary"} (default) for one row per pathway, or
#'   \code{"resamples"} for one row per resample and pathway.
#' @param stat For \code{what = "resamples"}, the statistic:
#'   \code{"count"} (default), \code{"next_probability"},
#'   \code{"divergence"}, \code{"G2"} or \code{"changes_prediction"}.
#' @param ... Ignored.
#' @return A data frame. For \code{what = "summary"}, one row per pathway
#'   with the columns described in \code{\link{bootstrap_pathways}()}. For
#'   \code{what = "resamples"}, one row per resample and pathway with
#'   columns \code{resample}, \code{pathway}, \code{stat} and \code{value}
#'   (\code{NA} where the statistic is undefined in that resample).
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_no_data}}{(error) \code{what =
#'     "resamples"} for a bootstrap run with \code{keep_resamples =
#'     FALSE}.}
#' }
#' @seealso \code{\link{bootstrap_pathways}()},
#'   \code{\link{plot_pathway_resamples}()}.
#' @examples
#' boot <- bootstrap_pathways(context_tree(trajectories, max_depth = 2),
#'                            iter = 50)
#' head(as.data.frame(boot))
#' head(as.data.frame(boot, what = "resamples", stat = "G2"))
#' @export
as.data.frame.transitiontrees_bootstrap <- function(x, row.names = NULL,
                                              optional = FALSE,
                                              what = c("summary", "resamples"),
                                              stat = c("count",
                                                       "next_probability",
                                                       "divergence", "G2",
                                                       "changes_prediction"),
                                              ...) {
  what <- match.arg(what)
  if (identical(what, "summary")) return(x$summary)
  stat <- match.arg(stat)
  M <- x[[paste0("M_", stat)]]
  if (is.null(M))
    .tt_abort("Resample matrices were not retained on this bootstrap ",
              "object. Re-run with keep_resamples = TRUE.",
              class = "transitiontrees_no_data")
  data.frame(resample = rep(seq_len(nrow(M)), times = ncol(M)),
             pathway  = rep(colnames(M), each = nrow(M)),
             stat     = stat,
             value    = as.numeric(M),
             stringsAsFactors = FALSE)
}

#' Plot bootstrap intervals of pathway informativeness
#'
#' @description
#' \code{plot()} draws a forest plot of each pathway's likelihood-ratio
#' statistic \eqn{G^2} against its parent context: the bootstrap mean with its
#' percentile interval, one row per pathway, colored and ordered by whether the
#' pathway is stable, informative, both or neither. A dashed line marks the
#' chi-square critical value used by the bootstrap; a pathway whose interval
#' lies wholly to its right departs reproducibly from its shorter history.
#'
#' @param x A \code{transitiontrees_bootstrap} object.
#' @param top Maximum number of pathways shown. Default 25.
#' @param min_stability Smallest \code{stability_rate} shown. Default
#'   \code{NULL}: the \code{stability_threshold} of the bootstrap.
#' @param ... Ignored.
#' @return A ggplot object: one row per pathway with the mean \eqn{G^2}
#'   and its interval.
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_empty_result}}{(error) no pathway with a
#'     \eqn{G^2} reaches \code{min_stability}.}
#' }
#' @seealso \code{\link{bootstrap_pathways}()},
#'   \code{\link{plot_pathway_resamples}()} for the resample
#'   distributions.
#' @examples
#' boot <- bootstrap_pathways(context_tree(trajectories, max_depth = 2),
#'                            iter = 50)
#' plot(boot, min_stability = 0.5)
#' @export
plot.transitiontrees_bootstrap <- function(x, top = 25L,
                                    min_stability = NULL, ...) {
  if (is.null(min_stability)) min_stability <- x$stability_threshold
  s <- x$summary
  s <- s[!is.na(s$stability_rate) &
         s$stability_rate >= min_stability &
         !is.na(s$mean_G2), , drop = FALSE]
  trust <- as.integer(s$stable) + as.integer(s$informative)
  trust[is.na(trust)] <- 0L
  s <- s[order(-trust, -s$stability_rate, -s$count), ,
         drop = FALSE]
  s <- utils::head(s, top)
  if (nrow(s) == 0L)
    .tt_abort("No pathways meet min_stability = ", min_stability, ".",
              class = "transitiontrees_empty_result")
  s$pathway <- factor(s$pathway, levels = rev(s$pathway))
  s$status <- factor(
    ifelse(s$stable & s$informative,        "stable + informative",
    ifelse(s$stable & !s$informative,        "stable only",
    ifelse(!s$stable & s$informative,        "informative only",
                                              "neither"))),
    levels = c("stable + informative", "stable only",
               "informative only", "neither"))
  ggplot2::ggplot(s,
    ggplot2::aes(x = .data$mean_G2, y = .data$pathway,
                 colour = .data$status)) +
    ggplot2::geom_vline(xintercept = x$g2_critical_value,
                        linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$ci_G2_lo, xmax = .data$ci_G2_hi),
      width = 0.25, orientation = "y", linewidth = 0.6) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(
      values = c(`stable + informative` = "#0072B2",
                 `stable only`          = "#888888",
                 `informative only`     = "#D55E00",
                 `neither`              = "#cccccc"),
      drop = FALSE, name = NULL) +
    ggplot2::labs(
      x = sprintf("bootstrap G^2 (mean and %d%% CI)",
                  round((1 - x$ci_level) * 100)),
      y = NULL,
      title = sprintf("Pathway bootstrap (%d resamples)", x$iter),
      subtitle = sprintf(
        "showing %d / %d pathways with stability_rate >= %.2f  |  dashed line: G^2 critical value at alpha = %.2f",
        nrow(s), nrow(x$summary), min_stability,
        x$alpha_g2)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40",
                                             size = 9),
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank())
}

#' Plot the bootstrap distribution of pathway statistics
#'
#' @description
#' \code{plot_pathway_resamples()} draws, for each selected pathway, a histogram
#' of one statistic across the bootstrap resamples. The spread of a histogram
#' shows how much the statistic depends on which sequences were sampled; a count
#' histogram with mass near zero marks a pathway that a different sample may not
#' contain.
#'
#' @param x A \code{transitiontrees_bootstrap} object fitted with
#'   \code{keep_resamples = TRUE}.
#' @param pathways Pathways to show, in arrow notation. Default
#'   \code{NULL}: the first \code{top} pathways of the summary table.
#' @param stat \code{"count"} (default), \code{"next_probability"},
#'   \code{"divergence"} or \code{"G2"}.
#' @param top Number of pathways when \code{pathways} is \code{NULL}.
#'   Default 6.
#' @param bins Number of histogram bins. Default 30.
#' @return A ggplot object with one histogram panel per pathway.
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_bad_input}}{(error) \code{x} is not a
#'     \code{transitiontrees_bootstrap}; the resamples were not kept; a
#'     requested pathway is not in the bootstrap.}
#' }
#' @seealso \code{\link{bootstrap_pathways}()},
#'   \code{\link{as.data.frame.transitiontrees_bootstrap}} for the values.
#' @examples
#' boot <- bootstrap_pathways(context_tree(trajectories, max_depth = 2),
#'                            iter = 50)
#' plot_pathway_resamples(boot, stat = "count", top = 4, bins = 15)
#' @export
plot_pathway_resamples <- function(x, pathways = NULL,
                                   stat = c("count", "next_probability",
                                            "divergence", "G2"),
                                   top = 6L, bins = 30L) {
  .tt_check_class(x, "transitiontrees_bootstrap", "x")
  stat <- match.arg(stat)
  M <- switch(stat,
              count            = x$M_count,
              next_probability = x$M_next_probability,
              divergence       = x$M_divergence,
              G2               = x$M_G2)
  if (is.null(M))
    .tt_abort("Resample matrices were not retained on this bootstrap ",
         "object. Re-run with keep_resamples = TRUE.",
              class = "transitiontrees_no_data")
  if (is.null(pathways)) {
    pathways <- utils::head(x$summary$pathway, top)
  }
  missing_pw <- setdiff(pathways, colnames(M))
  if (length(missing_pw) > 0L)
    .tt_abort("Unknown pathway(s): ",
         paste(missing_pw, collapse = ", "))
  rows <- lapply(pathways, function(p) {
    data.frame(pathway = p, value = M[, p],
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  df <- df[!is.na(df$value), , drop = FALSE]
  df$pathway <- factor(df$pathway, levels = pathways)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = bins, fill = "#0072B2",
                            colour = "white", alpha = 0.85) +
    ggplot2::facet_wrap(~ pathway, scales = "free_y") +
    ggplot2::labs(
      x = sprintf("bootstrap %s", stat), y = "resamples",
      title = sprintf("Bootstrap distribution of %s per pathway",
                       stat),
      subtitle = sprintf("%d resamples; %d pathways shown",
                         x$iter, length(pathways))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40",
                                             size = 9),
      strip.text    = ggplot2::element_text(face = "bold"))
}
