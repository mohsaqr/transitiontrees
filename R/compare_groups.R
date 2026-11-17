# ---- Multi-group comparison ----

#' @noRd
.cg_depth <- function(ctx) {
  if (identical(ctx, .ROOT)) return(0L)
  length(strsplit(ctx, " -> ", fixed = TRUE)[[1L]])
}

#' Shannon entropy in bits of a probability vector.
#' @noRd
.cg_entropy_bits <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0L) return(0)
  -sum(p * log2(p))
}

#' Count-weighted Jensen-Shannon divergence (bits) across the rows of a
#' K x |alphabet| count matrix. Groups with zero total at this context
#' are dropped; returns 0 when fewer than two groups are present.
#' @noRd
.cg_jsd_bits <- function(counts_mat) {
  tot  <- rowSums(counts_mat)
  keep <- tot > 0
  if (sum(keep) < 2L) return(0)
  cm <- counts_mat[keep, , drop = FALSE]
  tt <- tot[keep]
  w  <- tt / sum(tt)
  p  <- cm / tt                                   # row-wise MLE
  m  <- colSums(w * p)                            # mixture
  h_each <- sum(w * apply(p, 1L, .cg_entropy_bits))
  max(.cg_entropy_bits(m) - h_each, 0)
}

#' Expected counts under independence for a contingency table (the
#' margin-product / grand-total). Shared by the usage G^2 and the
#' Pearson-residual difference map so they cannot drift apart.
#' @noRd
.cg_expected <- function(tab) outer(rowSums(tab), colSums(tab)) / sum(tab)

#' G^2 homogeneity statistic for a pathway's prevalence across groups:
#' the K x 2 table [count_c,g ; opportunities_g - count_c,g].
#' @noRd
.cg_usage_g2 <- function(count_c, n_g) {
  tab <- cbind(count_c, pmax(n_g - count_c, 0))
  if (sum(tab) == 0) return(0)
  e  <- .cg_expected(tab)
  nz <- tab > 0 & e > 0
  2 * sum(tab[nz] * log(tab[nz] / e[nz]))
}

#' Per-context usage G^2 for one group assignment. The denominator for a
#' context is the number of windows at \emph{its own depth} (the sum of
#' all same-depth contexts' counts per group), not the total token count
#' --- so prevalence is "share of the depth-d windows" and is comparable
#' across groups regardless of mean sequence length.
#' @noRd
.cg_usage_vec <- function(cc, depths) {
  K        <- nrow(cc[[1L]])
  d_levels <- sort(unique(depths))
  tot_by_d <- lapply(d_levels, function(d)
    Reduce(`+`, lapply(cc[depths == d], rowSums), numeric(K)))
  names(tot_by_d) <- as.character(d_levels)
  vapply(seq_along(cc), function(i)
    .cg_usage_g2(rowSums(cc[[i]]), tot_by_d[[as.character(depths[i])]]),
    numeric(1))
}

#' Per-context K x |alphabet| count matrices for one group assignment.
#' @noRd
.cg_context_counts <- function(trajs_by_grp, grp_names, contexts, depths,
                               alphabet) {
  k      <- length(alphabet)
  uniq_d <- sort(unique(depths))
  ## one count table per group per needed depth (the same routine that
  ## builds the fitted node counts, so observed and permuted agree)
  per_grp <- lapply(grp_names, function(g) {
    tg <- trajs_by_grp[[g]]
    if (is.null(tg) || length(tg) == 0L) tg <- list(character(0))
    setNames(lapply(uniq_d, function(d)
      .ct_count_table(tg, depth = d, alphabet = alphabet)),
      as.character(uniq_d))
  })
  names(per_grp) <- grp_names
  out <- lapply(seq_along(contexts), function(i) {
    ctx <- contexts[i]; d <- as.character(depths[i])
    mat <- t(vapply(grp_names, function(g) {
      cv <- per_grp[[g]][[d]][[ctx]]
      if (is.null(cv)) rep(0, k) else cv
    }, numeric(k)))
    rownames(mat) <- grp_names
    mat
  })
  names(out) <- contexts
  out
}

#' Per-context statistics (JSD, usage G^2) and the two omnibus totals
#' for one group assignment.
#' @noRd
.cg_stats <- function(labels, pooled, grp_names, contexts, depths,
                      alphabet, root_ctx) {
  tbg <- split(pooled, factor(labels, levels = grp_names))
  cc  <- .cg_context_counts(tbg, grp_names, contexts, depths, alphabet)
  jsd <- vapply(cc, .cg_jsd_bits, numeric(1))
  use <- .cg_usage_vec(cc, depths)
  n_c <- vapply(cc, sum, numeric(1))
  list(jsd = jsd, usage = use,
       omni_beh = sum(n_c * jsd), omni_use = sum(use))
}

#' Test where groups of sequences differ, context by context
#'
#' @description
#' \code{compare_groups()} tests where groups of sequences, such as high and low
#' achievers or two phases of a course, differ, context by context and for the
#' trees as a whole. Groups can differ in two distinct ways: they can reach the
#' same history and then continue differently (a behavioral difference), or they
#' can reach some histories more often than others do (a usage difference). Both
#' are tested for every context represented in the groups' trees and for the
#' trees as a whole, so that a difference can be located in specific pathways
#' rather than reported as one number.
#'
#' The behavioral difference at a context is the Jensen-Shannon divergence (Lin,
#' 1991), in bits, between the groups' next-state distributions, each group
#' weighted by its count at that context. The usage difference is a \eqn{G^2}
#' homogeneity statistic on the context's prevalence: its count against the
#' other windows of the same length in each group. The omnibus statistics are
#' the count-weighted sum of the divergences and the sum of the \eqn{G^2}
#' values. All p-values come from one label permutation: sequences are pooled,
#' group labels are shuffled with group sizes fixed, and every statistic is
#' recomputed from raw counts. Per-context p-values are adjusted for false
#' discovery rate with the Benjamini-Hochberg (1995) procedure, separately for
#' each axis. When one person contributes several sequences, a \code{block}
#' identifier restricts the shuffle to within persons, so the null distribution
#' respects that nesting. The result is a
#' \code{transitiontrees_group_comparison} object whose tables are read with
#' \code{as.data.frame()}; \code{summary()} gives a one-row overview and
#' \code{plot()} the largest behavioral differences.
#'
#' \code{compare_groups()} handles two or more groups. For one overall distance
#' between two trees with a permutation test, use \code{\link{compare_trees}()};
#' it measures the symmetric Kullback-Leibler divergence, not the Jensen-Shannon
#' divergence used here. \code{\link{plot_difference}()} maps the per-context
#' differences of two groups next state by next state.
#'
#' @param group A \code{transitiontrees_group} of unweighted trees, from
#'   \code{context_tree(..., group = )}.
#' @param iter Number of label permutations. Default 999.
#' @param min_count Contexts with a total count across groups below this
#'   are dropped from the table. Default 1.
#' @param seed Seed for the permutations, or \code{NULL} to use the
#'   current random-number state. Default 1. A non-\code{NULL} seed sets
#'   the session's random-number state.
#' @param block Block identifiers for a stratified permutation, one per
#'   sequence in group-then-row order. Default \code{NULL}: the block
#'   stored by \code{context_tree(..., block = )} is used when present,
#'   and labels are shuffled freely otherwise.
#'
#' @return An object of class \code{transitiontrees_group_comparison}. Its
#'   tables are returned by \code{as.data.frame()}:
#'   \code{what = "pathways"} gives one row per tested context with
#'   columns \code{pathway}, \code{depth}, \code{count_total}, one
#'   \code{count_<group>} and one \code{modal_<group>} column per group
#'   (the most likely next state), \code{flips} (whether the groups' modal
#'   next states differ), \code{jsd_bits}, \code{jsd_p}, \code{jsd_padj},
#'   \code{usage_g2}, \code{usage_p} and \code{usage_padj} (\code{NA} for
#'   the root, which has no usage test), sorted by \code{jsd_bits};
#'   \code{what = "omnibus"} gives one row per axis; and
#'   \code{what = "distance"} gives one row per pair of groups with their
#'   \code{\link{tree_distance}()}. \code{summary()} returns a one-row
#'   overview and \code{plot()} the largest behavioral differences.
#'
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_bad_input}}{(error) \code{group} is not a
#'     \code{transitiontrees_group}, has fewer than two groups, or
#'     \code{block} does not have one entry per sequence.}
#'   \item{\code{transitiontrees_weighted_unsupported}}{(error) a member
#'     tree was fitted with weights; the permutation recounts the raw
#'     sequences and would drop them.}
#'   \item{\code{transitiontrees_degenerate_permutation}}{(warning) no
#'     block contains more than one group, so the stratified shuffle
#'     cannot move any label.}
#' }
#'
#' @references
#' Lin, J. (1991). Divergence measures based on the Shannon entropy.
#' \emph{IEEE Transactions on Information Theory}, 37(1), 145-151.
#' \doi{10.1109/18.61115}
#'
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
#' rate: A practical and powerful approach to multiple testing.
#' \emph{Journal of the Royal Statistical Society: Series B
#' (Methodological)}, 57(1), 289-300.
#' \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' @seealso \code{\link{summary.transitiontrees_group_comparison}},
#'   \code{\link{as.data.frame.transitiontrees_group_comparison}} and
#'   \code{\link{plot.transitiontrees_group_comparison}} for the result;
#'   \code{\link{compare_trees}()} for the pairwise test;
#'   \code{\link{plot_difference}()} for the next-state difference map.
#'
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' cmp <- compare_groups(grp, iter = 20)
#' cmp
#' summary(cmp)
#' @export
compare_groups <- function(group, iter = 999L, min_count = 1L, seed = 1L,
                           block = NULL) {
  .tt_check_class(group, "transitiontrees_group", "group")
  grp_names <- names(group)
  K <- length(grp_names)
  if (K < 2L) .tt_abort("'group' must contain at least two groups.")
  for (m in group) .pt_assert_unweighted(m, "compare_groups()")
  alphabet <- group[[1L]]$alphabet
  if (!is.null(seed)) set.seed(seed)

  trajs_by_grp <- lapply(group, function(m) .ct_traj(m$data))
  names(trajs_by_grp) <- grp_names
  pooled <- unlist(trajs_by_grp, recursive = FALSE, use.names = FALSE)
  labels <- rep(grp_names, vapply(trajs_by_grp, length, integer(1)))

  ## A block id carried on each member tree (from context_tree(block = ))
  ## is used by default. Filter it through the same survivor index the
  ## fit used (.ct_traj drops empty rows) so it stays aligned to `pooled`.
  if (is.null(block)) {
    mb <- lapply(seq_along(group), function(i) {
      b <- attr(group[[i]], "block")
      if (is.null(b)) NULL else b[attr(trajs_by_grp[[i]], "idx")]
    })
    if (!any(vapply(mb, is.null, logical(1))))
      block <- unlist(mb, use.names = FALSE)
  }

  ## Permutation scheme. Without 'block', labels are shuffled freely
  ## (sequences assumed exchangeable). With 'block' (one id per pooled
  ## sequence, in group-then-row order), the shuffle is stratified
  ## *within* each block, so nested/repeated-measures structure — e.g.
  ## several prompts from one student — is respected and any
  ## between-block confound is held fixed.
  if (is.null(block)) {
    perm_labels <- function() sample(labels)
  } else {
    if (length(block) != length(labels))
      .tt_abort("'block' must have one entry per sequence (length ",
           length(labels), "), aligned to the group-then-row order.")
    block <- as.character(block)
    blk_idx <- split(seq_along(labels), block)
    ## A block only contributes to the null if it holds >1 distinct
    ## label; if none do (e.g. a block id that is unique per sequence),
    ## no label can ever move and the p-values are a degenerate ~1.
    if (!any(vapply(blk_idx,
                    function(idx) length(unique(labels[idx])) > 1L,
                    logical(1))))
      .tt_warn("Stratified permutation is degenerate: no block contains ",
               "more than one group, so labels never move and p-values ",
               "collapse to ~1. Is 'block' unique per sequence?",
               class = "transitiontrees_degenerate_permutation")
    perm_labels <- function() {
      out <- labels
      ## permute the label *values* within each block; safe for
      ## single-sequence blocks (sample() of a length-1 character vector
      ## returns it, unlike sample() of a single index).
      ## Sequential on purpose: one sample() per block in block order keeps the
      ## RNG draw order fixed (the permutation null is pinned draw by draw).
      for (idx in blk_idx) out[idx] <- sample(labels[idx])
      out
    }
  }

  contexts <- unique(unlist(lapply(group, function(m) names(m$nodes)),
                            use.names = FALSE))
  if (!.ROOT %in% contexts) contexts <- c(.ROOT, contexts)
  depths   <- vapply(contexts, .cg_depth, integer(1))

  ## ---- observed ----
  obs_cc <- .cg_context_counts(
    split(pooled, factor(labels, levels = grp_names)),
    grp_names, contexts, depths, alphabet)
  obs <- list(
    jsd   = vapply(obs_cc, .cg_jsd_bits, numeric(1)),
    usage = .cg_usage_vec(obs_cc, depths))
  obs$omni_beh <- sum(vapply(obs_cc, sum, numeric(1)) * obs$jsd)
  obs$omni_use <- sum(obs$usage)

  ## ---- permutation null (label shuffle) ----
  ge_jsd  <- integer(length(contexts))
  ge_use  <- integer(length(contexts))
  ge_obeh <- 0L; ge_ouse <- 0L
  ## Loop kept on purpose: it accumulates exceedance counts in place, so
  ## memory stays O(contexts) instead of O(iter x contexts), and the
  ## permutations are drawn in a fixed RNG order.
  for (b in seq_len(iter)) {
    s <- .cg_stats(perm_labels(), pooled, grp_names, contexts, depths,
                   alphabet, .ROOT)
    ge_jsd  <- ge_jsd  + (s$jsd   >= obs$jsd   - 1e-12)
    ge_use  <- ge_use  + (s$usage >= obs$usage - 1e-12)
    ge_obeh <- ge_obeh + (s$omni_beh >= obs$omni_beh - 1e-12)
    ge_ouse <- ge_ouse + (s$omni_use >= obs$omni_use - 1e-12)
  }
  jsd_p   <- (1 + ge_jsd) / (iter + 1)
  usage_p <- (1 + ge_use) / (iter + 1)

  ## ---- per-group count / modal columns ----
  count_g <- vapply(grp_names, function(g)
    vapply(obs_cc, function(m) sum(m[g, ]), numeric(1)), numeric(length(contexts)))
  modal_g <- vapply(grp_names, function(g)
    vapply(obs_cc, function(m) {
      if (sum(m[g, ]) == 0) NA_character_ else alphabet[which.max(m[g, ])]
    }, character(1)), character(length(contexts)))
  if (length(contexts) == 1L) {                   # vapply drops to vector
    count_g <- matrix(count_g, nrow = 1L)
    modal_g <- matrix(modal_g, nrow = 1L)
  }
  colnames(count_g) <- grp_names
  colnames(modal_g) <- grp_names
  flips <- apply(modal_g, 1L, function(r) length(unique(r[!is.na(r)])) > 1L)

  pathway <- ifelse(contexts == .ROOT, .ROOT_LABEL, contexts)
  tab <- data.frame(pathway = pathway, depth = depths,
                    count_total = rowSums(count_g),
                    stringsAsFactors = FALSE)
  for (g in grp_names) tab[[paste0("count_", g)]] <- count_g[, g]
  for (g in grp_names) tab[[paste0("modal_", g)]] <- modal_g[, g]
  tab$flips      <- flips
  tab$jsd_bits   <- obs$jsd
  tab$jsd_p      <- jsd_p
  tab$jsd_padj   <- stats::p.adjust(jsd_p, method = "BH")
  ## The root has no prevalence test (it occupies 100% of positions by
  ## definition), so its usage statistic is structurally 0/p=1; report it
  ## as NA and keep it out of the usage FDR family so it does not pad the
  ## multiple-testing count.
  is_root        <- contexts == .ROOT
  tab$usage_g2   <- ifelse(is_root, NA_real_, obs$usage)
  tab$usage_p    <- ifelse(is_root, NA_real_, usage_p)
  tab$usage_padj <- NA_real_
  tab$usage_padj[!is_root] <- stats::p.adjust(usage_p[!is_root],
                                              method = "BH")

  tab <- tab[tab$count_total >= min_count, , drop = FALSE]
  tab <- tab[order(-tab$jsd_bits, -tab$usage_g2), , drop = FALSE]
  rownames(tab) <- NULL

  omnibus <- data.frame(
    axis      = c("behavioral", "usage"),
    statistic = c("count-weighted JSD (bits)", "sum G^2"),
    value     = c(obs$omni_beh, obs$omni_use),
    p_value   = c((1 + ge_obeh) / (iter + 1), (1 + ge_ouse) / (iter + 1)),
    stringsAsFactors = FALSE)

  ## ---- pairwise symmetric-KL distance matrix ----
  dm <- matrix(0, K, K, dimnames = list(grp_names, grp_names))
  if (K >= 2L) {
    pr <- utils::combn(K, 2L)
    for (j in seq_len(ncol(pr))) {
      d <- tree_distance(group[[pr[1, j]]], group[[pr[2, j]]])
      dm[pr[1, j], pr[2, j]] <- dm[pr[2, j], pr[1, j]] <- d
    }
  }

  structure(list(
    pathways = tab, omnibus = omnibus, distance_matrix = dm,
    groups = grp_names, iter = iter, seed = seed,
    n_contexts = nrow(tab), stratified = !is.null(block)),
    class = "transitiontrees_group_comparison")
}

#' Print a group comparison
#'
#' @description
#' \code{print()} prints the groups and permutation design, the omnibus
#' behavioral and usage tests, and the \code{n} contexts with the largest
#' behavioral differences.
#'
#' @param x A \code{transitiontrees_group_comparison} object.
#' @param n Number of contexts to print. Default 10.
#' @param digits Rounding digits for the printed tables. Default 3.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{compare_groups}()},
#'   \code{\link{as.data.frame.transitiontrees_group_comparison}} for every
#'   context.
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' print(compare_groups(grp, iter = 20), n = 3)
#' @export
print.transitiontrees_group_comparison <- function(x, n = 10L, digits = 3L, ...) {
  cat(sprintf("<transitiontrees_group_comparison>  %d groups, %d %spermutations\n",
              length(x$groups), x$iter,
              if (isTRUE(x$stratified)) "stratified " else ""))
  cat(sprintf("  groups   : %s\n", paste(x$groups, collapse = ", ")))
  cat(sprintf("  contexts : %d tested\n", x$n_contexts))
  cat("\nomnibus (permutation):\n")
  om <- x$omnibus
  om$value   <- round(om$value, digits)
  om$p_value <- round(om$p_value, digits)
  print.data.frame(om, row.names = FALSE)
  cat("\ntop pathways by behavioral divergence (JSD):\n")
  keep <- c("pathway", "depth", "count_total",
            grep("^modal_", names(x$pathways), value = TRUE),
            "flips", "jsd_bits", "jsd_padj", "usage_padj")
  show <- utils::head(x$pathways[, keep, drop = FALSE], n)
  show$jsd_bits   <- round(show$jsd_bits, digits)
  show$jsd_padj   <- round(show$jsd_padj, digits)
  show$usage_padj <- round(show$usage_padj, digits)
  print.data.frame(show, row.names = FALSE)
  if (nrow(x$pathways) > n)
    cat(sprintf("# ... %d more (use as.data.frame(x) for the full table)\n",
                nrow(x$pathways) - n))
  invisible(x)
}

#' Summarize a group comparison
#'
#' @description
#' \code{summary()} reduces a \code{\link{compare_groups}()} result to one row:
#' the omnibus behavioral and usage tests, and the number of contexts whose
#' adjusted p-value is below \code{alpha} on each axis or whose modal next state
#' differs between groups.
#'
#' @param object A \code{transitiontrees_group_comparison} object.
#' @param alpha False discovery rate cutoff used to count contexts.
#'   Default 0.05.
#' @param ... Ignored.
#' @return An object of class
#'   \code{summary.transitiontrees_group_comparison}: a one-row data frame
#'   with columns \code{groups} (the group names joined by \code{" vs "}),
#'   \code{n_contexts}, \code{iter}, \code{stratified},
#'   \code{behavioral_jsd} and \code{behavioral_p} (omnibus count-weighted
#'   Jensen-Shannon divergence, bits, and its permutation p-value),
#'   \code{usage_g2} and \code{usage_p} (omnibus usage statistic and its
#'   p-value), \code{alpha}, \code{n_behavioral} and \code{n_usage}
#'   (contexts with adjusted p-value below \code{alpha}) and
#'   \code{n_flips} (contexts whose modal next state differs). It prints as
#'   a short verdict.
#' @seealso \code{\link{compare_groups}()},
#'   \code{\link{as.data.frame.transitiontrees_group_comparison}}.
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' summary(compare_groups(grp, iter = 20), alpha = 0.1)
#' @export
summary.transitiontrees_group_comparison <- function(object, alpha = 0.05, ...) {
  .tt_check_class(object, "transitiontrees_group_comparison", "object")
  stopifnot("`alpha` must be a single number in (0, 1)" =
              is.numeric(alpha) && length(alpha) == 1L && !is.na(alpha) &&
              alpha > 0 && alpha < 1)
  p  <- object$pathways
  om <- object$omnibus
  out <- data.frame(
    groups         = paste(object$groups, collapse = " vs "),
    n_contexts     = as.integer(object$n_contexts),
    iter           = as.integer(object$iter),
    stratified     = isTRUE(object$stratified),
    behavioral_jsd = om$value[om$axis == "behavioral"],
    behavioral_p   = om$p_value[om$axis == "behavioral"],
    usage_g2       = om$value[om$axis == "usage"],
    usage_p        = om$p_value[om$axis == "usage"],
    alpha          = alpha,
    ## A pathway with an undefined adjusted p-value (NA) is not counted
    ## as significant; the count is of pathways known to pass.
    n_behavioral   = sum(p$jsd_padj < alpha, na.rm = TRUE),
    n_usage        = sum(p$usage_padj < alpha, na.rm = TRUE),
    n_flips        = sum(p$flips),
    stringsAsFactors = FALSE, row.names = NULL)
  class(out) <- c("summary.transitiontrees_group_comparison", "data.frame")
  out
}

#' Print a summary of a group comparison
#'
#' @description
#' \code{print()} prints the one-row summary of a \code{\link{compare_groups}()}
#' result as a short verdict: the omnibus tests and the counts of contexts that
#' differ.
#'
#' @param x A \code{summary.transitiontrees_group_comparison} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{summary.transitiontrees_group_comparison}}.
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' print(summary(compare_groups(grp, iter = 20)))
#' @export
print.summary.transitiontrees_group_comparison <- function(x, ...) {
  cat(sprintf("Group comparison: %s (%d contexts, %d %spermutations)\n",
              x$groups, x$n_contexts, x$iter,
              if (isTRUE(x$stratified)) "stratified " else ""))
  cat(sprintf("  behavioral omnibus: JSD = %.3f bits, p = %.3f\n",
              x$behavioral_jsd, x$behavioral_p))
  cat(sprintf("  usage omnibus     : G2 = %.2f, p = %.3f\n",
              x$usage_g2, x$usage_p))
  cat(sprintf("  pathways with a behavioral difference (adjusted p < %.2f): %d\n",
              x$alpha, x$n_behavioral))
  cat(sprintf("  pathways with a usage difference      (adjusted p < %.2f): %d\n",
              x$alpha, x$n_usage))
  cat(sprintf("  pathways whose modal next state differs between groups: %d\n",
              x$n_flips))
  invisible(x)
}

#' Convert a group comparison to a data frame
#'
#' @description
#' \code{as.data.frame()} returns one of the three tables of a
#' \code{\link{compare_groups}()} result: the per-context tests, the omnibus
#' tests, or the pairwise distances between groups.
#'
#' @param x A \code{transitiontrees_group_comparison} object.
#' @param row.names,optional Ignored.
#' @param what \code{"pathways"} (default), \code{"omnibus"} or
#'   \code{"distance"}.
#' @param ... Ignored.
#' @return A data frame. For \code{what = "pathways"}, one row per tested
#'   context, with the columns described in \code{\link{compare_groups}()}.
#'   For \code{what = "omnibus"}, one row per axis with columns
#'   \code{axis}, \code{statistic}, \code{value} and \code{p_value}. For
#'   \code{what = "distance"}, one row per unordered pair of groups with
#'   columns \code{group_a}, \code{group_b} and \code{distance} (the
#'   count-weighted symmetric Kullback-Leibler divergence of
#'   \code{\link{tree_distance}()}, nats).
#' @seealso \code{\link{compare_groups}()},
#'   \code{\link{summary.transitiontrees_group_comparison}}.
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' cmp <- compare_groups(grp, iter = 20)
#' as.data.frame(cmp)
#' as.data.frame(cmp, what = "omnibus")
#' as.data.frame(cmp, what = "distance")
#' @export
as.data.frame.transitiontrees_group_comparison <- function(x, row.names = NULL,
                                                       optional = FALSE,
                                                       what = c("pathways",
                                                                "omnibus",
                                                                "distance"),
                                                       ...) {
  what <- match.arg(what)
  switch(what,
    pathways = x$pathways,
    omnibus  = x$omnibus,
    distance = {
      dm <- x$distance_matrix
      g  <- rownames(dm)
      pr <- if (length(g) >= 2L) utils::combn(length(g), 2L) else
        matrix(integer(0), 2L, 0L)
      data.frame(group_a  = g[pr[1L, ]], group_b = g[pr[2L, ]],
                 distance = dm[cbind(pr[1L, ], pr[2L, ])],
                 stringsAsFactors = FALSE)
    })
}

#' Plot the context-level differences between groups
#'
#' @description
#' \code{plot()} draws a \code{\link{compare_groups}()} result in one of two
#' ways. The default lollipop chart ranks the contexts by the Jensen-Shannon
#' divergence between the groups' next-state distributions, with point size for
#' the context's count and color separating contexts whose adjusted p-value is
#' below \code{alpha}; it shows where, given the same history, the groups
#' continue differently. \code{style = "matrix"} draws the pairwise distances
#' between groups as a labeled heatmap.
#'
#' @param x A \code{transitiontrees_group_comparison} object.
#' @param style \code{"divergence"} (default) for the ranked contexts, or
#'   \code{"matrix"} for the between-group distance heatmap.
#' @param top Number of contexts shown by \code{"divergence"}. Default 15.
#' @param alpha False discovery rate cutoff for highlighting. Default 0.05.
#' @param ... Ignored.
#' @return A ggplot object: for \code{"divergence"}, one row per context
#'   with its divergence in bits; for \code{"matrix"}, one tile per pair of
#'   groups with its symmetric Kullback-Leibler distance printed.
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_empty_result}}{(error) for
#'     \code{"divergence"}, no context has a positive divergence.}
#' }
#' @seealso \code{\link{compare_groups}()}; \code{\link{plot_difference}()}
#'   for the next-state differences within each context.
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' cmp <- compare_groups(grp, iter = 20)
#' plot(cmp, alpha = 0.1)
#' plot(cmp, style = "matrix")
#' @export
plot.transitiontrees_group_comparison <- function(x, style = c("divergence",
                                                            "matrix"),
                                               top = 15L, alpha = 0.05, ...) {
  style <- match.arg(style)
  if (style == "matrix") {
    dm <- x$distance_matrix
    df <- expand.grid(a = rownames(dm), b = colnames(dm),
                      stringsAsFactors = FALSE)
    df$dist <- as.vector(dm)
    return(
      ggplot2::ggplot(df, ggplot2::aes(.data$a, .data$b, fill = .data$dist)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data$dist)),
                           size = 3.2) +
        ggplot2::scale_fill_gradient(low = "#deebf7", high = "#08519c",
                                     name = "sym KL") +
        ggplot2::labs(x = NULL, y = NULL,
                      title = "Between-group distance (symmetric KL)") +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(plot.title =
                         ggplot2::element_text(face = "bold")))
  }

  d <- utils::head(x$pathways[x$pathways$jsd_bits > 0, , drop = FALSE], top)
  if (nrow(d) == 0L)
    .tt_abort("No pathways with positive behavioral divergence to plot.", class = "transitiontrees_empty_result")
  d$pathway <- factor(d$pathway, levels = rev(d$pathway))
  d$sig <- ifelse(d$jsd_padj < alpha, sprintf("FDR < %.2f", alpha), "n.s.")
  ggplot2::ggplot(d, ggplot2::aes(x = .data$jsd_bits, y = .data$pathway,
                                  colour = .data$sig)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$jsd_bits,
                                       yend = .data$pathway),
                          colour = "grey75") +
    ggplot2::geom_point(ggplot2::aes(size = .data$count_total)) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(c("#D55E00", "#bbbbbb"),
                               c(sprintf("FDR < %.2f", alpha), "n.s.")),
      name = NULL) +
    ggplot2::scale_size_continuous(name = "count", range = c(2, 7)) +
    ggplot2::labs(
      x = "behavioral divergence between groups (JSD, bits)", y = NULL,
      title = sprintf("How %s differ by context",
                      paste(x$groups, collapse = " vs ")),
      subtitle = "given the same history, how differently do the groups continue?") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      panel.grid.major.y = ggplot2::element_blank())
}
