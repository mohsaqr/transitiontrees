# ---- Distribution and predictive-diagnostic plots ----

#' Plot the next-state distribution of each context
#'
#' @description
#' \code{plot_distributions()} draws the full next-state distribution
#' \eqn{P(\mathrm{next} \mid \mathrm{context})}{P(next | context)} of selected
#' contexts as small multiples, one bar chart per context with the most likely
#' next state in a darker shade. The display shows, for each history, whether
#' the next state is close to determined (one tall bar) or open (bars of similar
#' height), and which alternatives remain. By default the most frequent contexts
#' are shown.
#'
#' The same probabilities are shown for many contexts at once as a heatmap
#' by \code{\link{plot_pathways}()}, and along one history's chain of
#' shorter suffixes by \code{\link{plot_pruning}()}.
#'
#' @param tree A \code{transitiontrees} object.
#' @param contexts Contexts to show, in arrow notation. Default
#'   \code{NULL}: the \code{top} most frequent contexts.
#' @param top Number of contexts when \code{contexts} is \code{NULL}.
#'   Default 12.
#' @param min_count Contexts with a lower count are dropped. Default 1.
#'
#' @return A ggplot object with one bar-chart panel per context.
#'
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_bad_input}}{(error) \code{tree} is not a
#'     \code{transitiontrees} object.}
#'   \item{\code{transitiontrees_empty_result}}{(error) no context is
#'     left to plot.}
#' }
#'
#' @seealso \code{\link{plot_pathways}()}, \code{\link{plot_pruning}()};
#'   \code{\link{tree_pathways}()} for the table of contexts.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' plot_distributions(tree, top = 6)
#' plot_distributions(tree, contexts = c("Active", "Active -> Disengaged"))
#' @export
plot_distributions <- function(tree, contexts = NULL, top = 12L,
                               min_count = 1L) {
  .tt_check_class(tree, "transitiontrees", "tree")
  alpha <- tree$alphabet
  pw <- tree_pathways(tree, min_count = min_count, sort_by = "count")
  pw <- if (is.null(contexts))
    utils::head(pw, as.integer(top)) else
    pw[pw$pathway %in% as.character(contexts), , drop = FALSE]
  if (nrow(pw) == 0L)
    .tt_abort("No contexts to plot.", class = "transitiontrees_empty_result")

  long <- do.call(rbind, lapply(pw$pathway, function(p) {
    ctx  <- if (identical(p, .ROOT_LABEL)) .ROOT else p
    info <- tree$nodes[[ctx]]
    data.frame(pathway  = p,
               state    = factor(alpha, levels = alpha),
               prob     = info$prob,
               is_modal = seq_along(alpha) == which.max(info$prob),
               stringsAsFactors = FALSE)
  }))
  long$pathway <- factor(long$pathway, levels = pw$pathway)

  ggplot2::ggplot(long, ggplot2::aes(x = .data$state, y = .data$prob,
                                     fill = .data$is_modal)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::facet_wrap(~ pathway) +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "#08519c", `FALSE` = "#9ecae1"),
      guide = "none") +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "next state", y = "P(next | context)",
      title = "Per-context next-state distributions",
      subtitle = sprintf("%d contexts; modal next state highlighted",
                         nrow(pw))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text    = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank())
}

#' Plot how well a context tree predicts new sequences
#'
#' @description
#' \code{plot_predictive()} plots how much probability a context tree gives to
#' what happens next in sequences it has not seen, the measure of its value as a
#' predictive model. It scores every transition of held-out sequences with
#' \code{\link{score_positions}()} and shows where the model is confident and
#' where it is surprised.
#'
#' \code{type = "position"} plots the probability given to each observed
#' next state against its position in the sequence. \code{type = "ecdf"}
#' plots the empirical cumulative distribution of those probabilities,
#' showing how often the model gives high or low probability to what
#' actually happened. \code{type = "logloss"} plots the per-transition
#' log-loss \eqn{-\log_2 P(\mathrm{observed})}{-log2 P(observed)} in bits against position,
#' with the mean at each position and the overall mean; 0 bits is a
#' certain and correct prediction.
#'
#' For one number summarizing held-out prediction use
#' \code{\link{perplexity}()}, and to choose the settings that minimize it
#' use \code{\link{tune_tree}()}.
#'
#' @param tree A \code{transitiontrees} object.
#' @param newdata Held-out sequences in any format accepted by
#'   \code{\link{context_tree}()}.
#' @param type \code{"position"} (default), \code{"ecdf"} or
#'   \code{"logloss"}.
#'
#' @return A ggplot object of the scored transitions.
#'
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_bad_input}}{(error) \code{tree} is not a
#'     \code{transitiontrees} object.}
#'   \item{\code{transitiontrees_empty_result}}{(error) \code{newdata} has
#'     no transitions to score.}
#' }
#'
#' @seealso \code{\link{score_positions}()} for the scored transitions as
#'   a table; \code{\link{perplexity}()}; \code{\link{tune_tree}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' held_out <- head(engagement, 40)
#' plot_predictive(tree, held_out)
#' plot_predictive(tree, held_out, type = "ecdf")
#' plot_predictive(tree, held_out, type = "logloss")
#' @export
plot_predictive <- function(tree, newdata,
                            type = c("position", "ecdf", "logloss")) {
  .tt_check_class(tree, "transitiontrees", "tree")
  type <- match.arg(type)
  sp <- score_positions(tree, newdata)
  if (nrow(sp) == 0L)
    .tt_abort("No held-out positions to score.", class = "transitiontrees_empty_result")

  base_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      panel.grid.minor = ggplot2::element_blank())

  if (type == "position") {
    sp$sequence_id <- factor(sp$sequence_id)
    return(
      ggplot2::ggplot(sp, ggplot2::aes(x = .data$position,
                                       y = .data$predicted_prob)) +
        ggplot2::geom_line(ggplot2::aes(group = .data$sequence_id),
                           colour = "grey70", alpha = 0.5) +
        ggplot2::geom_point(ggplot2::aes(colour = .data$predicted_prob),
                            size = 2) +
        ggplot2::scale_colour_gradient(low = "#D55E00", high = "#0072B2",
                                       limits = c(0, 1), name = "P(obs)") +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(
          x = "position in sequence",
          y = "predicted probability of observed move",
          title = "Predictive confidence by position",
          subtitle = "low (orange) = the model was surprised here") +
        base_theme)
  }

  if (type == "logloss") {
    ## -log2 P(observed): per-position log-loss.
    sp$logloss <- -sp$log_lik / log(2)
    mean_ll <- mean(sp$logloss)
    ## per-position mean (how surprise evolves as sequences unfold) — a
    ## real trend, unlike connecting points across unrelated sequences.
    pos_mean <- stats::aggregate(logloss ~ position, data = sp, FUN = mean)
    return(
      ggplot2::ggplot(sp, ggplot2::aes(x = .data$position,
                                       y = .data$logloss)) +
        ggplot2::geom_jitter(ggplot2::aes(colour = .data$logloss),
                             width = 0.18, height = 0, size = 1.6,
                             alpha = 0.55) +
        ggplot2::geom_line(data = pos_mean,
                           ggplot2::aes(x = .data$position,
                                        y = .data$logloss),
                           colour = "grey20", linewidth = 0.9) +
        ggplot2::geom_point(data = pos_mean,
                            ggplot2::aes(x = .data$position,
                                         y = .data$logloss),
                            colour = "grey20", size = 1.8) +
        ggplot2::geom_hline(yintercept = mean_ll, linetype = "dashed",
                            colour = "#D55E00") +
        ggplot2::annotate("text", x = max(sp$position), y = mean_ll,
                          label = sprintf(" mean %.2f bits", mean_ll),
                          hjust = 1, vjust = -0.5, size = 3,
                          colour = "#D55E00") +
        ggplot2::scale_colour_gradient(low = "#0072B2", high = "#D55E00",
                                       name = "log-loss\n(bits)") +
        ggplot2::labs(
          x = "position in sequence",
          y = expression(-log[2] ~ P(observed ~ move)),
          title = "Predictive log-loss by position",
          subtitle = paste0("lower = better; points = positions, ",
                            "black line = per-position mean, ",
                            "dashed = overall mean")) +
        base_theme)
  }

  ggplot2::ggplot(sp, ggplot2::aes(x = .data$predicted_prob)) +
    ggplot2::stat_ecdf(geom = "step", colour = "#0072B2", linewidth = 1) +
    ggplot2::xlim(0, 1) +
    ggplot2::labs(
      x = "predicted probability of observed move",
      y = "cumulative fraction of positions",
      title = "Predictive probability ECDF",
      subtitle = sprintf("%d held-out positions scored", nrow(sp))) +
    base_theme
}

# ---- Pruning illustration (suffix-chain view) -----------------------

#' Suffix-chain pruning data for one pathway (internal)
#'
#' Walks the suffix chain of \code{pathway} -- the full context, then the
#' context with its oldest (leftmost) move dropped, and so on down to the
#' root -- and, at each link, evaluates the \code{G2} gain of a context
#' against its one-shorter parent under the same rule as
#' \code{\link{prune_tree}} (criterion \code{"G2"}). The requested pathway
#' must itself be a fitted context; the function errors rather than
#' silently starting from a shorter suffix. Each context is classified into
#' one of three states: \code{"informative"} (its own \code{G2} clears the
#' cutoff, so it adds memory over its parent), \code{"pruned"} (it and every
#' deeper context fail to diverge -- the redundant tail), or \code{"retained"}
#' (its own \code{G2} is below the cutoff but a deeper context diverges, so
#' \code{prune_tree()} keeps it only as a structural bridge).
#'
#' @return A long \code{data.frame}, one row per (context, state):
#'   \code{L} (memory length), \code{context}, \code{label} (the
#'   most-recent move), \code{state}, \code{prob}, \code{count},
#'   \code{g2}, \code{diverges} (own G2 clears the cutoff), \code{pruned}
#'   (redundant tail), \code{retained} (kept only for a deeper context),
#'   \code{status} (\code{"informative"} / \code{"retained"} /
#'   \code{"pruned"} / \code{"root"}).
#' @noRd
.suffix_chain <- function(tree, pathway, alpha = 0.05) {
  .tt_check_class(tree, "transitiontrees", "tree")
  stopifnot("`pathway` must be a single character string" =
              is.character(pathway) && length(pathway) == 1L)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1)
    .tt_abort("'alpha' must be a single number in (0, 1).")
  states <- tree$alphabet
  k      <- length(states)
  g2_cv  <- stats::qchisq(1 - alpha, df = k - 1L)

  toks <- if (pathway %in% c(.ROOT_LABEL, .ROOT)) character(0) else
    strsplit(pathway, " -> ", fixed = TRUE)[[1L]]
  n_tok <- length(toks)

  ## Contexts by memory length L = n_tok .. 1 (deepest first), then root.
  keys <- character(0); Ls <- integer(0)
  if (n_tok >= 1L) {
    Ls   <- seq.int(n_tok, 1L)
    keys <- vapply(Ls, function(d)
      paste(toks[(n_tok - d + 1L):n_tok], collapse = " -> "), character(1))
  }
  present <- keys %in% names(tree$nodes)
  ## The requested pathway must itself be a fitted context. If only a shorter
  ## suffix exists, the chain would silently start below the requested depth
  ## and the plot would answer a different memory question than the one asked,
  ## so refuse rather than quietly shorten.
  if (n_tok >= 1L && !present[1L]) {
    deepest <- if (any(present)) keys[which(present)[1L]] else .ROOT_LABEL
    .tt_abort(sprintf(
      "Context \"%s\" is not in the tree (max_depth = %d%s). The deepest fitted suffix is \"%s\"; call plot_pruning() with a pathway whose full context exists.",
      pathway, tree$max_depth,
      if (isTRUE(tree$pruned)) ", pruned" else "", deepest))
  }
  ## Suffix property: if the deepest context exists, so do all its suffixes.
  ## Guard against any gap regardless, so the chain is always contiguous from
  ## the requested pathway down to the root.
  if (n_tok >= 1L && !all(present)) {
    missing <- keys[!present]
    .tt_abort(sprintf(
      "Suffix chain of \"%s\" is not contiguous: context(s) %s are missing from the tree.",
      pathway, paste(sprintf('"%s"', missing), collapse = ", ")))
  }
  keys <- c(keys, .ROOT)
  Ls   <- c(Ls, 0L)

  nodes <- lapply(keys, function(kk) tree$nodes[[kk]])
  m     <- length(keys)

  ## G2 of each context vs its parent (the next, one-shorter context).
  g2 <- c(vapply(seq_len(m - 1L), function(i)
    .ct_g2(nodes[[i]]$counts, nodes[[i + 1L]]$prob), numeric(1)),
    NA_real_)
  ## A context is INFORMATIVE if its own G2 clears the cutoff (it predicts
  ## differently from its one-shorter parent). It is PRUNED if it, and every
  ## deeper context in the chain, fail to diverge (the redundant tail). The
  ## remainder are RETAINED: their own G2 is below the cutoff, but a deeper
  ## context diverges, so prune_tree() keeps them structurally as a bridge to
  ## that descendant. Only INFORMATIVE contexts genuinely add memory.
  diverges    <- g2 > g2_cv                          # NA for the root
  pruned_tail <- c(cumprod(!diverges[-m]) == 1, NA)  # contiguous tail from leaf
  status <- ifelse(is.na(diverges), "root",
             ifelse(diverges, "informative",
              ifelse(pruned_tail, "pruned", "retained")))
  pruned   <- status == "pruned"                     # back-compat logical
  retained <- status == "retained"
  pruned[is.na(diverges)]   <- NA
  retained[is.na(diverges)] <- NA
  labels <- ifelse(keys == .ROOT, .ROOT_LABEL,
                   vapply(keys, function(kk)
                     utils::tail(strsplit(kk, " -> ", fixed = TRUE)[[1L]],
                                 1L), character(1)))
  display <- ifelse(keys == .ROOT, .ROOT_LABEL, keys)

  rows <- lapply(seq_len(m), function(i)
    data.frame(L        = Ls[i],
               context  = display[i],
               label    = labels[i],
               state    = factor(states, levels = states),
               prob     = as.numeric(nodes[[i]]$prob),
               count    = sum(nodes[[i]]$counts),
               g2       = g2[i],
               diverges = diverges[i],
               pruned   = pruned[i],
               retained = retained[i],
               status   = status[i],
               stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

#' Plot how much memory one pathway needs
#'
#' @description
#' \code{plot_pruning()} plots how much memory one pathway needs: how many past
#' states are needed to predict what follows it. A context tree keeps a longer
#' history only when it predicts the next state differently from its shorter
#' suffix, and the plot makes that decision visible for one pathway. It follows
#' the pathway's suffix chain (the full history, then the history without its
#' earliest state, and so on down to the root) and draws the next-state
#' distribution at every link as its own panel, deepest history on the left and
#' the root on the right.
#'
#' Each context is classified by the \eqn{G^2} rule of
#' \code{\link{prune_tree}()}: the likelihood-ratio statistic of its
#' next-state counts against its one-shorter parent is compared with the
#' chi-square critical value at \code{alpha}. A context whose own
#' \eqn{G^2} exceeds it \emph{adds memory} and is drawn solid; a context
#' below it that is kept only because a deeper context exceeds it is drawn
#' at medium opacity; a context below it with no deeper context exceeding
#' it is \emph{pruned} and drawn faded. Each panel title gives the
#' context, the decision and \eqn{G^2}. The distributions are the node
#' values reported by \code{\link{tree_pathways}()} and
#' \code{\link{query_pathway}()}.
#'
#' Pass an unpruned tree so that the whole chain is visible. To prune the
#' whole tree use \code{\link{prune_tree}()}; to compare pruning rules use
#' \code{\link{compare_pruning}()}.
#'
#' @param tree A \code{transitiontrees} object, usually unpruned.
#' @param pathway One fitted context in arrow notation
#'   (\code{"A -> B -> C"}, earliest state on the left).
#' @param alpha Significance level of the \eqn{G^2} test. Default 0.05.
#'
#' @return A ggplot object with one panel per context of the suffix
#'   chain.
#'
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_bad_input}}{(error) \code{tree} is not a
#'     \code{transitiontrees} object; \code{alpha} is not in (0, 1);
#'     \code{pathway} is not a context of the tree.}
#' }
#'
#' @seealso \code{\link{prune_tree}()}, \code{\link{compare_pruning}()},
#'   \code{\link{plot_distributions}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' plot_pruning(tree, "Average -> Active -> Active")
#' @export
plot_pruning <- function(tree, pathway, alpha = 0.05) {
  .tt_check_class(tree, "transitiontrees", "tree")
  long <- .suffix_chain(tree, pathway, alpha = alpha)
  states <- tree$alphabet
  pal <- .okabe_ito(length(states))

  ## One panel per context along the chain, deepest (full memory) on the
  ## left, root on the right. The decision is encoded by opacity so it never
  ## competes with the next-state fill colours; the full context, the
  ## decision, and G2 live in the panel title. Three states are kept
  ## distinct: informative contexts add memory over their parent; retained
  ## contexts do not, but survive because a deeper context does; pruned
  ## contexts are the redundant tail.
  node <- long[!duplicated(long$context), , drop = FALSE]
  ord  <- node$context[order(-node$L)]
  tag  <- c(informative = "adds memory",
            retained    = "kept for a deeper context",
            pruned      = "pruned (redundant)")
  mk_lab <- function(ctx, L, status, g2)
    if (status == "root")
      sprintf("L=%d   %s\n(root  -  no parent)", L, ctx)
    else
      sprintf("L=%d   %s\n%s   -   G2 = %.1f", L, ctx, tag[[status]], g2)
  labmap <- stats::setNames(
    mapply(mk_lab, node$context, node$L, node$status, node$g2), node$context)
  long$panel <- factor(labmap[long$context], levels = labmap[ord])

  ## Highlight the modal next state in each context. Opacity encodes the
  ## decision: informative solid, retained mid, pruned faded.
  long$modal <- stats::ave(long$prob, long$context,
                           FUN = function(p) seq_along(p) == which.max(p)) == 1
  op_map <- c(informative = 1, retained = 0.7, pruned = 0.45, root = 1)
  long$opacity <- unname(op_map[long$status])

  ggplot2::ggplot(long, ggplot2::aes(x = .data$state, y = .data$prob,
                                     fill = .data$state)) +
    ggplot2::geom_col(ggplot2::aes(alpha = .data$opacity,
                                   colour = .data$modal),
                      width = 0.85, linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$prob)),
                       vjust = -0.4, size = 2.7, colour = "grey25") +
    ggplot2::facet_wrap(~ .data$panel, nrow = 1) +
    ggplot2::scale_fill_manual(values = pal, name = "next state") +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "grey20",
                                            `FALSE` = NA), guide = "none") +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_y_continuous(limits = c(0, 1.08),
                                breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(
      x = NULL, y = "P(next move | context)",
      title = sprintf("How much memory does \"%s\" need?", pathway),
      subtitle = paste0("each panel = one context along the suffix chain ",
                        "(deepest on the left, root on the right); solid = adds memory, ",
                        "mid = kept for a deeper context, faded = pruned")) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      strip.text    = ggplot2::element_text(size = 8.5, lineheight = 1.1),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position = "right")
}

#' Okabe-Ito colorblind-safe qualitative palette (internal)
#' @noRd
.okabe_ito <- function(n) {
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
           "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")
  if (n <= length(pal)) pal[seq_len(n)] else
    grDevices::hcl.colors(n, "Dark 3")
}

#' Map the next-state differences between two groups
#'
#' @description
#' \code{plot_difference()} maps the next-state differences between two groups
#' of a \code{transitiontrees_group}: where the groups share a history and still
#' continue from it differently. The map has one row per context the groups
#' share and one column per next state. With \code{depth = 1} the map is the
#' difference between the two first-order transition matrices; deeper contexts
#' show whether the difference depends on a longer history.
#'
#' By default each cell is the Pearson residual of the first group's count
#' of that next state against the count expected if both groups continued
#' alike from that context, \eqn{(O - E) / \sqrt{E}}{(O - E) / sqrt(E)}, with \eqn{E} from
#' the margins of the context's two-group table. Residuals take the
#' context's count into account, and an absolute value above about 2 is
#' notable. \code{measure = "probability"} shows the plain difference of
#' next-state probabilities instead. Rows are ordered by their total
#' absolute difference, largest at the top. With a
#' \code{\link{compare_groups}()} result as \code{comparison}, contexts
#' with a significant behavioral difference are starred.
#'
#' \code{layout = "tree"} draws a different quantity: the context tree of
#' the pooled sequences, each node and branch colored by which group
#' reaches that context more often, as a share of the group's
#' observations. The map describes differences; the permutation tests are
#' in \code{\link{compare_groups}()}.
#'
#' @param group A \code{transitiontrees_group} object.
#' @param groups The two groups to compare, first minus second. Default
#'   \code{NULL}: the two groups of a two-group object; required for more
#'   groups.
#' @param depth Only contexts of this depth. Default \code{NULL}: all
#'   shared contexts.
#' @param min_count Contexts with a count below this in either group are
#'   dropped. Default 1.
#' @param comparison A \code{transitiontrees_group_comparison} from
#'   \code{\link{compare_groups}()}, used to star contexts with
#'   \code{jsd_padj < alpha}. Default \code{NULL}.
#' @param alpha False discovery rate cutoff for the stars. Default 0.05.
#' @param annotate Print the value in each cell (\code{"tile"} only).
#'   Default \code{TRUE}.
#' @param layout \code{"tile"} (default) for the heatmap, or \code{"tree"}
#'   for the pooled tree colored by prevalence difference.
#' @param measure For \code{"tile"}: \code{"residual"} (default) or
#'   \code{"probability"}.
#'
#' @return A ggplot object: a heatmap of contexts by next states, or a
#'   tree colored by prevalence difference.
#'
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_bad_input}}{(error) \code{group} is not a
#'     \code{transitiontrees_group}; \code{groups} is missing for more than
#'     two groups or does not name two of them; no context is shared at
#'     \code{depth} and \code{min_count}; \code{comparison} is not a
#'     \code{transitiontrees_group_comparison}.}
#'   \item{\code{transitiontrees_empty_result}}{(error) for
#'     \code{layout = "tree"}, the pooled tree is empty.}
#' }
#'
#' @seealso \code{\link{compare_groups}()} for the tests;
#'   \code{\link{plot.transitiontrees_group}} for each group's tree.
#'
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' plot_difference(grp)
#' plot_difference(grp, measure = "probability")
#' plot_difference(grp, comparison = compare_groups(grp, iter = 20),
#'                 alpha = 0.1)
#' plot_difference(grp, layout = "tree")
#' @export
plot_difference <- function(group, groups = NULL, depth = NULL,
                            min_count = 1L, comparison = NULL,
                            alpha = 0.05, annotate = TRUE,
                            layout = c("tile", "tree"),
                            measure = c("residual", "probability")) {
  .tt_check_class(group, "transitiontrees_group", "group")
  layout  <- match.arg(layout)
  measure <- match.arg(measure)
  gn <- names(group)
  if (is.null(groups)) {
    if (length(gn) != 2L)
      .tt_abort("The object has ", length(gn), " groups; name the two to ",
           "subtract via 'groups = c(\"a\", \"b\")'.")
    groups <- gn[1:2]
  }
  if (length(groups) != 2L || !all(groups %in% gn))
    .tt_abort("'groups' must name two groups of the object.")
  a <- group[[groups[1L]]]; b <- group[[groups[2L]]]
  alphabet <- a$alphabet

  if (layout == "tree")
    return(.diff_tree(a, b, groups, alphabet))

  ctx <- intersect(names(a$nodes), names(b$nodes))
  if (!is.null(depth))
    ctx <- ctx[vapply(ctx, .cg_depth, integer(1)) == as.integer(depth)]
  ctx <- ctx[vapply(ctx, function(cc)
    a$nodes[[cc]]$n >= min_count && b$nodes[[cc]]$n >= min_count,
    logical(1))]
  if (length(ctx) == 0L)
    .tt_abort("No shared contexts to compare (try a lower 'min_count' or ",
         "another 'depth').")

  ## significance stars from the permutation comparison, if supplied
  sig_label <- function(p) p
  if (!is.null(comparison)) {
    .tt_check_class(comparison, "transitiontrees_group_comparison", "comparison")
    cp <- comparison$pathways
    sig <- cp$pathway[cp$jsd_padj < alpha]
    sig_label <- function(p) ifelse(p %in% sig, paste0(p, " *"), p)
  }

  pretty <- function(cc) {
    lab <- if (identical(cc, .ROOT)) .ROOT_LABEL else cc
    sig_label(lab)
  }
  ## per-cell value: Pearson residual of group 1 against the no-group
  ## null (observed vs expected from the context's margins), or the raw
  ## probability difference.
  cell_vals <- function(cc) {
    if (measure == "probability")
      return(a$nodes[[cc]]$prob - b$nodes[[cc]]$prob)
    o <- rbind(a$nodes[[cc]]$counts, b$nodes[[cc]]$counts)
    e <- .cg_expected(o)
    r <- (o[1L, ] - e[1L, ]) / sqrt(e[1L, ])
    r[!is.finite(r)] <- 0
    r
  }
  long <- do.call(rbind, lapply(ctx, function(cc) {
    v <- cell_vals(cc)
    data.frame(pathway    = pretty(cc),
               next_state = factor(alphabet, levels = alphabet),
               val        = v, abs_tot = sum(abs(v)),
               stringsAsFactors = FALSE)
  }))
  ## most-different contexts at the top
  row_order <- unique(long$pathway[order(long$abs_tot)])
  long$pathway <- factor(long$pathway, levels = row_order)

  fmt   <- if (measure == "probability") "%+.2f" else "%+.1f"
  thresh <- if (measure == "probability") 0.01 else 0.5
  fill_name <- if (measure == "probability")
    sprintf("P(%s) - P(%s)", groups[1L], groups[2L]) else
    sprintf("Pearson\nresidual\n(%s)", groups[1L])
  M <- max(abs(long$val)); if (M == 0) M <- 1

  p <- ggplot2::ggplot(long, ggplot2::aes(.data$next_state, .data$pathway,
                                          fill = .data$val)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4)
  if (isTRUE(annotate))
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ifelse(abs(.data$val) >= thresh,
                                  sprintf(fmt, .data$val), "")),
      size = 3)
  sub <- if (measure == "probability")
    sprintf("red = more likely in %s, blue = more likely in %s",
            groups[1L], groups[2L]) else
    sprintf(paste0("observed vs expected (no-group null): red = %s does ",
                   "the move MORE than expected, blue = less; |r| > 2 notable"),
            groups[1L])
  p +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-M, M), name = fill_name) +
    ggplot2::labs(
      x = "next state", y = NULL,
      title = sprintf("Difference map: %s vs %s", groups[1L], groups[2L]),
      subtitle = sprintf("%s%s", sub,
        if (is.null(comparison)) "" else
          sprintf("  |  * = FDR < %.2f", alpha))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      axis.text.y   = ggplot2::element_text(family = "mono"),
      panel.grid    = ggplot2::element_blank())
}

#' Difference rendered on the horizontal context-tree layout: a pooled
#' backbone whose nodes and edges are colored by which group reaches a
#' context more (prevalence early - late).
#' @noRd
.diff_tree <- function(a, b, groups, alphabet) {
  trajs_a <- .ct_traj(a$data); trajs_b <- .ct_traj(b$data)
  pooled <- context_tree(c(trajs_a, trajs_b),
                         max_depth = a$max_depth, min_count = a$nmin)
  if (length(pooled$nodes) == 0L)
    .tt_abort("Pooled tree is empty; nothing to plot.", class = "transitiontrees_empty_result")

  contexts <- names(pooled$nodes)
  depths   <- vapply(contexts, .cg_depth, integer(1))
  cc <- .cg_context_counts(stats::setNames(list(trajs_a, trajs_b), groups),
                           groups, contexts, depths, alphabet)
  Ng <- vapply(groups, function(g) sum(cc[[.ROOT]][g, ]), numeric(1))
  delta <- vapply(contexts, function(ctx)
    sum(cc[[ctx]][groups[1L], ]) / Ng[1L] -
    sum(cc[[ctx]][groups[2L], ]) / Ng[2L], numeric(1))
  names(delta) <- contexts

  layout  <- .ct_horizontal_layout(pooled)
  e_paths <- .ct_horizontal_edge_paths(layout, pooled$edges)
  layout$delta  <- delta[layout$context]
  e_paths$delta <- delta[pooled$edges$child][e_paths$edge_id]

  is_root <- layout$context == .ROOT
  body_df <- layout[!is_root, , drop = FALSE]
  root_df <- layout[ is_root, , drop = FALSE]
  is_leaf <- !(layout$context %in% unique(pooled$edges$parent))
  leaf_df <- layout[is_leaf & !is_root, , drop = FALSE]
  leaf_df$label <- leaf_df$context
  internal_df <- layout[!is_leaf & !is_root, , drop = FALSE]
  internal_df$abbr <- vapply(internal_df$context, .pt_last_state,
                             character(1))

  M <- max(abs(layout$delta), na.rm = TRUE); if (!is.finite(M) || M == 0) M <- 1
  div_fill <- ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0,
    limits = c(-M, M),
    name = sprintf("prevalence\n%s - %s", groups[1L], groups[2L]))
  div_col <- ggplot2::scale_colour_gradient2(
    low = "#2166AC", mid = "grey75", high = "#B2182B", midpoint = 0,
    limits = c(-M, M), guide = "none")

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = e_paths,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge_id,
                   linewidth = .data$edge_weight, colour = .data$delta),
      lineend = "round", linejoin = "round") +
    ggplot2::geom_point(
      data = body_df,
      ggplot2::aes(x = .data$x, y = .data$y, size = .data$n,
                   fill = .data$delta),
      shape = 21, colour = "grey25", stroke = 0.2, na.rm = TRUE) +
    ggplot2::geom_text(
      data = internal_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$abbr),
      vjust = 1, nudge_y = -0.05, size = 2.7, colour = "grey20",
      fontface = "bold", na.rm = TRUE) +
    ggplot2::geom_text(
      data = leaf_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = 0, nudge_x = 0.10, size = 3, colour = "grey20",
      family = "mono") +
    ggplot2::geom_point(
      data = root_df, ggplot2::aes(x = .data$x, y = .data$y),
      shape = 16, size = 12, colour = "black") +
    ggplot2::geom_text(
      data = root_df, ggplot2::aes(x = .data$x, y = .data$y),
      label = .ROOT_LABEL, hjust = 1, nudge_x = -0.14, size = 3,
      colour = "grey25", fontface = "bold") +
    div_fill + div_col +
    ggplot2::scale_size_continuous(range = c(4, 14), name = "context count") +
    ggplot2::scale_linewidth_continuous(range = c(0.3, 2.5),
                                        name = "edge flow") +
    ggplot2::scale_x_continuous(
      breaks = seq.int(0L, pooled$max_depth, by = 1L),
      expand = ggplot2::expansion(mult = c(0.17, 0.32))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      legend.position = "right",
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank()) +
    ggplot2::labs(
      x = "depth",
      title = sprintf("Difference tree: %s vs %s", groups[1L], groups[2L]),
      subtitle = sprintf(
        "node & edge colour = which phase reaches the context more (red = %s, blue = %s); size = pooled count",
        groups[1L], groups[2L]))
}
