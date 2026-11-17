# ---- Pathway-centric visualisations ----

#' Plot pathways as a heatmap of next-state probabilities
#'
#' @description
#' \code{plot_pathways()} draws the pathway table of a context tree as a
#' heatmap: one row per pathway, one column per next state, and the fitted
#' probability \eqn{P(\mathrm{next} \mid \mathrm{pathway})}{P(next | pathway)}
#' as fill. The pathway table lists, for every history, how often it occurs and
#' what tends to follow it. The most likely next state of each row is printed in
#' bold, other probabilities of at least 0.05 in plain text, and pathways whose
#' most likely next state differs from that of their shorter suffix are marked
#' with a leading \code{>}. One display then shows which pathways are common,
#' which lead almost always to one state, which leave the next state open, and
#' where a longer history changes the prediction that a first-order model would
#' make.
#'
#' Rows are the \code{top} pathways by count, by divergence from the
#' shorter suffix, or by depth (the root is not shown). For the full
#' distribution of a few contexts as bar charts use
#' \code{\link{plot_distributions}()}; to rank pathways by how much their
#' history matters use \code{\link{plot_divergence}()}.
#'
#' @param tree A \code{transitiontrees} object.
#' @param top Maximum number of pathways shown. Default 20.
#' @param sort_by Order of the rows: \code{"count"} (default),
#'   \code{"divergence"} or \code{"depth"} (deepest first).
#' @param min_count Pathways with a lower count are dropped. Default 5.
#' @param show_flips Mark pathways whose most likely next state differs
#'   from their parent's. Default \code{TRUE}.
#' @param title Plot title. Default \code{NULL}: a title naming the
#'   ordering.
#' @param ... Ignored.
#' @return A ggplot object: a heatmap of pathways by next states.
#'
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_bad_input}}{(error) \code{tree} is not a
#'     \code{transitiontrees} object.}
#'   \item{\code{transitiontrees_empty_result}}{(error) no pathway reaches
#'     \code{min_count}.}
#' }
#'
#' @seealso \code{\link{tree_pathways}()} for the table;
#'   \code{\link{plot_distributions}()}, \code{\link{plot_divergence}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' plot_pathways(tree, top = 12)
#' plot_pathways(tree, sort_by = "divergence", top = 12)
#' @export
plot_pathways <- function(tree,
                          top = 20L,
                          sort_by = c("count", "divergence", "depth"),
                          min_count = 5L,
                          show_flips = TRUE,
                          title = NULL, ...) {
  .tt_check_class(tree, "transitiontrees", "tree")
  sort_by <- match.arg(sort_by)

  pw <- tree_pathways(tree, min_count = min_count, sort_by = sort_by,
                          decreasing = TRUE)
  pw <- pw[pw$pathway != .ROOT_LABEL, , drop = FALSE]
  pw <- utils::head(pw, n = as.integer(top))
  if (nrow(pw) == 0L)
    .tt_abort("No pathways meet the threshold.", class = "transitiontrees_empty_result")

  ## Build the long form: one row per (pathway, next_state)
  states <- tree$alphabet
  long <- do.call(rbind, lapply(seq_len(nrow(pw)), function(i) {
    ctx <- pw$pathway[i]
    info <- tree$nodes[[ctx]]
    data.frame(
      pathway    = ctx,
      next_state = factor(states, levels = states),
      prob       = info$prob,
      stringsAsFactors = FALSE
    )
  }))

  ## Display label: prefix flippers with a caret
  pw$label <- ifelse(isTRUE(show_flips) & !is.na(pw$changes_prediction) &
                       pw$changes_prediction,
                     sprintf("> %s", pw$pathway),
                     pw$pathway)
  ## Reverse ordering so highest-ranked appears at top of the y axis
  pw <- pw[seq.int(nrow(pw), 1L), , drop = FALSE]
  long$label <- factor(pw$label[match(long$pathway, pw$pathway)],
                        levels = pw$label)

  ## Highlight modal cells (vapply pins the result to logical)
  long$is_modal <- vapply(seq_len(nrow(long)), function(i) {
    info <- tree$nodes[[long$pathway[i]]]
    states[which.max(info$prob)] == as.character(long$next_state[i])
  }, logical(1))

  if (is.null(title)) {
    sort_label <- switch(sort_by,
                         count      = "by frequency",
                         divergence = "by divergence from shorter history",
                         depth      = "by length (deepest first)")
    title <- sprintf("Top %d pathways %s", nrow(pw), sort_label)
  }

  ggplot2::ggplot(long, ggplot2::aes(x = .data$next_state,
                                     y = .data$label,
                                     fill = .data$prob)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      data = long[long$is_modal, , drop = FALSE],
      ggplot2::aes(label = sprintf("%.2f", .data$prob)),
      colour = "white", fontface = "bold", size = 3.2) +
    ggplot2::geom_text(
      data = long[!long$is_modal & long$prob >= 0.05, , drop = FALSE],
      ggplot2::aes(label = sprintf("%.2f", .data$prob)),
      colour = "grey20", size = 3.0) +
    ggplot2::scale_fill_gradient(
      low = "#deebf7", high = "#08519c",
      name = "P(next | pathway)",
      limits = c(0, 1)) +
    ggplot2::labs(
      x = "Next state", y = NULL,
      title = title,
      subtitle = sprintf(
        "%s; n_min = %d; %d flips marked with leading >",
        switch(sort_by,
               count      = "ranked by pathway frequency",
               divergence = "ranked by divergence from (k-1)-suffix",
               depth      = "ranked by depth"),
        min_count, sum(pw$changes_prediction, na.rm = TRUE))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(family = "mono"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}


#' Plot the pathways whose history changes the prediction most
#'
#' @description
#' \code{plot_divergence()} ranks pathways by the Kullback-Leibler divergence,
#' in bits, between their next-state distribution and that of their shorter
#' suffix, the context without the earliest state. In a variable-order model a
#' longer history matters only where it changes what is expected next. Each
#' pathway is a lollipop whose point size is its count; points are colored by
#' whether the most likely next state differs from the suffix's, and those
#' pathways are annotated with the new most likely state and its probability.
#' The chart identifies the histories that carry information a shorter memory
#' would lose.
#'
#' The pathways are those returned by \code{\link{divergent_pathways}()}.
#' A large divergence on a small count can be a sampling accident; check
#' reproducibility with \code{\link{bootstrap_pathways}()}, and the
#' pruning decision for one pathway with \code{\link{plot_pruning}()}.
#'
#' @param tree A \code{transitiontrees} object.
#' @param top Number of pathways shown. Default 15.
#' @param min_count Pathways with a lower count are dropped. Default 5.
#' @param title Plot title. Default \code{NULL}: a title giving the number
#'   of pathways.
#' @param ... Ignored.
#' @return A ggplot object: one lollipop per pathway, divergence (bits) on
#'   the x-axis.
#'
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_bad_input}}{(error) \code{tree} is not a
#'     \code{transitiontrees} object.}
#'   \item{\code{transitiontrees_empty_result}}{(error) no pathway reaches
#'     \code{min_count}.}
#' }
#'
#' @seealso \code{\link{divergent_pathways}()} for the table;
#'   \code{\link{plot_pathways}()}, \code{\link{plot_pruning}()},
#'   \code{\link{bootstrap_pathways}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' plot_divergence(tree, top = 10)
#' @export
plot_divergence <- function(tree, top = 15L, min_count = 5L,
                            title = NULL, ...) {
  .tt_check_class(tree, "transitiontrees", "tree")
  pw <- divergent_pathways(tree, top = top, min_count = min_count)
  if (nrow(pw) == 0L) .tt_abort("No pathways meet the threshold.", class = "transitiontrees_empty_result")
  pw$pathway <- factor(pw$pathway, levels = rev(pw$pathway))

  if (is.null(title))
    title <- sprintf("Top %d pathways by KL divergence", nrow(pw))

  ggplot2::ggplot(pw, ggplot2::aes(x = .data$divergence,
                                   y = .data$pathway)) +
    ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = .data$pathway),
                          colour = "grey70") +
    ggplot2::geom_point(ggplot2::aes(size = .data$count,
                                     colour = .data$changes_prediction)) +
    ggplot2::geom_text(
      data = pw[!is.na(pw$changes_prediction) & pw$changes_prediction, ],
      ggplot2::aes(label = sprintf("%.0f%% %s",
                                   100 * .data$next_probability,
                                   .data$likely_next)),
      hjust = -0.15, size = 3.2, colour = "#D55E00") +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = "#0072B2", `TRUE` = "#D55E00"),
      labels = c("modal next-state same",
                 "modal next-state flips"),
      name = NULL,
      na.value = "grey60") +
    ggplot2::scale_size_continuous(name = "pathway count") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.30))) +
    ggplot2::labs(
      x = "KL divergence (bits)  vs (k-1)-suffix",
      y = NULL,
      title = title) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}
