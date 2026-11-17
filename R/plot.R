# ---- plot.transitiontrees() dispatcher + dendrogram + horizontal phylogram ----

#' @noRd
.pt_state_palette <- function(state_levels) {
  ## Base-R qualitative palettes from grDevices::hcl.colors().
  ## "Set 2" is the colour-blind-safe 8-class default; falls through
  ## to higher-cardinality palettes for richer alphabets.
  n <- length(state_levels)
  pal_name <- if (n <= 8L) "Set 2"
              else if (n <= 12L) "Set 3"
              else "Dark 3"
  pal <- grDevices::hcl.colors(max(2L, n), palette = pal_name)
  setNames(pal[seq_len(n)], state_levels)
}

#' @noRd
.ct_radial_layout <- function(tree) {
  ## Radial coordinates for each node:
  ##   r     = depth
  ##   theta = leaf-rank for leaves (evenly spaced 0..2pi),
  ##           mean of children's theta for internal nodes.
  ## Result includes (x, y) = (r * cos theta, r * sin theta).
  if (length(tree$nodes) == 0L) return(NULL)

  children_by_parent <- .pt_children_of(tree)
  children_of <- function(parent)
    children_by_parent[[parent]] %||% character(0)

  ## DFS to collect leaves left-to-right.
  leaves <- character(0)
  walk <- function(ctx) {
    ch <- children_of(ctx)
    if (length(ch) == 0L) {
      leaves[[length(leaves) + 1L]] <<- ctx
    } else {
      lapply(ch, walk)
    }
    invisible(NULL)
  }
  walk(.ROOT)

  theta_env <- new.env(hash = TRUE, parent = emptyenv())

  if (length(leaves) <= 1L) {
    ## Degenerate trees (root only, or root + single leaf) get a fixed angle.
    Map(function(nm) assign(nm, 0, envir = theta_env), names(tree$nodes))
  } else {
    leaf_theta <- 2 * pi * (seq_along(leaves) - 0.5) / length(leaves)
    Map(function(nm, t) assign(nm, t, envir = theta_env),
        leaves, leaf_theta)

    assign_theta <- function(ctx) {
      ch <- children_of(ctx)
      if (length(ch) == 0L) return(theta_env[[ctx]])
      child_thetas <- vapply(ch, assign_theta, numeric(1))
      t <- mean(child_thetas)
      assign(ctx, t, envir = theta_env)
      t
    }
    assign_theta(.ROOT)
  }

  do.call(rbind, lapply(names(tree$nodes), function(ctx) {
    d  <- tree$nodes[[ctx]]$depth
    th <- theta_env[[ctx]]
    data.frame(context = ctx, depth = d, r = d, theta = th,
               x = d * cos(th), y = d * sin(th),
               n = tree$nodes[[ctx]]$n,
               stringsAsFactors = FALSE)
  }))
}

#' @noRd
.ct_radial_edge_paths <- function(layout, edges, n_arc = 24L) {
  ## Phylogram-style elbow edges with an edge-weight column equal to
  ## the child's count, so callers can map it to linewidth.
  if (nrow(edges) == 0L)
    return(data.frame(edge_id = integer(0), order = integer(0),
                      x = numeric(0), y = numeric(0),
                      edge_weight = numeric(0)))

  by_ctx <- setNames(seq_len(nrow(layout)), layout$context)
  paths  <- Map(function(p_ctx, c_ctx, eid) {
    rp <- layout$r[[by_ctx[[p_ctx]]]];  tp <- layout$theta[[by_ctx[[p_ctx]]]]
    rc <- layout$r[[by_ctx[[c_ctx]]]];  tc <- layout$theta[[by_ctx[[c_ctx]]]]
    nc <- layout$n[[by_ctx[[c_ctx]]]]
    arc_t <- seq.int(tp, tc, length.out = n_arc)
    data.frame(
      edge_id     = eid,
      order       = seq_len(n_arc + 1L),
      x           = c(rp * cos(arc_t), rc * cos(tc)),
      y           = c(rp * sin(arc_t), rc * sin(tc)),
      edge_weight = nc
    )
  }, edges$parent, edges$child, seq_len(nrow(edges)))
  do.call(rbind, paths)
}

#' @noRd
.ct_horizontal_layout <- function(tree) {
  ## Cartesian coordinates for a left-to-right phylogram:
  ##   x = depth (root at x=0, leaves at x=max_depth)
  ##   y = leaf-rank for leaves (evenly spaced in [0, 1]),
  ##       mean of children's y for internal nodes.
  if (length(tree$nodes) == 0L) return(NULL)

  children_by_parent <- .pt_children_of(tree)
  children_of <- function(parent)
    children_by_parent[[parent]] %||% character(0)

  ## DFS to collect leaves top-to-bottom.
  leaves <- character(0)
  walk <- function(ctx) {
    ch <- children_of(ctx)
    if (length(ch) == 0L) {
      leaves[[length(leaves) + 1L]] <<- ctx
    } else {
      lapply(ch, walk)
    }
    invisible(NULL)
  }
  walk(.ROOT)

  y_env <- new.env(hash = TRUE, parent = emptyenv())

  if (length(leaves) <= 1L) {
    Map(function(nm) assign(nm, 0.5, envir = y_env), names(tree$nodes))
  } else {
    leaf_y <- seq.int(0, 1, length.out = length(leaves))
    ## Reverse so the first leaf is at the top of the plot canvas.
    leaf_y <- rev(leaf_y)
    Map(function(nm, t) assign(nm, t, envir = y_env), leaves, leaf_y)

    assign_y <- function(ctx) {
      ch <- children_of(ctx)
      if (length(ch) == 0L) return(y_env[[ctx]])
      child_ys <- vapply(ch, assign_y, numeric(1))
      t <- mean(child_ys)
      assign(ctx, t, envir = y_env)
      t
    }
    assign_y(.ROOT)
  }

  do.call(rbind, lapply(names(tree$nodes), function(ctx) {
    d <- tree$nodes[[ctx]]$depth
    y <- y_env[[ctx]]
    data.frame(context = ctx, depth = d, x = d, y = y,
               n = tree$nodes[[ctx]]$n,
               stringsAsFactors = FALSE)
  }))
}

#' @noRd
.ct_horizontal_edge_paths <- function(layout, edges, n_pt = 40L) {
  ## Smooth S-curve edges: a cosine smoothstep carries each branch from
  ## (x_parent, y_parent) to (x_child, y_child), leaving the parent and
  ## arriving at the child horizontally, so branches read as flowing
  ## curves rather than mechanical right-angle elbows. Each edge is an
  ## n_pt-vertex polyline; edge_weight = child's count.
  if (nrow(edges) == 0L)
    return(data.frame(edge_id = integer(0), order = integer(0),
                      x = numeric(0), y = numeric(0),
                      edge_weight = numeric(0)))

  by_ctx <- setNames(seq_len(nrow(layout)), layout$context)
  t  <- seq(0, 1, length.out = n_pt)
  sm <- (1 - cos(pi * t)) / 2          # 0 -> 1 smoothstep
  paths <- Map(function(p_ctx, c_ctx, eid) {
    xp <- layout$x[[by_ctx[[p_ctx]]]]; yp <- layout$y[[by_ctx[[p_ctx]]]]
    xc <- layout$x[[by_ctx[[c_ctx]]]]; yc <- layout$y[[by_ctx[[c_ctx]]]]
    nc <- layout$n[[by_ctx[[c_ctx]]]]
    data.frame(
      edge_id     = eid,
      order       = seq_len(n_pt),
      x           = xp + t  * (xc - xp),
      y           = yp + sm * (yc - yp),
      edge_weight = nc
    )
  }, edges$parent, edges$child, seq_len(nrow(edges)))
  do.call(rbind, paths)
}

#' @noRd
.plot_dendrogram <- function(x,
                              point_size_range = c(5, 16),
                              edge_size_range  = c(0.3, 2.5),
                              ...) {
  if (length(x$nodes) == 0L)
    .tt_abort("Cannot plot an empty tree.", class = "transitiontrees_empty_result")

  layout  <- .ct_radial_layout(x)
  e_paths <- .ct_radial_edge_paths(layout, x$edges)

  state_levels   <- x$alphabet
  last_state_str <- vapply(layout$context, .pt_last_state, character(1))
  layout$state   <- factor(last_state_str, levels = state_levels)
  layout$abbr    <- as.character(layout$state)
  pal <- .pt_state_palette(state_levels)

  ## Root is rendered as a separate layer (last, on top, fixed size,
  ## fixed black fill) so it never gets buried under depth-1 children
  ## drawn on top of the (0, 0) origin.
  is_root  <- layout$context == .ROOT
  body_df  <- layout[!is_root, , drop = FALSE]
  root_df  <- layout[ is_root, , drop = FALSE]

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = e_paths,
      ggplot2::aes(x = .data$x, y = .data$y,
                   group = .data$edge_id,
                   linewidth = .data$edge_weight),
      colour = "grey60", lineend = "round"
    ) +
    ggplot2::geom_point(
      data = body_df,
      ggplot2::aes(x = .data$x, y = .data$y,
                   size = .data$n, fill = .data$state),
      shape = 21, colour = "grey25", stroke = 0.15, na.rm = TRUE
    ) +
    ggplot2::geom_text(
      data = body_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$abbr),
      size = 2.4, colour = "grey15", fontface = "bold", na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = root_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 16, size = max(point_size_range) * 0.8,
      colour = "black"
    ) +
    ggplot2::scale_size_continuous(range = point_size_range,
                                   name = "context count") +
    ggplot2::scale_linewidth_continuous(range = edge_size_range,
                                         name = "edge flow") +
    ggplot2::scale_fill_manual(values = pal, name = "State",
                               na.translate = FALSE,
                               drop = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title    = "Context tree (radial)",
      subtitle = sprintf(
        "%d nodes, depth <= %d;  node size = count;  edge thickness = flow",
        length(x$nodes), x$max_depth)
    )
}

#' @noRd
.plot_horizontal <- function(x,
                              point_size_range = c(4, 14),
                              edge_size_range  = c(0.3, 2.5),
                              show_prediction  = TRUE,
                              ...) {
  if (length(x$nodes) == 0L)
    .tt_abort("Cannot plot an empty tree.", class = "transitiontrees_empty_result")

  layout  <- .ct_horizontal_layout(x)
  e_paths <- .ct_horizontal_edge_paths(layout, x$edges)

  ## Node colour = the MOST-RECENT move in the context (the rightmost
  ## token), i.e. "what just happened". This is the intuitive reading and
  ## makes each branch off a depth-1 hub share one colour (every "X ->
  ## Specify" ends in Specify, so the whole Specify branch is one colour).
  state_levels     <- x$alphabet
  recent_state_str <- vapply(layout$context, function(ctx) {
    if (identical(ctx, .ROOT)) return(NA_character_)
    parts <- strsplit(ctx, " -> ", fixed = TRUE)[[1L]]
    parts[[length(parts)]]
  }, character(1))
  layout$state <- factor(recent_state_str, levels = state_levels)
  pal <- .pt_state_palette(state_levels)

  is_root  <- layout$context == .ROOT
  body_df  <- layout[!is_root, , drop = FALSE]
  root_df  <- layout[ is_root, , drop = FALSE]
  ## Labels: full pathway in arrow notation drawn UNDER every non-root
  ## node. Simple by default (context only); set show_prediction = TRUE to
  ## add the modal prediction "(state pct%)" on a second line.
  body_df$label <- if (isTRUE(show_prediction)) {
    vapply(body_df$context, function(ctx) {
      pr   <- x$nodes[[ctx]]$prob
      best <- which.max(pr)
      sprintf("%s\n(%s %.0f%%)", ctx, x$alphabet[best], 100 * pr[best])
    }, character(1))
  } else {
    body_df$context
  }

  ## Per-node label clearance. The size aesthetic maps count -> area
  ## (point diameter ~ sqrt of the rescaled count), so the large depth-1
  ## hubs need their label pushed further below the glyph than the small
  ## leaf nodes. Offset each label by an amount that grows with the
  ## node's plotted radius instead of a single fixed nudge.
  n_rng <- if (nrow(body_df) > 0L) range(body_df$n) else c(0, 0)
  size_aes <- if (nrow(body_df) > 0L && diff(n_rng) > 0)
    point_size_range[1L] +
      sqrt((body_df$n - n_rng[1L]) / diff(n_rng)) * diff(point_size_range)
  else rep(mean(point_size_range), nrow(body_df))
  body_df$label_y <- body_df$y - (0.012 + size_aes * 0.0035)

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = e_paths,
      ggplot2::aes(x = .data$x, y = .data$y,
                   group = .data$edge_id,
                   linewidth = .data$edge_weight),
      colour = "grey60", lineend = "round", linejoin = "round"
    ) +
    ggplot2::geom_point(
      data = body_df,
      ggplot2::aes(x = .data$x, y = .data$y,
                   size = .data$n, fill = .data$state),
      shape = 21, colour = "grey25", stroke = 0.2, na.rm = TRUE
    ) +
    ggplot2::geom_text(
      data = body_df,
      ggplot2::aes(x = .data$x, y = .data$label_y, label = .data$label),
      hjust = 0.5, vjust = 1, size = 2.5,
      lineheight = 0.9, colour = "grey20", family = "mono", na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = root_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 16, size = max(point_size_range) * 0.85,
      colour = "black"
    ) +
    ggplot2::geom_text(
      data = root_df,
      ggplot2::aes(x = .data$x, y = .data$y), label = .ROOT_LABEL,
      hjust = 1, nudge_x = -0.14, size = 3,
      colour = "grey25", fontface = "bold"
    ) +
    ggplot2::scale_size_continuous(range = point_size_range,
                                   name = "context count") +
    ggplot2::scale_linewidth_continuous(range = edge_size_range,
                                         name = "edge flow") +
    ggplot2::scale_fill_manual(values = pal, name = "Most recent move",
                               na.translate = FALSE,
                               drop = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = seq.int(0L, x$max_depth, by = 1L),
      expand = ggplot2::expansion(mult = c(0.17, 0.32))
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(colour = "grey40", size = 9),
      legend.position  = "right",
      axis.title.y     = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_blank(),
      axis.ticks.y     = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x        = "depth",
      title    = "Context tree",
      subtitle = "colour = most recent move; size = count; ( ) = predicted next move"
    )
}

#' Plot a context tree
#'
#' @description
#' \code{plot()} draws a fitted context tree. A context tree organizes the
#' histories that precede each state: the root holds the overall distribution of
#' states; each child extends its parent's history by one earlier state; and
#' every node carries the distribution of the next state given that history. A
#' drawing of the tree shows which histories occur often enough to be modeled,
#' how far back the model's memory reaches along each branch, and where branches
#' stop because a longer history was too rare or, after pruning, added nothing
#' to the prediction.
#'
#' Four drawings are available. \code{"horizontal"} (the default) is a
#' left-to-right tree with depth on the x-axis; each node is labeled with its
#' full history and, by default, its most likely next state and that state's
#' probability, and is filled by the most recent state of the history.
#' \code{"dendrogram"} is a radial tree with the root at the center, each node
#' filled by the state its depth adds to the history. \code{"icicle"} is a
#' sunburst: rings are depths, the angular width of an arc is proportional to
#' the summed count of the deepest histories beneath it, and the fill is the
#' state the ring adds to the history. \code{"interactive"} is a
#' \pkg{visNetwork} widget with draggable nodes and a tooltip that gives each
#' node's count and full next-state distribution. In the node-link drawings,
#' node size is the count of the history and edge width the count of the child
#' history. The first three drawings are ggplot objects; the interactive one is
#' an htmlwidget.
#'
#' The context tree reads backwards in time: every history ends at the
#' moment of prediction. For the forward view, sequences followed from
#' their first state, use \code{\link{plot_trajectories}()}. The
#' next-state probabilities of many pathways are compared more easily in
#' \code{\link{plot_pathways}()} or \code{\link{plot_distributions}()},
#' and differences between two groups in \code{\link{plot_difference}()}.
#'
#' @param x A \code{transitiontrees} object.
#' @param style \code{"horizontal"} (default), \code{"dendrogram"},
#'   \code{"icicle"} or \code{"interactive"}.
#' @param point_size_range Smallest and largest node size. Default
#'   \code{NULL}: \code{c(4, 14)} for \code{"horizontal"}, \code{c(5, 16)}
#'   for \code{"dendrogram"}, \code{c(10, 45)} pixels for
#'   \code{"interactive"}. Not used by \code{"icicle"}.
#' @param edge_size_range Thinnest and thickest edge. Default \code{NULL}:
#'   \code{c(0.3, 2.5)} for the static drawings, \code{c(1, 10)} pixels for
#'   \code{"interactive"}. Not used by \code{"icicle"}.
#' @param ... Options of the chosen drawing. \code{"horizontal"}:
#'   \code{show_prediction} (default \code{TRUE}) adds the most likely next
#'   state to each label. \code{"icicle"}: \code{label_min_fraction}
#'   (default 0.03) and \code{label_abbrev_fraction} (default 0.08), the
#'   arc shares below which a label is omitted or abbreviated;
#'   \code{abbrev_len} (default 2) or \code{abbrev_fun} (default
#'   \code{NULL}) for the abbreviation; \code{gap_size} (default 5) and
#'   \code{gap_size_within} (default 0.3), the widths of the white lines
#'   between first-level wedges and between arcs; \code{jitter} (default
#'   0.04) and \code{jitter_seed} (default 1), a small random variation of
#'   the outer edge (0 for none). \code{"interactive"}: \code{width} and
#'   \code{height} of the widget.
#'
#' @return A ggplot object for \code{"horizontal"}, \code{"dendrogram"}
#'   and \code{"icicle"}; a \pkg{visNetwork} htmlwidget for
#'   \code{"interactive"}.
#'
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_empty_result}}{(error) the tree has no
#'     nodes.}
#'   \item{\code{transitiontrees_missing_package}}{(error)
#'     \code{"interactive"} without the \pkg{visNetwork} package
#'     installed.}
#' }
#'
#' @seealso \code{\link{plot.transitiontrees_group}} for one plot per
#'   group; \code{\link{plot_trajectories}()} for the forward tree;
#'   \code{\link{plot_pathways}()}, \code{\link{plot_distributions}()} and
#'   \code{\link{plot_pruning}()} for next-state distributions.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' plot(tree)
#' plot(tree, show_prediction = FALSE)
#' plot(tree, style = "dendrogram")
#' plot(tree, style = "icicle")
#' if (requireNamespace("visNetwork", quietly = TRUE)) {
#'   widget <- plot(tree, style = "interactive")
#' }
#' @export
plot.transitiontrees <- function(x,
                          style = c("horizontal", "dendrogram",
                                     "icicle", "interactive"),
                          point_size_range = NULL,
                          edge_size_range  = NULL,
                          ...) {
  style <- match.arg(style)
  if (style == "interactive" &&
      !requireNamespace("visNetwork", quietly = TRUE)) {
    .tt_abort("plot(style = 'interactive') requires 'visNetwork'. ",
         "Install it, or use style = 'dendrogram' / 'horizontal' / ",
         "'icicle'.",
              class = "transitiontrees_missing_package")
  }
  switch(style,
    icicle      = .plot_icicle(x, ...),
    interactive = .plot_interactive(x,
                   point_size_range = point_size_range %||% c(10, 45),
                   edge_size_range  = edge_size_range  %||% c(1, 10),
                   ...),
    dendrogram  = .plot_dendrogram(x,
                   point_size_range = point_size_range %||% c(5, 16),
                   edge_size_range  = edge_size_range  %||% c(0.3, 2.5),
                   ...),
    horizontal  = .plot_horizontal(x,
                   point_size_range = point_size_range %||% c(4, 14),
                   edge_size_range  = edge_size_range  %||% c(0.3, 2.5),
                   ...))
}

#' Plot the context tree of each group
#'
#' @description
#' \code{plot()} draws every tree of a \code{transitiontrees_group} with
#' \code{\link{plot.transitiontrees}}, one figure per group with the group name
#' as caption, so that the memory structure of the groups can be compared side
#' by side. Each figure is printed as it is made. To test or map the differences
#' between groups use \code{\link{compare_groups}()} and
#' \code{\link{plot_difference}()}.
#'
#' @param x A \code{transitiontrees_group} object.
#' @param ... Passed to \code{\link{plot.transitiontrees}}, such as
#'   \code{style}.
#' @return Invisibly, a named list with one plot per group, in group
#'   order: ggplot objects, or htmlwidgets for
#'   \code{style = "interactive"}.
#' @seealso \code{\link{plot.transitiontrees}};
#'   \code{\link{plot_difference}()} and \code{\link{compare_groups}()}
#'   for the differences between groups.
#' @examples
#' grp <- context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'                     action = "Action", group = "Achiever", max_depth = 1)
#' plot(grp, show_prediction = FALSE)
#' @export
plot.transitiontrees_group <- function(x, ...) {
  plots <- lapply(names(x), function(nm) {
    p <- plot(x[[nm]], ...)
    if (inherits(p, "ggplot"))
      p <- p + ggplot2::labs(caption = paste0("group: ", nm))
    print(p)
    p
  })
  invisible(setNames(plots, names(x)))
}
