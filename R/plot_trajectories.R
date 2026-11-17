# ---- Forward trajectory tree ----

#' Forward trajectory-tree node table (internal, testable)
#'
#' Builds the forward prefix tree: every prefix occurring at least
#' \code{min_count} times becomes a node, connected to its one-shorter
#' prefix (the root for first moves). \code{count} is the prefix
#' frequency; \code{eprob} is the model's \eqn{P(\text{last move} \mid
#' \text{history})} read via \code{query_pathway()} (the root distribution
#' for first moves). Prefixes whose parent did not survive the count
#' filter are dropped so the returned tree is always connected to the root.
#'
#' @return A \code{data.frame} with columns \code{node}, \code{parent},
#'   \code{depth}, \code{count}, \code{last}, \code{eprob}.
#' @noRd
.trajectory_data <- function(tree, min_count = 4L) {
  .tt_check_class(tree, "transitiontrees", "tree")
  states <- tree$alphabet
  seqs   <- .ct_traj(tree$data)
  if (length(seqs) == 0L)
    .tt_abort("Tree carries no sequences to draw.")

  prefixes <- unlist(lapply(seqs, function(s)
    vapply(seq_along(s), function(k) paste(s[seq_len(k)], collapse = " -> "),
           character(1))), use.names = FALSE)
  pc   <- table(prefixes)
  keep <- names(pc)[pc >= as.integer(min_count)]
  if (length(keep) == 0L)
    .tt_abort(sprintf("No prefix occurs at least %d times; lower min_count.",
                 as.integer(min_count)),
              class = "transitiontrees_empty_result")

  parent_of <- function(p) {
    m <- strsplit(p, " -> ", fixed = TRUE)[[1L]]
    if (length(m) == 1L) .ROOT_LABEL else paste(m[-length(m)], collapse = " -> ")
  }
  last_of <- function(p)
    utils::tail(strsplit(p, " -> ", fixed = TRUE)[[1L]], 1L)

  node  <- c(.ROOT_LABEL, keep)
  depth <- c(0L, lengths(strsplit(keep, " -> ", fixed = TRUE)))
  count <- c(length(seqs), as.integer(pc[keep]))
  par   <- c(NA_character_, vapply(keep, parent_of, character(1)))
  last  <- c(NA_character_, vapply(keep, last_of,   character(1)))
  names(depth) <- names(count) <- names(par) <- names(last) <- node
  ## drop any prefix whose parent was filtered out (keep the tree connected)
  ok   <- node == .ROOT_LABEL | par %in% node
  node <- node[ok]; depth <- depth[node]; count <- count[node]
  par  <- par[node]; last <- last[node]

  root_prob <- stats::setNames(as.numeric(tree$nodes[[.ROOT]]$prob), states)
  pred_of <- function(n) {
    if (n == .ROOT_LABEL) return(NA_real_)
    p <- par[[n]]
    if (p == .ROOT_LABEL) return(unname(root_prob[[last[[n]]]]))
    hist <- utils::tail(strsplit(p, " -> ", fixed = TRUE)[[1L]], tree$max_depth)
    ## query_pathway() falls back to the longest fitted suffix and always
    ## reaches the root, so it returns a full distribution for any history
    ## built from the tree's own sequences. A state outside the alphabet
    ## has no probability and yields NA (named lookup of a missing name).
    pr <- query_pathway(tree, paste(hist, collapse = " -> "))
    unname(as.numeric(pr[last[[n]]]))
  }
  eprob <- vapply(node, pred_of, numeric(1))

  data.frame(node = node, parent = unname(par[node]), depth = depth[node],
             count = count[node], last = last[node], eprob = eprob[node],
             stringsAsFactors = FALSE, row.names = NULL)
}

#' Draw a forward trajectory tree
#'
#' @description
#' \code{plot_trajectories()} draws a forward trajectory tree of the sequences a
#' tree was fitted on: it starts from a common root and follows each sequence
#' state by state from its beginning, so every node is a prefix (the first
#' state, the first two states, and so on) and its position on the x-axis is the
#' move number. A context tree reads sequences backwards, each node a history
#' that ends at the moment of prediction; the trajectory tree reads them
#' forwards. It shows how sequences fan out from their first state, which
#' opening paths are common, and where they diverge.
#'
#' The tree is colored in one of two ways. With
#' \code{measure = "frequency"} the color of a node is the number of
#' sequences that share the prefix. With \code{measure = "predictability"}
#' it is the probability the fitted model gives to the node's last state
#' after the preceding prefix, \eqn{P(\mathrm{state} \mid
#' \mathrm{history})}{P(state | history)}, with the history truncated to the tree's maximum
#' depth and read through \code{\link{query_pathway}()}; the root
#' distribution is used for first states. Edge width is the number of
#' sequences in both cases, and darker color is a higher value. A pruned
#' tree shows predictability under the pruned model.
#'
#' Only prefixes shared by at least \code{min_count} sequences are drawn.
#' For the backward context tree use \code{\link{plot.transitiontrees}};
#' for the memory needed by one pathway use \code{\link{plot_pruning}()}.
#'
#' @param tree A \code{transitiontrees} object that carries its sequences.
#' @param measure \code{"frequency"} (default) or \code{"predictability"}.
#' @param min_count Smallest number of sequences sharing a prefix for it
#'   to be drawn. Default 4.
#'
#' @return A ggplot object: the prefix tree, move number on the x-axis.
#'
#' @section Conditions:
#' \describe{
#'   \item{\code{transitiontrees_bad_input}}{(error) \code{tree} is not a
#'     \code{transitiontrees} object, or carries no sequences.}
#'   \item{\code{transitiontrees_empty_result}}{(error) no prefix is shared
#'     by \code{min_count} sequences.}
#' }
#'
#' @seealso \code{\link{plot.transitiontrees}} for the backward context
#'   tree; \code{\link{plot_pruning}()}; \code{\link{query_pathway}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' plot_trajectories(tree, min_count = 10)
#' plot_trajectories(prune_tree(tree), measure = "predictability",
#'                   min_count = 10)
#' @export
plot_trajectories <- function(tree,
                              measure = c("frequency", "predictability"),
                              min_count = 4L) {
  .tt_check_class(tree, "transitiontrees", "tree")
  measure <- match.arg(measure)

  states <- tree$alphabet
  D     <- .trajectory_data(tree, min_count = min_count)
  node  <- D$node
  par   <- stats::setNames(D$parent, D$node)
  last  <- stats::setNames(D$last,   D$node)
  depth <- stats::setNames(D$depth,  D$node)
  count <- stats::setNames(D$count,  D$node)
  eprob <- stats::setNames(D$eprob,  D$node)

  ## --- vertical layout: leaves stack, parents centre on children ------
  kids <- split(node[node != .ROOT_LABEL], par[node != .ROOT_LABEL])
  yenv <- new.env(); leaf_y <- 0
  assign_y <- function(nd) {
    k <- kids[[nd]]
    y <- if (is.null(k)) { leaf_y <<- leaf_y + 1; leaf_y }
         else mean(vapply(k, assign_y, numeric(1)))
    assign(nd, y, envir = yenv); y
  }
  assign_y(.ROOT_LABEL)
  ypos <- vapply(node, function(n) get(n, envir = yenv), numeric(1))

  L <- data.frame(node = node, x = depth[node], y = ypos[node],
                  depth = depth[node], count = count[node],
                  last = last[node], eprob = eprob[node],
                  stringsAsFactors = FALSE)
  E <- L[L$node != .ROOT_LABEL, , drop = FALSE]
  E$px <- L$x[match(par[E$node], L$node)]
  E$py <- L$y[match(par[E$node], L$node)]

  ## cosine-smoothed curved edges
  tt <- seq(0, 1, length.out = 40); sm <- (1 - cos(pi * tt)) / 2
  elb <- do.call(rbind, lapply(seq_len(nrow(E)), function(i)
    data.frame(grp = i, x = E$px[i] + tt * (E$x[i] - E$px[i]),
               y = E$py[i] + sm * (E$y[i] - E$py[i]))))

  NB <- L[L$node != .ROOT_LABEL, , drop = FALSE]
  RT <- L[L$node == .ROOT_LABEL, , drop = FALSE]

  ## value driving colour: count or predictability
  nodeval <- if (measure == "frequency") NB$count else NB$eprob
  lim     <- if (measure == "frequency") NULL else c(0, 1)
  num     <- if (measure == "frequency") sprintf("n=%d", NB$count) else
    sprintf("%.0f%%", 100 * NB$eprob)
  leg     <- if (measure == "frequency") "sequences" else "P(move | history)"
  title   <- if (measure == "frequency")
    "Trajectory tree - by frequency" else
    "Trajectory tree - by predictability"
  sub     <- if (measure == "frequency")
    "move inside, count below; colour & width = sequences (higher = darker)" else
    "colour = P(move | history) from the model; width = sequences"

  ## fixed-size rectangular node glyphs (corners computed by hand below)
  HWX <- 0.17; HWY <- 0.26
  ## Pill-shaped (rounded-rectangle) node glyphs, traced in base R:
  ## four quarter-circle corners; geom_polygon joins them with the
  ## straight edges. Corner radius = the smaller half-extent, so the
  ## short sides round into full caps (a capsule).
  RAD <- min(HWX, HWY)
  .corner <- function(cx, cy, ox, oy, a0) {
    a <- seq(a0, a0 + pi / 2, length.out = 6L)
    cbind(cx + ox + RAD * cos(a), cy + oy + RAD * sin(a))
  }
  sqf <- function(df) do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    cx <- df$x[i]; cy <- df$y[i]
    ox <- HWX - RAD; oy <- HWY - RAD               # corner-arc centres
    pts <- rbind(.corner(cx, cy,  ox, -oy, -pi / 2),   # bottom-right
                 .corner(cx, cy,  ox,  oy,  0),        # top-right
                 .corner(cx, cy, -ox,  oy,  pi / 2),   # top-left
                 .corner(cx, cy, -ox, -oy,  pi))       # bottom-left
    data.frame(id = df$node[i], x = pts[, 1L], y = pts[, 2L],
               val = df$valcol[i], stringsAsFactors = FALSE)
  }))

  d <- NB; d$valcol <- nodeval; poly <- sqf(d)
  rt <- RT; rt$valcol <- NA_real_
  vn  <- (nodeval - min(nodeval)) / (max(nodeval) - min(nodeval) + 1e-12)
  txt <- ifelse(vn > 0.45, "white", "grey15")
  ecol <- nodeval[match(E$node, NB$node)][elb$grp]
  ew   <- E$count[elb$grp]

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = elb, ggplot2::aes(x = .data$x, y = .data$y, group = .data$grp,
                               colour = ecol, linewidth = ew),
      lineend = "round") +
    ggplot2::scale_linewidth_continuous(range = c(0.4, 3.6),
                                        name = "flow (sequences)") +
    ggplot2::geom_polygon(
      data = poly, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                                fill = .data$val),
      colour = "grey30", linewidth = 0.3) +
    ggplot2::geom_polygon(
      data = sqf(transform(rt, valcol = NA_real_)),
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
      fill = "grey15") +
    ggplot2::geom_text(data = NB, ggplot2::aes(x = .data$x, y = .data$y,
                                               label = .data$last),
                       colour = txt, size = 2.5, fontface = "bold") +
    ggplot2::geom_text(data = data.frame(NB, num = num),
                       ggplot2::aes(x = .data$x, y = .data$y, label = .data$num),
                       vjust = 1, nudge_y = -0.42, size = 2.5,
                       colour = "grey30") +
    ggplot2::scale_fill_viridis_c(option = "C", direction = -1, limits = lim,
                                  name = leg) +
    ggplot2::scale_colour_viridis_c(option = "C", direction = -1, limits = lim,
                                    guide = "none") +
    ggplot2::scale_x_continuous(breaks = 0:max(L$depth)) +
    ggplot2::labs(x = "move number (time ->)", y = NULL,
                  title = title, subtitle = sub) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9))
}
