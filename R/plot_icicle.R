# ---- Static icicle / partition visualisation ----

#' @noRd
.pt_arc_text_angle <- function(x, y) {
  ## Tangent angle (in degrees) for text running along the arc of a
  ## circular sunburst, with an upside-down flip so text stays
  ## right-side-up. (x, y) are the layout-assigned node centroids.
  ang <- atan2(y, x) * 180 / pi - 90
  ifelse(ang < -90, ang + 180,
    ifelse(ang > 90, ang - 180, ang))
}

#' @noRd
.plot_icicle <- function(tree,
                          label_min_fraction    = 0.03,
                          label_abbrev_fraction = 0.08,
                          gap_size              = 5,
                          gap_size_within       = 0.3,
                          abbrev_len            = 2L,
                          abbrev_fun            = NULL,
                          jitter                = 0.04,
                          jitter_seed           = 1L,
                          ...) {
  ## Two-tier (default) vs three-tier label policy:
  ## - If `label_abbrev_fraction` is NULL (default), the icicle is
  ##   binary: any slice >= label_min_fraction gets the full state
  ##   name; anything smaller gets no label.
  ## - If `label_abbrev_fraction` is supplied AND > label_min_fraction,
  ##   slices in [label_min_fraction, label_abbrev_fraction) get an
  ##   abbreviation (via `abbrev_fun` or `substr` to `abbrev_len`);
  ##   slices >= label_abbrev_fraction still get the full name.
  if (is.null(label_abbrev_fraction))
    label_abbrev_fraction <- label_min_fraction
  stopifnot(
    is.numeric(label_min_fraction),    label_min_fraction    >= 0,
    is.numeric(label_abbrev_fraction), label_abbrev_fraction >= 0,
    is.numeric(gap_size),              gap_size        >= 0,
    is.numeric(gap_size_within),       gap_size_within >= 0,
    is.numeric(abbrev_len),            abbrev_len >= 1L,
    is.null(abbrev_fun) || is.function(abbrev_fun),
    is.numeric(jitter),                jitter >= 0
  )
  ## Abbreviation function: by default a simple substr to abbrev_len
  ## characters. Pass abbrev_fun = base::abbreviate (vowel-dropping
  ## smart shortening with uniqueness guarantees) or any custom
  ## function taking a character vector and returning short codes.
  .abbrev <- if (!is.null(abbrev_fun)) abbrev_fun
             else function(s) substr(s, 1L, as.integer(abbrev_len))

  node_names <- names(tree$nodes)
  is_leaf <- !(node_names %in% unique(tree$edges$parent))

  ## Each tile's fill encodes the *state* it represents — the last
  ## token of the pathway. Root has no last state.
  last_state <- vapply(node_names, .pt_last_state, character(1))
  state_levels <- tree$alphabet
  state_factor <- factor(last_state, levels = state_levels)

  pt_depth <- vapply(tree$nodes, function(n) n$depth, integer(1))
  pt_count <- vapply(tree$nodes, function(n) as.numeric(n$n), numeric(1))

  ## Label policy (size-based, not depth-based):
  ##   arc_fraction < label_min_fraction    → no label  (too thin to read)
  ##   arc_fraction < label_abbrev_fraction → abbreviation (abbrev_len chars)
  ##   otherwise                            → full state name
  ##
  ## arc_fraction is the *actual* visual arc fraction in the circular
  ## partition: each leaf contributes its count as its weight; each
  ## internal node's visual arc is the sum of its descendant leaves'
  ## counts. Using the node's own count would overstate the fraction
  ## for internal nodes (parent counts > sum of leaf-descendant
  ## counts whenever some observations end before reaching the leaf).
  total_leaf_count <- sum(pt_count[is_leaf])
  children_by_parent <- if (nrow(tree$edges) > 0L)
                          split(tree$edges$child, tree$edges$parent)
                        else list()
  .subtree_leaf <- function(ctx) {
    ch <- children_by_parent[[ctx]]
    if (is.null(ch) || length(ch) == 0L) return(pt_count[[ctx]])
    sum(vapply(ch, .subtree_leaf, numeric(1)))
  }
  subtree_leaf_w <- vapply(node_names, .subtree_leaf, numeric(1))
  arc_fraction <- subtree_leaf_w / max(total_leaf_count, 1)
  ## Pre-compute abbreviations for every alphabet state once so the
  ## abbrev_fun sees the whole vocabulary at once (matters for
  ## base::abbreviate which needs the full set to enforce
  ## uniqueness). Then index by state.
  state_abbrev <- as.character(.abbrev(state_levels))
  abbrev_lookup <- setNames(state_abbrev, state_levels)
  abbreviated   <- ifelse(is.na(last_state) | last_state == .ROOT_LABEL,
                          NA_character_,
                          abbrev_lookup[last_state])
  label_text <- ifelse(
    pt_depth == 0L | last_state == .ROOT_LABEL | is.na(last_state) |
      arc_fraction < label_min_fraction,
    NA_character_,
    ifelse(arc_fraction < label_abbrev_fraction,
           abbreviated,
           last_state)
  )

  ## Per-node frame used as the plot's default data, so introspection
  ## (e.g. `plot(...)$data`) still exposes `label_text` / `arc_fraction`.
  nodes <- data.frame(
    name         = node_names,
    state        = unname(state_factor),
    count        = unname(pt_count),
    pt_depth     = unname(pt_depth),
    label_text   = unname(label_text),
    arc_fraction = unname(arc_fraction),
    stringsAsFactors = FALSE,
    row.names    = NULL
  )

  pal <- .pt_state_palette(state_levels)

  ## ---- Circular-partition layout, computed in base R ----
  ## Each node owns an angular interval whose width is its subtree's
  ## leaf weight (internal arcs = sum of descendant leaves, matching the
  ## classic "partition" layout). Children are laid out consecutively
  ## inside their parent's interval, in edge order, so the result is
  ## deterministic. Depth maps to radius: ring d spans [d, d + 1].
  total <- max(total_leaf_count, 1)
  node_a0 <- setNames(numeric(length(node_names)), node_names)
  node_a1 <- node_a0
  assign_angles <- function(ctx, a0) {
    node_a0[[ctx]] <<- a0
    node_a1[[ctx]] <<- a0 + subtree_leaf_w[[ctx]]
    ch <- children_by_parent[[ctx]]
    if (!is.null(ch) && length(ch) > 0L) {
      w_ch <- subtree_leaf_w[ch]
      offs <- a0 + c(0, cumsum(w_ch)[-length(w_ch)])
      invisible(Map(assign_angles, ch, offs))
    }
  }
  assign_angles(.ROOT, 0)

  theta0 <- node_a0 / total * 2 * pi
  theta1 <- node_a1 / total * 2 * pi
  ## Each depth ring is `ring_thickness` wide, with a `ring_gap` of empty
  ## space before the next ring (ring d spans [d * (thickness + gap),
  ## + thickness]). Slot spacing = thickness + gap.
  ring_thickness <- 0.68
  ring_gap       <- 0.10
  rr0    <- as.numeric(pt_depth) * (ring_thickness + ring_gap)
  rr1    <- rr0 + ring_thickness

  ## Optional radial jitter on leaf tiles: each leaf's outer radius is
  ## perturbed by up to `jitter` of a ring width, giving the outer edge
  ## an organic, slightly-irregular look. RNG state is saved & restored
  ## so we don't clobber the caller's seed.
  if (jitter > 0 && any(is_leaf)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
                  get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed))
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(jitter_seed))
    leaf_idx <- which(is_leaf)
    rr1[leaf_idx] <- rr1[leaf_idx] +
      stats::runif(length(leaf_idx), -jitter, jitter)
  }

  ## One arc-bar polygon per node: outer arc (r1, a0 -> a1) followed by
  ## the inner arc (r0, a1 -> a0). Sampling density scales with the
  ## angular span so wide arcs stay smooth. The circular-partition
  ## convention is theta_cartesian = pi/2 - theta_partition, i.e.
  ## (x, y) = (r * sin(theta), r * cos(theta)), so the arc text-angle
  ## helper applies unchanged.
  .arc_poly <- function(i) {
    span <- theta1[i] - theta0[i]
    if (span <= 0) return(NULL)
    npt     <- max(2L, ceiling(span / (2 * pi) * 360))
    ang_out <- seq(theta0[i], theta1[i], length.out = npt)
    ang_in  <- rev(ang_out)
    data.frame(
      id    = node_names[i],
      state = unname(state_factor[i]),
      x     = c(rr1[i] * sin(ang_out), rr0[i] * sin(ang_in)),
      y     = c(rr1[i] * cos(ang_out), rr0[i] * cos(ang_in)),
      stringsAsFactors = FALSE
    )
  }
  poly_df <- do.call(rbind, lapply(seq_along(node_names), .arc_poly))

  ## Label positions: centroid angle, mid-ring radius.
  mid_ang <- (theta0 + theta1) / 2
  mid_r   <- (rr0 + rr1) / 2
  lab_df  <- data.frame(
    x     = mid_r * sin(mid_ang),
    y     = mid_r * cos(mid_ang),
    label = label_text,
    stringsAsFactors = FALSE
  )
  lab_df <- lab_df[!is.na(lab_df$label), , drop = FALSE]

  ## Radial separator segments at each depth-1 boundary, extending from
  ## the inner edge of the depth-1 ring to the outer edge of the deepest
  ## ring -- wedge dividers that "stay" out through every deeper ring.
  is_d1 <- pt_depth == 1L
  if (sum(is_d1) > 1L) {
    angles <- sort(unique(c(theta0[is_d1], theta1[is_d1]) %% (2 * pi)))
    r_in   <- ring_thickness + ring_gap        # inner edge of the depth-1 ring
    r_out  <- max(rr1, na.rm = TRUE)
    sep_df <- data.frame(
      x    = r_in  * sin(angles),
      y    = r_in  * cos(angles),
      xend = r_out * sin(angles),
      yend = r_out * cos(angles)
    )
  } else {
    sep_df <- data.frame(x = numeric(0), y = numeric(0),
                         xend = numeric(0), yend = numeric(0))
  }

  ## The per-node `nodes` frame is the plot's default data: it carries
  ## `label_text` and `arc_fraction` for introspection (each geom below
  ## supplies its own data, so this is never drawn directly).
  ggplot2::ggplot(data = nodes) +
    ggplot2::geom_polygon(
      data = poly_df,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                   fill = .data$state),
      colour    = "white",
      linewidth = gap_size_within,
      na.rm     = TRUE
    ) +
    ggplot2::geom_segment(
      data = sep_df,
      ggplot2::aes(x = .data$x, y = .data$y,
                   xend = .data$xend, yend = .data$yend),
      colour    = "white",
      linewidth = gap_size,
      lineend   = "butt"
    ) +
    ggplot2::geom_text(
      data = lab_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                   angle = .pt_arc_text_angle(.data$x, .data$y)),
      size = 2.8, colour = "grey15", fontface = "bold",
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(values = pal, name = "State",
                               na.value = "grey92",
                               na.translate = FALSE,
                               drop = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40"),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title    = "Context tree (sunburst)",
      subtitle = sprintf(
        "%d nodes, depth <= %d;  arc width = count;  fill = state;  rings = depth",
        length(tree$nodes), tree$max_depth)
    )
}
