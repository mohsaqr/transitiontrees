# ---- Core: context_tree() and friends ----

#' Internal sentinel for the root context. This value is the pinned
#' object key (parity contract) and never changes.
#' @noRd
.ROOT <- "<root>"

#' User-facing display label for the root context, shown in pathway
#' tables, plots, print methods, and accepted by \code{query_pathway()}.
#' Distinct from the internal sentinel \code{.ROOT}.
#' @noRd
.ROOT_LABEL <- "(start)"

#' Is this a supported network/transition model object?
#'
#' transitiontrees can take fitted transition/sequence-network objects
#' directly. Detected by known class, \emph{or} structurally: any
#' list-like object carrying a 2-D sequence frame in a
#' \code{$data}/\code{$sequences}/\code{$seqdata} slot. The structural
#' arm means a new object that follows the same convention works with
#' no code change here. Whether it can actually be fitted still depends
#' on a usable sequence frame being present -- see
#' \code{.ct_unwrap_netobject()}.
#' @noRd
.ct_is_netobject <- function(x) {
  if (inherits(x, c("netobject", "cograph_network", "tna", "tna_data")))
    return(TRUE)
  ## Structural duck-typing for any other family model: a list (but
  ## not a bare data.frame) exposing a rectangular sequence slot.
  if (is.list(x) && !is.data.frame(x)) {
    ## Loop with early return over three or four candidate slot names.
    for (nm in c("data", "sequences", "seqdata", "sequence_data")) {
      v <- x[[nm]]
      if (!is.null(v) && length(dim(v)) == 2L &&
          NROW(v) > 0L && NCOL(v) > 0L)
        return(TRUE)
    }
  }
  FALSE
}

#' Unwrap a network/transition object into its sequence frame.
#'
#' Network/transition objects optionally carry a
#' \code{$data} slot: the wide, trailing-NA-padded character data.frame
#' \code{context_tree()} consumes -- one row per session, one column per
#' step. \code{$nodes$label} holds the canonical alphabet (the network's
#' node set), surfaced so the fitted tree's symbol set matches the
#' network even when a code never appears within a counted window.
#'
#' The sequence frame is extracted "as is" by scanning every place a
#' network object is known to keep it -- the \code{$data}
#' slot first (the documented handoff), then other plausible
#' slots/attributes (\code{$sequences}, \code{$seqdata}, an embedded
#' netobject, or a \code{"data"}/\code{"sequences"} attribute) -- so a
#' caller can hand over any such object without knowing where the
#' upstream stashed it.
#'
#' A sequence slot may hold an \strong{integer-coded} frame plus a
#' \code{$nodes} id/label table (a common way such objects store
#' sequences). The
#' "reject numeric matrices" rule is about disambiguating a
#' \emph{top-level} square transition matrix -- it does \strong{not}
#' apply here, because a network object's sequence slot is
#' contractually sequences (its transition matrix lives in
#' \code{$weights}). So a numeric frame found in the slot is decoded,
#' not rejected: integer codes are mapped through \code{$nodes}
#' (\code{label[match(code, id)]}) when a label table exists, else
#' cast to character; \code{NA} is preserved as end-of-sequence.
#'
#' Only a \emph{pure graph} projection -- nodes/edges/weights with the
#' sequences provably absent everywhere (e.g. an aggregated transition
#' network whose constructor nulls the
#' data slot) -- cannot yield a variable-order tree: the original
#' sequences are unrecoverable from edge weights. Only then does this
#' error, with explicit guidance, rather than fabricating sequences.
#' @return list(data = <data.frame>, alphabet = <character or NULL>).
#' @noRd
.ct_seq_frame_ok <- function(d) {
  ## Any non-empty rectangular frame in a contractual sequence slot
  ## counts — numeric/integer included (states may be integer-coded).
  (is.data.frame(d) || is.matrix(d)) &&
    length(dim(d)) == 2L && nrow(d) > 0L && ncol(d) > 0L
}

#' Resolve a code -> label decoder across the supported conventions.
#'
#' Different objects disagree on where the state-label map lives:
#' some use a \code{$nodes} id/label table; others
#' use a positional \code{$labels} vector (code k = state k) and also
#' stamp \code{attr(data, "labels")} / \code{"alphabet")}. We try them
#' in that order and fall back to no decoding.
#' @return list(ids = <codes> or NULL, labels = <character> or NULL).
#'   \code{ids = NULL} with non-NULL \code{labels} means positional.
#' @noRd
.ct_resolve_labels <- function(x, hit) {
  nodes <- x$nodes
  if (is.data.frame(nodes) && !is.null(nodes$label) && !is.null(nodes$id))
    return(list(ids = nodes$id, labels = as.character(nodes$label)))
  ## Loop with early return: the first label source that matches wins.
  for (src in list(x$labels, x$alphabet,
                   attr(hit, "labels"), attr(hit, "alphabet"),
                   attr(hit, "levels")))
    if (is.character(src) && length(src) > 0L)
      return(list(ids = NULL, labels = src))   # positional
  list(ids = NULL, labels = NULL)
}

#' Decode a (possibly integer-coded) sequence frame to a character
#' data.frame, remapping codes through whatever label map the family
#' object carries.
#' @noRd
.ct_decode_seq_frame <- function(x, hit) {
  m <- as.matrix(hit)
  if (is.numeric(m)) {
    lk  <- .ct_resolve_labels(x, hit)
    if (!is.null(lk$labels)) {
      codes <- as.vector(m)
      if (is.null(lk$ids)) {
        ## Positional: code k = state k, 1-based (a common convention).
        ## A 0 (or any out-of-range code) would index labels[0] / past
        ## the end and silently drop or recycle, corrupting the frame.
        bad <- !is.na(codes) &
          (codes < 1L | codes > length(lk$labels) | codes != floor(codes))
        if (any(bad))
          .tt_abort("Integer-coded sequence frame has positional state ",
               "code(s) outside 1..", length(lk$labels), " (the ",
               "label set): ", paste(sort(unique(codes[bad])),
                                      collapse = ", "), ". transitiontrees ",
               "reads positional $labels as 1-based (code k = state ",
               "k); recode to 1-based or supply a $nodes id/label ",
               "table.")
        idx <- codes
      } else {
        idx <- match(codes, lk$ids)                  # id table
        bad <- !is.na(codes) & is.na(idx)
        if (any(bad))
          .tt_abort("Integer-coded sequence frame has code(s) absent from ",
               "the $nodes id table: ",
               paste(sort(unique(codes[bad])), collapse = ", "), ".")
      }
      out <- matrix(lk$labels[idx], nrow = nrow(m), ncol = ncol(m),
                    dimnames = dimnames(m))
    } else {
      out <- m
      storage.mode(out) <- "character"               # 1,2 -> "1","2"
    }
  } else {
    out <- m
    storage.mode(out) <- "character"                 # already char
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

.ct_unwrap_netobject <- function(x) {
  cls <- paste(class(x), collapse = "/")

  ## Candidate sequence locations, in priority order. $data is the
  ## documented handoff; the rest are defensive so the family object
  ## is accepted "as is" wherever the upstream kept the sequences.
  ## `sequence_data` is where tna::prepare_data() keeps the wide frame.
  cand <- list(x$data, x$sequences, x$seqdata, x$sequence_data,
               attr(x, "data"), attr(x, "sequences"))
  ## An embedded netobject (some wrappers nest the original).
  emb <- Filter(function(z) inherits(z, "netobject"), x)
  if (length(emb)) cand <- c(cand, lapply(emb, `[[`, "data"))

  hit <- Find(.ct_seq_frame_ok, cand)
  if (is.null(hit))
    .tt_abort("This ", cls, " object carries no sequence data anywhere ",
         "(no usable $data / $sequences / $seqdata / embedded ",
         "netobject). It looks like a pure graph - an aggregated ",
         "transition network - and transitiontrees fits on raw sequences, ",
         "not on aggregated transitions: the original sequences ",
         "cannot be recovered from edge weights (the same reason ",
         "numeric transition matrices are rejected). Pass a ",
         "sequence-bearing object instead, e.g. a wide-format ",
         "transition/network object that still carries its raw ",
         "sequences, or the original wide sequence data.frame.")

  d <- .ct_decode_seq_frame(x, hit)
  ## Default alphabet = the object's declared node/label set, so the
  ## fitted tree's symbol space matches the source model.
  alpha <- .ct_resolve_labels(x, hit)$labels
  list(data = d, alphabet = alpha)
}

#' Is this a grouped network/transition object (a named list of fittable objects)?
#'
#' Recognizes the \code{netobject_group} and \code{group_tna} classes,
#' plus any \emph{named} list whose every
#' element is itself a single fittable network/transition object. A bare
#' ragged list of character vectors is
#' \strong{not} a group (its elements are not such objects), and a
#' single network object is not a group (it is detected as a single).
#' @noRd
.ct_is_group <- function(x) {
  if (inherits(x, c("netobject_group", "group_tna", "transitiontrees_group")))
    return(TRUE)
  if (is.list(x) && !is.data.frame(x) && !.ct_is_netobject(x) &&
      length(x) > 0L) {
    nm <- names(x)
    if (!is.null(nm) && all(nzchar(nm)) &&
        all(vapply(x, .ct_is_netobject, logical(1))))
      return(TRUE)
  }
  FALSE
}

#' Normalize a grouped object to a plain named list of elements.
#' @noRd
.ct_group_elements <- function(x) {
  els <- unclass(x)
  attr(els, "group") <- NULL
  if (is.null(names(els)) || any(!nzchar(names(els))))
    names(els) <- paste0("group", seq_along(els))
  els
}

#' Shared alphabet across a group's elements: the explicit \code{alphabet}
#' if supplied, else the union of each element's declared symbol set.
#' @noRd
.ct_group_alphabet <- function(elements, alphabet) {
  if (!is.null(alphabet)) return(alphabet)
  als <- lapply(elements, function(e)
    if (.ct_is_netobject(e)) .ct_unwrap_netobject(e)$alphabet else NULL)
  als <- sort(unique(unlist(als, use.names = FALSE)))
  if (length(als)) als else NULL
}

#' Split one dataset into a named list of subsets by a grouping factor.
#'
#' \code{group} is either a character scalar naming a column of a network
#' object's \code{$metadata}, or a vector with one entry per input
#' sequence. Returns the subsets, the resolved grouping-variable name,
#' and the overall (shared) alphabet so every subtree spans the same
#' symbols.
#' @noRd
.ct_group_split_by <- function(data, group) {
  var <- NA_character_
  if (.ct_is_netobject(data)) {
    base <- .ct_unwrap_netobject(data)$data
    if (is.character(group) && length(group) == 1L) {
      meta <- data$metadata
      if (is.null(meta) || is.null(meta[[group]]))
        .tt_abort("group = \"", group, "\" is not a column of the object's ",
             "$metadata.")
      g <- meta[[group]]; var <- group
    } else g <- group
  } else {
    if (is.character(group) && length(group) == 1L)
      .tt_abort("A single-string 'group' names a $metadata column and is ",
           "only valid for a network object; for a matrix / data.frame ",
           "/ list, pass a grouping vector with one entry per sequence.")
    base <- data
    g    <- group
  }
  n <- .ct_data_n_rows(base)
  if (length(g) != n)
    .tt_abort("'group' must have one entry per input sequence (got ",
         length(g), ", expected ", n, ").")
  g   <- as.factor(g)
  idx <- split(seq_len(n), g)
  subsets <- lapply(idx, function(ix)
    if (is.list(base) && !is.data.frame(base)) base[ix]
    else base[ix, , drop = FALSE])
  ## Shared alphabet: honour a network object's declared node set (so a
  ## valid state absent from the observed rows is not lost), exactly as a
  ## single-object netobject fit does; otherwise derive it from the
  ## observed sequences.
  declared <- if (.ct_is_netobject(data))
                .ct_unwrap_netobject(data)$alphabet else NULL
  overall  <- if (!is.null(declared)) declared
              else .ct_alphabet(.ct_traj(.ct_coerce(base)))
  list(subsets = subsets, var = var, alphabet = overall, idx = idx, n = n)
}

#' Validate a top-level \code{weights} vector and split it across the
#' groups produced by \code{\link{.ct_group_split_by}} (same row index
#' as \code{$subsets}). \code{NULL} weights yield a same-length list of
#' \code{NULL}s so each per-group fit is unweighted.
#' @noRd
.ct_split_weights <- function(weights, parts) {
  if (is.null(weights)) {
    out <- vector("list", length(parts$subsets))
    names(out) <- names(parts$subsets)
    return(out)
  }
  if (!is.numeric(weights))
    .tt_abort("'weights' must be a numeric vector.")
  if (length(weights) != parts$n)
    .tt_abort("'weights' must have length equal to number of input ",
         "sequences (got ", length(weights), ", expected ",
         parts$n, ").")
  if (anyNA(weights))
    .tt_abort("'weights' must not contain NA.")
  if (any(weights < 0))
    .tt_abort("'weights' must be non-negative.")
  lapply(parts$idx, function(ix) weights[ix])
}

#' Assemble a list of fitted trees into a \code{transitiontrees_group}.
#' @noRd
.ct_as_group <- function(trees, group_var) {
  structure(trees, class = c("transitiontrees_group", "list"),
            group = if (is.null(group_var)) NA_character_ else group_var)
}

#' @noRd
.ct_coerce <- function(data) {
  if (.ct_is_netobject(data)) {
    ## network/transition object: take its $data
    ## sequence frame (errors with guidance if it is a pure graph with
    ## no sequences). context_tree() normally unwraps upstream so
    ## weights/row-count see the frame; this branch keeps the
    ## documented dispatch list self-consistent.
    return(.ct_unwrap_netobject(data)$data)
  }
  if (is.matrix(data) && !is.numeric(data)) {
    return(as.data.frame(data, stringsAsFactors = FALSE))
  }
  if (is.list(data) && !is.data.frame(data)) {
    ## list of character vectors → ragged trajectories; keep as list
    return(lapply(data, as.character))
  }
  if (is.data.frame(data)) return(data)
  .tt_abort("'data' must be a wide data.frame, character/logical matrix, ",
       "list of character vectors, or a sequence-bearing ",
       "transition/network object.")
}

#' Refuse weighted trees on surfaces that resample/permute raw sequences
#'
#' Resampling and permutation rebuild counts from \code{tree$data} without
#' the per-sequence weights, so a weighted fit would silently become an
#' unweighted analysis. Rather than return wrong numbers, these surfaces
#' error.
#' @noRd
.pt_assert_unweighted <- function(tree, what) {
  if (!is.null(tree$weights))
    .tt_abort(what, " does not support weighted trees: resampling rebuilds ",
         "counts from the raw sequences and would drop the weights. ",
         "Refit without weights, or expand each sequence by its integer ",
         "weight, before calling ", what, ".",
              class = "transitiontrees_weighted_unsupported")
  invisible(TRUE)
}

#' @noRd
.ct_data_n_rows <- function(data) {
  if (is.matrix(data) || is.data.frame(data)) return(nrow(data))
  if (is.list(data)) return(length(data))
  NA_integer_
}

#' Symbols treated as "no observation" in a sequence: R's \code{NA} and
#' the empty string, plus the two codes TraMineR uses in its state
#' sequences -- \code{"\%"} (void / end-of-sequence) and \code{"*"}
#' (missing). A column of such codes therefore never becomes a literal
#' state in the alphabet.
#' @noRd
.CT_MISSING <- c("", "%", "*")

#' Drop missing / empty cells from one sequence, leaving observed
#' states in order. Missing codes (see \code{.CT_MISSING}) and \code{NA}
#' are treated identically whether they fall at the end of a row
#' (end-of-sequence padding) or internally (a gap in the recording): the
#' gap is closed, never turned into a literal state. This is the single
#' cleaning rule shared by the list and the matrix/data.frame paths so
#' identical data yields identical trees regardless of input shape.
#' @noRd
.ct_clean_seq <- function(x) {
  x <- as.character(x)
  x[!is.na(x) & !(x %in% .CT_MISSING)]
}

#' @noRd
.ct_traj <- function(data) {
  ## Attaches an "idx" attribute mapping surviving sequences back to
  ## their original row index — needed to realign per-sequence weights
  ## after empty/short rows are dropped.
  if (is.list(data) && !is.data.frame(data)) {
    out <- lapply(data, .ct_clean_seq)
  } else {
    m <- as.matrix(data)
    storage.mode(m) <- "character"
    out <- lapply(seq_len(nrow(m)), function(i) .ct_clean_seq(m[i, ]))
  }
  keep <- which(lengths(out) > 0L)
  res <- out[keep]
  attr(res, "idx") <- keep
  res
}

#' @noRd
.ct_alphabet <- function(trajs) {
  sort(unique(unlist(trajs, use.names = FALSE)))
}

#' @noRd
.ct_count_table <- function(trajs, depth, alphabet, weights = NULL) {
  ## Vectorised k-gram counter. Returns a named list keyed by context
  ## (string with " -> " separator); the .ROOT key denotes the
  ## marginal distribution. Values are numeric vectors indexed by
  ## alphabet position (integer when weights = NULL).
  k <- depth
  na_size <- length(alphabet)
  weighted <- !is.null(weights)

  if (k == 0L) {
    pooled <- unlist(trajs, use.names = FALSE)
    if (!weighted) {
      counts <- tabulate(match(pooled, alphabet), nbins = na_size)
    } else {
      pooled_w <- unlist(mapply(function(t, w) rep(w, length(t)),
                                 trajs, weights, SIMPLIFY = FALSE),
                          use.names = FALSE)
      idx <- match(pooled, alphabet)
      counts <- vapply(seq_len(na_size),
                       function(j) sum(pooled_w[idx == j], na.rm = TRUE),
                       numeric(1))
    }
    return(setNames(list(counts), .ROOT))
  }

  per_traj <- lapply(seq_along(trajs), function(i) {
    traj <- trajs[[i]]
    L <- length(traj)
    if (L < k + 1L) return(NULL)
    starts <- seq_len(L - k)
    ctx <- vapply(starts, function(t)
      paste(traj[t:(t + k - 1L)], collapse = " -> "),
      character(1))
    list(ctx = ctx,
         nxt = traj[starts + k],
         w   = if (weighted) rep(weights[i], length(ctx)) else NULL)
  })
  ctx_vec <- unlist(lapply(per_traj, `[[`, "ctx"), use.names = FALSE)
  nxt_vec <- unlist(lapply(per_traj, `[[`, "nxt"), use.names = FALSE)
  if (weighted)
    w_vec <- unlist(lapply(per_traj, `[[`, "w"), use.names = FALSE)
  if (length(ctx_vec) == 0L) return(list())

  nxt_idx <- match(nxt_vec, alphabet)
  keep    <- !is.na(nxt_idx)
  ctx_vec <- ctx_vec[keep]
  nxt_idx <- nxt_idx[keep]
  if (weighted) w_vec <- w_vec[keep]

  unique_ctx <- unique(ctx_vec)
  out <- lapply(unique_ctx, function(ctx_str) {
    msk <- ctx_vec == ctx_str
    if (!weighted) {
      tabulate(nxt_idx[msk], nbins = na_size)
    } else {
      vapply(seq_len(na_size),
             function(j) sum(w_vec[msk & nxt_idx == j]),
             numeric(1))
    }
  })
  setNames(out, unique_ctx)
}

#' @noRd
.ct_kl <- function(p, q) {
  ## D_KL(p || q) in nats; defines 0 log 0 := 0 and returns Inf when
  ## p_i > 0 and q_i == 0.
  msk <- p > 0
  if (any(p[msk] > 0 & q[msk] == 0)) return(Inf)
  sum(p[msk] * log(p[msk] / q[msk]))
}

#' @noRd
.ct_g2 <- function(N_child, p_parent) {
  ## Likelihood-ratio G^2 statistic for the child's count vector tested
  ## against the parent's distribution. Df = (alphabet - 1) for a leaf.
  n <- sum(N_child)
  if (n == 0) return(0)
  expected <- n * p_parent
  msk <- N_child > 0
  2 * sum(N_child[msk] * log(N_child[msk] / pmax(expected[msk],
                                                  .Machine$double.eps)))
}

#' Fit one tree per group and wrap as a transitiontrees_group (internal).
#'
#' Handles both a grouped family object and an explicit \code{group =}
#' (a metadata column name or a per-sequence vector). Shares one alphabet
#' so the trees are comparable; carries an aligned \code{block} attribute
#' on each member when supplied.
#' @noRd
.ct_fit_grouped <- function(data, group, weights, block, alphabet,
                            max_depth, min_count, smoothing) {
  if (.ct_is_group(data)) {
    if (!is.null(group))
      .tt_abort("Pass either a grouped object or 'group =', not both.")
    if (!is.null(weights))
      .tt_abort("'weights' is not supported with a grouped object; fit ",
           "each group separately.")
    elements <- .ct_group_elements(data)
    alpha    <- .ct_group_alphabet(elements, alphabet)
    trees <- lapply(elements, function(e)
      context_tree(e, max_depth = max_depth, min_count = min_count,
                   smoothing = smoothing, alphabet = alpha))
    gv <- attr(data, "group")
    return(.ct_as_group(trees,
      if (is.character(gv) && length(gv) == 1L) gv else NA_character_))
  }
  parts      <- .ct_group_split_by(data, group)
  alpha      <- if (is.null(alphabet)) parts$alphabet else alphabet
  ## Split a supplied weights vector across groups by the same row index
  ## so each per-group fit stays weighted (was silently dropped).
  w_by_group <- .ct_split_weights(weights, parts)
  trees <- Map(function(d, w)
    context_tree(d, max_depth = max_depth, min_count = min_count,
                 smoothing = smoothing, alphabet = alpha, weights = w),
    parts$subsets, w_by_group)
  ## Carry the block id on each member tree, aligned to that member's own
  ## rows, so compare_groups() can use it after the same row-drop as the fit.
  if (!is.null(block)) {
    if (length(block) != parts$n)
      .tt_abort("'block' must have one entry per input sequence (got ",
           length(block), ", expected ", parts$n, ").")
    trees <- Map(function(tr, ix) { attr(tr, "block") <- block[ix]; tr },
                 trees, parts$idx)
  }
  .ct_as_group(trees, parts$var)
}

#' Validate and normalize a per-sequence weight vector (internal).
#'
#' Returns \code{NULL} (unweighted) when \code{weights} is \code{NULL} or
#' uniform; otherwise the weights reindexed to the surviving trajectories
#' (\code{idx}). Errors on the usual invalid shapes.
#' @noRd
.ct_validate_weights <- function(weights, n_rows_in, idx) {
  if (is.null(weights)) return(NULL)
  if (!is.numeric(weights))
    .tt_abort("'weights' must be a numeric vector.")
  if (length(weights) != n_rows_in)
    .tt_abort("'weights' must have length equal to number of input ",
         "sequences (got ", length(weights), ", expected ", n_rows_in,
         ").")
  if (anyNA(weights))
    .tt_abort("'weights' must not contain NA.")
  if (any(weights < 0))
    .tt_abort("'weights' must be non-negative.")
  if (sum(weights) <= 0)
    .tt_abort("'weights' must include at least one positive value; an ",
         "all-zero weight vector leaves nothing to fit.")
  weights <- weights[idx]
  ## Uniform weights == an unweighted fit (a global scale cancels in every
  ## probability), so normalise them away.
  if (length(weights) && isTRUE(all.equal(weights,
                                          rep(weights[1L], length(weights)))))
    return(NULL)
  weights
}

#' Build the node list top-down (internal).
#'
#' Returns \code{list(nodes, max_depth)}; \code{max_depth} is lowered if a
#' depth yields no node clearing \code{min_count}.
#' @noRd
.ct_build_nodes <- function(trajs, max_depth, min_count, alphabet, sm,
                            weights) {
  nodes <- list()
  for (d in seq.int(0L, max_depth)) {
    counts_d <- .ct_count_table(trajs, depth = d, alphabet = alphabet,
                                weights = weights)
    keep <- vapply(counts_d, sum, numeric(1)) >= min_count
    if (d == 0L && .ROOT %in% names(counts_d)) keep[.ROOT] <- TRUE
    counts_d <- counts_d[keep]
    if (length(counts_d) == 0L) { max_depth <- d - 1L; break }
    for (ctx in names(counts_d)) {
      counts      <- counts_d[[ctx]]
      parent_prob <- if (d == 0L) NULL else .pt_parent_prob(nodes, ctx)
      nodes[[ctx]] <- list(
        depth  = d,
        counts = counts,
        prob   = .ct_smooth_dispatch(sm, counts, parent_prob),
        n      = sum(counts))
    }
  }
  list(nodes = nodes, max_depth = max_depth)
}

#' Build the parent/child edge frame from a node list (internal).
#'
#' A node "x1 -> ... -> xk" has parent "x2 -> ... -> xk" (drop the
#' leftmost / oldest move); the root carries the \code{.ROOT} sentinel.
#' @noRd
.ct_build_edges <- function(nodes) {
  edges <- do.call(rbind, lapply(names(nodes), function(ctx) {
    if (identical(ctx, .ROOT)) return(NULL)
    parts  <- strsplit(ctx, " -> ", fixed = TRUE)[[1L]]
    parent <- if (length(parts) == 1L) .ROOT else
      paste(parts[-1L], collapse = " -> ")
    if (!parent %in% names(nodes)) return(NULL)
    data.frame(parent = parent, child = ctx, symbol = parts[1L],
               stringsAsFactors = FALSE)
  }))
  if (is.null(edges))
    edges <- data.frame(parent = character(0), child = character(0),
                        symbol = character(0), stringsAsFactors = FALSE)
  edges
}

#' Fit a variable-order transition tree
#'
#' @description
#' \code{context_tree()} fits a variable-order transition tree to a set of
#' sequences. The tree, known in the literature as a prediction suffix tree
#' (Ron, Singer & Tishby, 1996) or a variable-length Markov chain
#' (\enc{Bühlmann}{Buhlmann} & Wyner, 1999), describes how the next state of a
#' process depends on the states that precede it, while letting the amount of
#' history that matters differ from one history to another. Each node of the
#' tree is a context, or pathway: the most recent states, read from oldest to
#' newest (\code{"A -> B -> C"}). The root holds the marginal distribution of
#' states; its children condition on the last state, their children on the last
#' two states, and so on. The tree shows how much of the past is needed to
#' predict what comes next, and along which pathways a longer history changes
#' the prediction.
#'
#' It counts every context of up to \code{max_depth} states that occurs at least
#' \code{min_count} times and estimates, for each, the conditional distribution
#' of the next state, smoothed so that unobserved continuations keep a small
#' probability. Sequences may be given as a wide table or matrix (one row per
#' sequence, one column per time point), a list of character vectors, a long
#' event log (one row per event; name the state column in \code{action}), or a
#' fitted network object from \pkg{tna} or \pkg{Nestimate} that carries its
#' sequences. With \code{group}, or with a grouped network object, one tree is
#' fitted per group over a shared set of states. With the default floor
#' smoothing, the fitted counts and probabilities match those of the \pkg{PST}
#' package (Gabadinho & Ritschard, 2016). It returns a \code{transitiontrees}
#' object, or a \code{transitiontrees_group} with one tree per group.
#'
#' The tree is returned unpruned: every context that clears
#' \code{min_count} is kept. Use \code{\link{prune_tree}()} to retain only
#' the contexts whose prediction differs from that of their shorter suffix,
#' \code{\link{tune_tree}()} to choose \code{max_depth}, \code{min_count}
#' and the smoothing scheme by cross-validation, and
#' \code{\link{tree_pathways}()} or \code{as.data.frame()} to read the
#' fitted contexts as a table. Group differences are tested with
#' \code{\link{compare_groups}()}, and the stability of individual pathways
#' with \code{\link{bootstrap_pathways}()}.
#'
#' @param data Sequences: a wide \code{data.frame} or character matrix (one
#'   row per sequence, one column per time point), a list of character
#'   vectors, a long event log (with \code{action}), or a \pkg{tna} or
#'   \pkg{Nestimate} model object that carries its sequence data.
#'   \code{NA}, \code{""} and the \pkg{TraMineR} codes \code{"\%"} and
#'   \code{"*"} mark positions without an observation and are dropped.
#'   Numeric matrices are not accepted; convert them to character first.
#' @param max_depth Longest context, in states, that the tree may contain.
#'   Default 5.
#' @param min_count Minimum number of occurrences for a context to become a
#'   node. Default 5.
#' @param smoothing Smoothing scheme for the next-state probabilities: a
#'   method name or a list of the method name and its hyperparameters.
#'   Methods are \code{"floor"} (default; \code{ymin = 0.001},
#'   \code{rule = "interpolate"} or \code{"cap"}), \code{"laplace"}
#'   (\code{alpha = 1}), \code{"kneser_ney"} (\code{discount = 0.75}),
#'   \code{"witten_bell"} and \code{"jelinek_mercer"} (\code{lambda = 0.5}).
#'   See \code{\link{smooth_tree}()}.
#' @param alphabet Character vector of states. \code{NULL} (default) uses
#'   the states observed in \code{data} (for a network object, its node
#'   labels). Supply it when new data may contain states unseen here.
#' @param weights Non-negative weight per sequence, one per row or list
#'   element of \code{data}. \code{NULL} (default) fits unweighted.
#' @param group Grouping for one tree per group: a vector with one value per
#'   sequence, the name of a column of \code{data} (long format) or of a
#'   network object's metadata. \code{NULL} (default) fits a single tree.
#' @param block Optional block (stratum) identifier in the same forms as
#'   \code{group}, kept with the fitted group and used by
#'   \code{\link{compare_groups}()} for a stratified permutation test.
#'   Default \code{NULL}.
#' @param action Name of the state column of a long event log. Giving it
#'   reshapes \code{data} with \code{\link{prepare_input}()} before
#'   fitting. Default \code{NULL} (\code{data} is already in sequence form).
#' @param actor,time,order,session,time_threshold Long-format arguments
#'   passed to \code{\link{prepare_input}()}: the column identifying whose
#'   events form a sequence, the timestamp column, an ordering column, a
#'   session column, and the gap in seconds that starts a new session
#'   (default 900). Each defaults to \code{NULL}.
#'
#' @return An object of class \code{transitiontrees}, the fitted tree. Read
#'   it through its accessors: \code{as.data.frame()} or
#'   \code{\link{tree_pathways}()} for one row per context,
#'   \code{\link{state_frequencies}()} for the state counts,
#'   \code{\link{model_fit}()} for the log-likelihood, information criteria
#'   and perplexity, and \code{\link{n_nodes}()} for the number of
#'   contexts. With \code{group}, or a grouped input, an object of class
#'   \code{transitiontrees_group}: one \code{transitiontrees} per group, in
#'   group order, summarized by \code{summary()} with one row per group and
#'   flattened by \code{as.data.frame()} to one row per group and context.
#'
#' @section Conditions:
#' \code{transitiontrees_bad_input} when \code{max_depth} or
#' \code{min_count} is not a single valid number, \code{data} has an
#' unsupported type or no usable sequence, a named column is absent,
#' \code{weights}, \code{group} or \code{block} do not have one entry per
#' sequence, a network object carries no sequence data, or the smoothing
#' specification is invalid.
#'
#' @references
#' \enc{Bühlmann}{Buhlmann}, P., & Wyner, A. J. (1999). Variable length
#' Markov chains. \emph{The Annals of Statistics}, 27(2), 480-513.
#' \doi{10.1214/aos/1018031204}
#'
#' Gabadinho, A., & Ritschard, G. (2016). Analyzing state sequences with
#' probabilistic suffix trees: The PST R package. \emph{Journal of
#' Statistical Software}, 72(3). \doi{10.18637/jss.v072.i03}
#'
#' Ron, D., Singer, Y., & Tishby, N. (1996). The power of amnesia: Learning
#' probabilistic automata with variable memory length. \emph{Machine
#' Learning}, 25(2-3), 117-149. \doi{10.1023/A:1026490906255}
#'
#' @seealso \code{\link{prune_tree}()} and \code{\link{smooth_tree}()} to
#'   refine a fitted tree; \code{\link{tune_tree}()} to choose its
#'   settings; \code{\link{tree_pathways}()} to read it;
#'   \code{\link{plot.transitiontrees}()} to draw it;
#'   \code{\link{compare_groups}()} to compare groups.
#'
#' @examples
#' # Wide sequences: one learner per row, one time point per column
#' tree <- context_tree(trajectories, max_depth = 3)
#' tree
#' as.data.frame(tree)
#'
#' # A long event log, split into sessions at gaps over 15 minutes
#' context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'              action = "Action", max_depth = 2)
#'
#' # One tree per achievement group over a shared set of states
#' groups <- context_tree(group_regulation_long, actor = "Actor",
#'                        time = "Time", action = "Action",
#'                        group = "Achiever", max_depth = 2)
#' summary(groups)
#'
#' @export
context_tree <- function(data,
                         max_depth = 5L,
                         min_count = 5L,
                         smoothing = "floor",
                         alphabet  = NULL,
                         weights   = NULL,
                         group     = NULL,
                         block     = NULL,
                         actor = NULL, time = NULL, action = NULL,
                         order = NULL, session = NULL,
                         time_threshold = 900) {
  if (!is.numeric(max_depth) || length(max_depth) != 1L ||
      is.na(max_depth) || max_depth < 0)
    .tt_abort("'max_depth' must be a single non-negative integer.")
  if (!is.numeric(min_count) || length(min_count) != 1L ||
      is.na(min_count) || min_count < 1)
    .tt_abort("'min_count' must be a single integer >= 1.")

  ## ---- Long-format `group` / `block` as column names ------------------
  ## In long mode, `group` and `block` may simply name columns of `data`;
  ## prepare_input() carries them through the reshape as per-sequence
  ## metadata aligned to the wide rows (one value per sequence), so the
  ## caller never hand-builds or re-aligns a vector. This works for every
  ## reshape path (time/session sessionization included), because the
  ## alignment is produced by the reshape itself, not re-derived from the
  ## actor key.
  long_mode <- !is.null(action) && is.data.frame(data)
  is_str <- function(x) is.character(x) && length(x) == 1L
  grp_col <- if (long_mode && is_str(group)) group else NULL
  blk_col <- if (long_mode && is_str(block)) block else NULL
  for (cn in c(grp_col, blk_col))
    if (!cn %in% names(data))
      .tt_abort("'", cn, "' is not a column of 'data'.")

  ## ---- Long-format input: reshape to a wide sequence frame first ------
  ## When `action` (the state column) is named, `data` is treated as a
  ## long, one-row-per-event table and reshaped with prepare_input() --
  ## grouping by `actor`, ordering by `time`/`order`, and splitting an
  ## actor's events into sessions when the time gap exceeds
  ## `time_threshold`. Everything below then fits on the wide frame.
  data <- .ct_maybe_reshape(data, actor, time, action, order, session,
                            time_threshold,
                            meta = unique(c(grp_col, blk_col)))
  if (!is.null(grp_col)) group <- attr(data, "meta")[[grp_col]]
  if (!is.null(blk_col)) block <- attr(data, "meta")[[blk_col]]

  ## ---- Grouped fit -> a transitiontrees_group (one tree per group) ----------
  ## A grouped family object (netobject_group / group_tna / a named list
  ## of family objects) dispatches over its elements; an explicit
  ## group = (a metadata column name or a per-sequence vector) splits a
  ## single dataset. Both share one alphabet so the trees are
  ## comparable, and return a transitiontrees_group.
  if (.ct_is_group(data) || !is.null(group))
    return(.ct_fit_grouped(data, group, weights, block, alphabet,
                           max_depth, min_count, smoothing))

  if (.ct_is_netobject(data)) {
    ## Unwrap before row-count / weight detection so they see the
    ## sequence frame, not the 18-slot network list. The network's
    ## node set becomes the default alphabet unless the caller
    ## overrides it explicitly.
    no <- .ct_unwrap_netobject(data)
    data <- no$data
    if (is.null(alphabet)) alphabet <- no$alphabet
  }
  n_rows_in <- .ct_data_n_rows(data)

  data  <- .ct_coerce(data)
  trajs <- .ct_traj(data)
  if (length(trajs) == 0L)
    .tt_abort("No usable sequences after coercion.")

  weights <- .ct_validate_weights(weights, n_rows_in, attr(trajs, "idx"))

  alphabet  <- if (is.null(alphabet)) .ct_alphabet(trajs) else alphabet
  max_depth <- as.integer(max_depth)
  min_count <- as.integer(min_count)
  sm        <- .pt_resolve_smoothing(smoothing)

  ## Smooth top-down so each node's parent prob is available when the node
  ## itself is computed (kneser_ney/witten_bell/jelinek_mercer interpolate
  ## with the parent); a depth with no surviving node lowers max_depth.
  built     <- .ct_build_nodes(trajs, max_depth, min_count, alphabet, sm,
                               weights)
  nodes     <- built$nodes
  max_depth <- built$max_depth
  edges     <- .ct_build_edges(nodes)

  structure(
    list(
      nodes      = nodes,
      edges      = edges,
      alphabet   = alphabet,
      max_depth  = max_depth,
      nmin       = min_count,
      n_seq      = length(trajs),
      ## Weighted token total (= sum of root counts); equals the
      ## unweighted token count when no weights are in play, and matches
      ## the observation count logLik() scores over.
      n_obs      = sum(nodes[[.ROOT]]$counts),
      smoothing  = sm,
      pruned     = FALSE,
      pruning    = NULL,
      data       = trajs,
      weights    = weights
    ),
    class = "transitiontrees"
  )
}

#' Convert a fitted transition tree to a data frame of contexts
#'
#' @description
#' \code{as.data.frame()} returns the contexts of a fitted tree as a table with
#' one row per context, the same table that \code{\link{tree_pathways}()}
#' returns. Each row gives the context, its depth and count, the most likely
#' next state and its probability, and how much the context's prediction departs
#' from that of its shorter suffix.
#'
#' @param x A \code{transitiontrees} object.
#' @param row.names,optional Not used.
#' @param ... Passed to \code{\link{tree_pathways}()}: \code{min_count},
#'   \code{sort_by} and \code{decreasing}.
#'
#' @return A \code{data.frame} with one row per context and columns
#'   \code{pathway}, \code{depth}, \code{count}, \code{likely_next},
#'   \code{next_probability}, \code{divergence} and
#'   \code{changes_prediction}, described in \code{\link{tree_pathways}()}.
#'
#' @seealso \code{\link{tree_pathways}()}, \code{\link{common_pathways}()},
#'   \code{\link{divergent_pathways}()}, \code{\link{tree_dependence}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' as.data.frame(tree)
#' as.data.frame(tree, sort_by = "divergence")
#' @export
as.data.frame.transitiontrees <- function(x, row.names = NULL,
                                    optional = FALSE, ...) {
  tree_pathways(x, ...)
}

#' Print a group of fitted transition trees
#'
#' @description
#' \code{print()} prints one line per group of a \code{transitiontrees_group}:
#' the number of contexts, the depth of the tree, and the numbers of sequences
#' and observations it was fitted on.
#'
#' @param x A \code{transitiontrees_group} object.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.transitiontrees_group}()} for the fit of
#'   each group; \code{\link{context_tree}()}.
#'
#' @examples
#' groups <- context_tree(group_regulation_long, actor = "Actor",
#'                        time = "Time", action = "Action",
#'                        group = "Achiever", max_depth = 2)
#' print(groups)
#' @export
print.transitiontrees_group <- function(x, ...) {
  gv <- attr(x, "group")
  by <- if (!is.null(gv) && !is.na(gv)) sprintf(" by '%s'", gv) else ""
  cat(sprintf("<transitiontrees_group>  %d groups%s\n", length(x), by))
  for (nm in names(x)) {
    t <- x[[nm]]
    cat(sprintf("  %-14s %3d nodes, depth <= %d, %d seq, %d obs%s\n",
                nm, length(t$nodes), t$max_depth, t$n_seq, t$n_obs,
                if (isTRUE(t$pruned)) "  [pruned]" else ""))
  }
  invisible(x)
}

#' Summarize a group of fitted transition trees
#'
#' @description
#' \code{summary()} summarizes a \code{transitiontrees_group} with one row per
#' group: the size of the data each tree was fitted on, the size and depth of
#' the tree, and its in-sample fit (log-likelihood, number of parameters, AIC,
#' BIC and perplexity, as returned by \code{\link{model_fit}()}). The in-sample
#' fit describes each tree; to test whether the groups differ, use
#' \code{\link{compare_groups}()}.
#'
#' @param object A \code{transitiontrees_group} object.
#' @param ... Not used.
#'
#' @return An object of class \code{summary.transitiontrees_group}, a
#'   \code{data.frame} with one row per group and columns \code{group},
#'   \code{n_seq}, \code{n_obs}, \code{n_nodes}, \code{max_depth},
#'   \code{pruned}, \code{logLik}, \code{df}, \code{nobs}, \code{AIC},
#'   \code{BIC} and \code{perplexity}.
#'
#' @seealso \code{\link{context_tree}()}, \code{\link{model_fit}()},
#'   \code{\link{compare_groups}()}.
#'
#' @examples
#' groups <- context_tree(group_regulation_long, actor = "Actor",
#'                        time = "Time", action = "Action",
#'                        group = "Achiever", max_depth = 2)
#' summary(groups)
#' @export
summary.transitiontrees_group <- function(object, ...) {
  .tt_check_class(object, "transitiontrees_group", "object")
  size <- data.frame(
    group     = names(object),
    n_seq     = vapply(object, function(t) as.integer(t$n_seq), integer(1)),
    n_obs     = vapply(object, function(t) as.integer(t$n_obs), integer(1)),
    n_nodes   = vapply(object, function(t) length(t$nodes), integer(1)),
    max_depth = vapply(object, function(t) as.integer(t$max_depth),
                       integer(1)),
    pruned    = vapply(object, function(t) isTRUE(t$pruned), logical(1)),
    stringsAsFactors = FALSE, row.names = NULL)
  fit <- model_fit(object)
  out <- merge(size, fit, by = "group", sort = FALSE)
  out <- out[match(names(object), out$group), , drop = FALSE]
  rownames(out) <- NULL
  stopifnot("one summary row per group" = nrow(out) == length(object))
  class(out) <- c("summary.transitiontrees_group", "data.frame")
  out
}

#' Print the summary of a group of transition trees
#'
#' @description
#' \code{print()} prints the table returned by
#' \code{\link{summary.transitiontrees_group}()} under a one-line header, with
#' numeric columns rounded.
#'
#' @param x A \code{summary.transitiontrees_group} object.
#' @param digits Number of significant digits for numeric columns.
#'   Default 3.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.transitiontrees_group}()}.
#'
#' @examples
#' groups <- context_tree(group_regulation_long, actor = "Actor",
#'                        time = "Time", action = "Action",
#'                        group = "Achiever", max_depth = 2)
#' print(summary(groups), digits = 2)
#' @export
print.summary.transitiontrees_group <- function(x, digits = 3L, ...) {
  cat(sprintf("<transitiontrees_group summary>  %d groups\n\n", nrow(x)))
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  print.data.frame(df, row.names = FALSE, digits = digits)
  invisible(x)
}

#' Convert a group of transition trees to one data frame of contexts
#'
#' @description
#' \code{as.data.frame()} stacks the context table of each tree in a
#' \code{transitiontrees_group}, as returned by \code{\link{tree_pathways}()},
#' under a leading \code{group} column, so that the contexts of all groups can
#' be read, filtered or plotted together.
#'
#' @param x A \code{transitiontrees_group} object.
#' @param row.names,optional Not used.
#' @param ... Passed to \code{\link{tree_pathways}()}: \code{min_count},
#'   \code{sort_by} and \code{decreasing}.
#'
#' @return A \code{data.frame} with one row per group and context: a
#'   \code{group} column followed by the columns of
#'   \code{\link{tree_pathways}()}.
#'
#' @seealso \code{\link{tree_pathways}()}; \code{\link{compare_groups}()}
#'   to test group differences context by context.
#'
#' @examples
#' groups <- context_tree(group_regulation_long, actor = "Actor",
#'                        time = "Time", action = "Action",
#'                        group = "Achiever", max_depth = 2)
#' head(as.data.frame(groups, min_count = 50))
#' @export
as.data.frame.transitiontrees_group <- function(x, row.names = NULL,
                                         optional = FALSE, ...) {
  parts <- lapply(names(x), function(nm) {
    d <- tree_pathways(x[[nm]], ...)
    if (nrow(d) == 0L) return(NULL)
    cbind(group = nm, d, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  if (is.null(out))
    return(cbind(group = character(0), .pt_empty_pathways_df(),
                 stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}

#' Print a fitted transition tree
#'
#' @description
#' \code{print()} prints a header (numbers of contexts, depth and states, the
#' data the tree was fitted on, the smoothing scheme and \code{min_count}, and
#' the pruning criterion when the tree is pruned), followed by the tree drawn as
#' text. Each line of the drawing is one context: the state it adds to its
#' parent's history (the older state), its count \code{n}, and the most likely
#' next state with its probability. Long trees are cut after \code{max_lines}
#' lines; \code{as.data.frame()} gives every context.
#'
#' @param x A \code{transitiontrees} object.
#' @param max_lines Maximum number of tree lines to print. Default 25.
#' @param digits Decimal places for probabilities. Default 2.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.transitiontrees}()},
#'   \code{\link{as.data.frame.transitiontrees}()},
#'   \code{\link{plot.transitiontrees}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' print(tree, max_lines = 10)
#' @export
print.transitiontrees <- function(x, max_lines = 25L, digits = 2L, ...) {
  alpha <- x$alphabet
  sm_label <- .pt_smoothing_label(x$smoothing)
  cat(sprintf(
    "<transitiontrees>  %d nodes, depth <= %d, %d states  [%s]\n",
    length(x$nodes), x$max_depth, length(alpha),
    if (isTRUE(x$pruned)) "pruned" else "unpruned"
  ))
  cat(sprintf("  alphabet : %s\n", paste(alpha, collapse = ", ")))
  cat(sprintf("  fit on   : %d sequences, %d observations\n",
              x$n_seq, x$n_obs))
  cat(sprintf("  smoothing: %s   min_count = %d\n", sm_label, x$nmin))
  if (isTRUE(x$pruned) && !is.null(x$pruning))
    cat(sprintf("  pruned by: %s   alpha = %s\n",
                x$pruning$criterion, format(x$pruning$alpha)))
  local_root <- attr(x, "local_root")
  if (!is.null(local_root))
    cat(sprintf("  subtree of: %s\n",
                if (identical(local_root, .ROOT)) .ROOT_LABEL else local_root))
  if (length(x$nodes) == 0L) return(invisible(x))

  render <- function(parent, prefix) {
    children <- sort(x$edges$child[x$edges$parent == parent])
    for (i in seq_along(children)) {
      child  <- children[i]; last <- i == length(children)
      branch <- if (last) "`-- " else "|-- "
      info   <- x$nodes[[child]]
      modal  <- alpha[which.max(info$prob)]
      label  <- if (identical(child, .ROOT)) .ROOT_LABEL else
        strsplit(child, " -> ", fixed = TRUE)[[1L]][[1L]]
      lab <- sprintf("%-8s  n=%-5s  -> %s (%.*f)",
                     label, format(as.integer(info$n)), modal, digits,
                     max(info$prob))
      cat(prefix, branch, lab, "\n", sep = "")
      render(child, paste0(prefix, if (last) "    " else "|   "))
    }
  }
  root <- x$nodes[[.ROOT]]
  if (!is.null(root)) {
    modal <- alpha[which.max(root$prob)]
    cat(sprintf("%-8s  n=%-5s  -> %s (%.*f)\n",
                .ROOT_LABEL, format(as.integer(root$n)), modal, digits,
                max(root$prob)))
  }
  out <- utils::capture.output(render(.ROOT, ""))
  if (length(out) > max_lines) {
    cat(paste(out[seq_len(max_lines)], collapse = "\n"), "\n")
    cat(sprintf("... %d more nodes (use as.data.frame(x) or summary(x))\n",
                length(out) - max_lines))
  } else if (length(out) > 0L) {
    cat(paste(out, collapse = "\n"), "\n")
  }
  invisible(x)
}

#' Summarize a fitted transition tree
#'
#' @description
#' \code{summary()} summarizes a fitted tree by its size (numbers of states and
#' contexts, depth, pruned or not) and its table of contexts ordered by depth
#' and, within depth, by decreasing count, so that the table reads from the root
#' downward.
#'
#' @param object A \code{transitiontrees} object.
#' @param ... Not used.
#'
#' @return An object of class \code{summary.transitiontrees}, which prints
#'   a header and the first rows of the context table (columns as in
#'   \code{\link{tree_pathways}()}). For the full table use
#'   \code{as.data.frame()} on the tree; for fit statistics use
#'   \code{\link{model_fit}()}.
#'
#' @seealso \code{\link{print.summary.transitiontrees}()},
#'   \code{\link{tree_pathways}()}, \code{\link{model_fit}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' summary(tree)
#' @export
summary.transitiontrees <- function(object, ...) {
  tbl <- tree_pathways(object, min_count = 1L)
  if (nrow(tbl) > 0L) {
    tbl <- tbl[order(tbl$depth, -tbl$count), , drop = FALSE]
    rownames(tbl) <- NULL
  }
  structure(
    list(table = tbl, n_states = length(object$alphabet),
         max_depth = object$max_depth, n_nodes = length(object$nodes),
         pruned = isTRUE(object$pruned)),
    class = "summary.transitiontrees"
  )
}

#' Print the summary of a transition tree
#'
#' @description
#' \code{print()} prints the object returned by
#' \code{\link{summary.transitiontrees}()}: a one-line header followed by the
#' first \code{n} rows of the context table.
#'
#' @param x A \code{summary.transitiontrees} object.
#' @param n Number of context rows to print. Default 10.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{summary.transitiontrees}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 2)
#' print(summary(tree), n = 5)
#' @export
print.summary.transitiontrees <- function(x, n = 10L, ...) {
  cat(sprintf("<transitiontrees summary>  %d nodes, depth <= %d, %d states  [%s]\n\n",
              x$n_nodes, x$max_depth, x$n_states,
              if (x$pruned) "pruned" else "unpruned"))
  tbl <- x$table
  print(utils::head(tbl, n), row.names = FALSE)
  if (nrow(tbl) > n)
    cat(sprintf("# ... %d more rows (use as.data.frame(tree) for the full table)\n",
                nrow(tbl) - n))
  invisible(x)
}
