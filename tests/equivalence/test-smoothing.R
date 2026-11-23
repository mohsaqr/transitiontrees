# smooth_tree(), compare_smoothing(), context_tree(smoothing = ...): the smoothers other than PST's floor rule.
# PST has none of them, so each is written out from its published definition (see ref_smoother() in
# helper-reference.R) and applied to PST's counts, shallow to deep, each context backing off to its parent's
# smoothed distribution. The floor rule PST uses (rule = "interpolate") is compared with PST in
# test-context_tree.R; the "cap" variant is written out here.
# The option "kneser_ney" is checked twice: against the rule the package documents (interpolated absolute
# discounting that backs off to the parent's smoothed distribution, Ney, Essen & Kneser 1994), and against
# interpolated Kneser-Ney itself, with continuation counts in the lower orders (Kneser & Ney 1995; Chen & Goodman
# 1999), which the package cites. The second disagrees: a declared defect, the owner's decision pending.

smoothers <- list(
  laplace = list(spec = "laplace", method = "laplace", parameter = 1),
  laplace_half = list(spec = list("laplace", alpha = 0.5), method = "laplace", parameter = 0.5),
  kneser_ney = list(spec = "kneser_ney", method = "absolute_discounting", parameter = 0.75),
  kneser_ney_half = list(spec = list("kneser_ney", discount = 0.5), method = "absolute_discounting", parameter = 0.5),
  witten_bell = list(spec = "witten_bell", method = "witten_bell", parameter = NULL),
  jelinek_mercer = list(spec = "jelinek_mercer", method = "jelinek_mercer", parameter = 0.5),
  jelinek_mercer_03 = list(spec = list("jelinek_mercer", lambda = 0.3), method = "jelinek_mercer", parameter = 0.3),
  floor_cap = list(spec = list("floor", ymin = 0.01, rule = "cap"), method = "floor_cap", parameter = 0.01))

state_probs <- function(tree) {
  do.call(rbind, lapply(tree_pathways(tree)$pathway, \(pathway) {
    probs <- query_pathway(tree, pathway, exact = TRUE)
    data.frame(pathway = pathway, state = names(probs), prob = unname(probs), stringsAsFactors = FALSE)
  }))
}

ref_state_probs <- function(smoothed) {
  do.call(rbind, lapply(names(smoothed), \(pathway) data.frame(pathway = pathway, state = names(smoothed[[pathway]]),
                                                               prob = unname(smoothed[[pathway]]), stringsAsFactors = FALSE)))
}

reference_tree <- function(pst, smoother) {
  smoothed <- ref_smooth_tree(pst, smoother$method, smoother$parameter)
  lapply(smoothed, \(probs) stats::setNames(probs, pst@alphabet))
}

test_that("every smoother, fitted or applied to a fitted tree, gives the distributions of its definition", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()
  rows <- lapply(names(inputs), \(config) {
    input <- inputs[[config]]
    pst <- pst_fit_input(input, ymin = 0)
    base <- context_tree(input$sequences, max_depth = input$max_depth, min_count = input$min_count,
                         smoothing = list("floor", ymin = 0), alphabet = input$alphabet)
    lapply(names(smoothers), \(name) {
      smoother <- smoothers[[name]]
      fitted <- context_tree(input$sequences, max_depth = input$max_depth, min_count = input$min_count,
                             smoothing = smoother$spec, alphabet = input$alphabet)
      list(fitted = cbind(config = config, smoother = name, state_probs(fitted)),
           resmoothed = cbind(config = config, smoother = name, state_probs(smooth_tree(base, smoothing = smoother$spec))),
           theirs = cbind(config = config, smoother = name, ref_state_probs(reference_tree(pst, smoother))))
    })
  })
  flat <- unlist(rows, recursive = FALSE)
  stack <- \(part) do.call(rbind, lapply(flat, `[[`, part))
  # kneser_ney here is the documented rule, interpolated absolute discounting (Ney, Essen & Kneser 1994); the
  # other definitions: Chen & Goodman 1999 (additive), Jelinek & Mercer 1980, Witten & Bell 1991.
  reference <- "independent: smoothers' definitions on PST's counts (kneser_ney = absolute discounting, Ney, Essen & Kneser 1994)"
  expect_reference(stack("fitted"), stack("theirs"), reference = reference,
                   covers = c("context_tree()", "smooth_tree()"), quantity = "smoothed next-state probability (fitted)",
                   tolerance = 1e-12, key = c("config", "smoother", "pathway", "state"))
  expect_reference(stack("resmoothed"), stack("theirs"), reference = reference, covers = "smooth_tree()",
                   quantity = "smoothed next-state probability (smooth_tree() of a fitted tree)", tolerance = 1e-12,
                   key = c("config", "smoother", "pathway", "state"))
})

test_that("smooth_tree(smoothing = 'floor') of an unsmoothed tree equals PST's ymin = 0.001 fit", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()
  ours <- do.call(rbind, lapply(names(inputs), \(config) {
    input <- inputs[[config]]
    base <- context_tree(input$sequences, max_depth = input$max_depth, min_count = input$min_count,
                         smoothing = list("floor", ymin = 0), alphabet = input$alphabet)
    cbind(config = config, state_probs(smooth_tree(base, smoothing = list("floor", ymin = 0.001))))
  }))
  theirs <- do.call(rbind, lapply(names(inputs), \(config) {
    cbind(config = config, subset(pst_state_table(pst_fit_input(inputs[[config]], ymin = 0.001)), select = -count))
  }))
  expect_reference(ours, theirs, reference = "PST", covers = "smooth_tree()",
                   quantity = "floor-smoothed probability (ymin = 0.001)", tolerance = 1e-10,
                   key = c("config", "pathway", "state"))
})

test_that("compare_smoothing() reports each smoother's nodes and in-sample perplexity", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()[c("c1", "c3", "c4", "c7")]
  methods <- c("floor", "laplace", "kneser_ney", "witten_bell", "jelinek_mercer")
  ours <- do.call(rbind, lapply(names(inputs), \(config) {
    input <- inputs[[config]]
    cbind(config = config, compare_smoothing(input$sequences, max_depth = input$max_depth,
                                             min_count = input$min_count, alphabet = input$alphabet))
  }))
  # floor: PST's own fit (ymin = 0.001, the default) and its logLik; the others: their definitions on PST's
  # counts, each position predicted by its longest context in the tree.
  floor_rows <- do.call(rbind, lapply(names(inputs), \(config) {
    pst <- pst_fit_input(inputs[[config]], ymin = 0.001)
    data.frame(config = config, smoothing = "floor", n_nodes = pst_n_nodes(pst),
               perplexity = exp(-as.numeric(pst_logLik(pst)) / pst_nobs(pst)))
  }))
  expect_reference(subset(ours, smoothing == "floor"), floor_rows,
                   reference = "independent: exp(-logLik / nobs) from PST's logLik() and nobs()",
                   covers = "compare_smoothing()", quantity = "floor: nodes and perplexity", tolerance = 1e-10,
                   measure = "relative", key = c("config", "smoothing"))
  defaults <- list(laplace = 1, kneser_ney = 0.75, witten_bell = NULL, jelinek_mercer = 0.5)
  other_rows <- do.call(rbind, lapply(names(inputs), \(config) {
    input <- inputs[[config]]
    pst <- pst_fit_input(input, ymin = 0)
    do.call(rbind, lapply(names(defaults), \(method) {
      rule <- if (method == "kneser_ney") "absolute_discounting" else method
      probs <- lapply(ref_smooth_tree(pst, rule, defaults[[method]]), \(p) stats::setNames(p, input$alphabet))
      position_probs <- ref_position_probs(input$sequences, probs, input$alphabet, input$max_depth)
      data.frame(config = config, smoothing = method, n_nodes = pst_n_nodes(pst),
                 perplexity = exp(-mean(log(position_probs))))
    }))
  }))
  # The number of nodes is the tree's, whatever the smoother: PST's count of nodes and leaves of the same fit.
  nodes <- do.call(rbind, lapply(names(inputs), \(config) {
    data.frame(config = config, smoothing = methods, n_nodes = pst_n_nodes(pst_fit_input(inputs[[config]], ymin = 0.001)))
  }))
  expect_reference(subset(ours, select = c(config, smoothing, n_nodes)), nodes, reference = "PST",
                   covers = "compare_smoothing()", quantity = "nodes of the tree per smoother", tolerance = 0,
                   key = c("config", "smoothing"))
  expect_reference(subset(ours, smoothing != "floor"), other_rows,
                   reference = "independent: perplexity of each position's longest context, smoothers' definitions on PST's counts (kneser_ney = absolute discounting)",
                   covers = "compare_smoothing()", quantity = "nodes and perplexity per smoother", tolerance = 1e-10,
                   measure = "relative", key = c("config", "smoothing"))
})

test_that("smoothing = 'kneser_ney' against interpolated Kneser-Ney with continuation counts", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()
  discounts <- c(kneser_ney = 0.75, kneser_ney_half = 0.5)
  rows <- unlist(lapply(names(inputs), \(config) {
    input <- inputs[[config]]
    pst <- pst_fit_input(input, ymin = 0)
    base <- context_tree(input$sequences, max_depth = input$max_depth, min_count = input$min_count,
                         smoothing = list("floor", ymin = 0), alphabet = input$alphabet)
    lapply(names(discounts), \(name) {
      spec <- list("kneser_ney", discount = discounts[[name]])
      fitted <- context_tree(input$sequences, max_depth = input$max_depth, min_count = input$min_count,
                             smoothing = spec, alphabet = input$alphabet)
      kneser_ney <- lapply(ref_kneser_ney_tree(pst, input$sequences, discounts[[name]]),
                           \(probs) stats::setNames(probs, input$alphabet))
      list(fitted = cbind(config = config, smoother = name, state_probs(fitted)),
           resmoothed = cbind(config = config, smoother = name, state_probs(smooth_tree(base, smoothing = spec))),
           theirs = cbind(config = config, smoother = name, ref_state_probs(kneser_ney)))
    })
  }), recursive = FALSE)
  stack <- \(part) do.call(rbind, lapply(rows, `[[`, part))
  reference <- "independent: interpolated Kneser-Ney, continuation counts (Kneser & Ney 1995; Chen & Goodman 1999)"
  defect <- paste("smoothing = 'kneser_ney' is interpolated absolute discounting, not Kneser-Ney (no continuation",
                  "counts in the lower orders); owner decision pending (implement or rename); FINDINGS.md")
  expect_reference(stack("fitted"), stack("theirs"), reference = reference,
                   covers = c("context_tree()", "smooth_tree()"),
                   quantity = "Kneser-Ney next-state probability (fitted)", tolerance = 1e-12,
                   key = c("config", "smoother", "pathway", "state"), defect = defect)
  expect_reference(stack("resmoothed"), stack("theirs"), reference = reference, covers = "smooth_tree()",
                   quantity = "Kneser-Ney next-state probability (smooth_tree() of a fitted tree)", tolerance = 1e-12,
                   key = c("config", "smoother", "pathway", "state"), defect = defect)
})
