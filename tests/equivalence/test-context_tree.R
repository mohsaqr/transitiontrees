# context_tree(): the fitted variable-order tree. PST::pstree() (Gabadinho & Ritschard 2016) fits the same
# model on the same sequences: the same set of contexts, the same count of every next state after each context,
# and, under the floor smoother (PST's ymin), the same probabilities. markovchain's MLE is the order-1 case.
# (Moved from local_testing_and_equivalence/test-equiv-PST.R and test-equiv-markovchain.R.)

# The package's side, one row per context and next state. API gap: no exported accessor returns the per-state
# counts of a context, so they are read from the fitted object's nodes (FINDINGS.md); the probabilities come
# from query_pathway() and the contexts from tree_pathways().
tt_state_table <- function(tree) {
  contexts <- tree_pathways(tree)$pathway
  do.call(rbind, lapply(contexts, \(pathway) {
    key <- if (identical(pathway, "(start)")) "<root>" else pathway
    probs <- query_pathway(tree, pathway, exact = TRUE)
    data.frame(pathway = pathway, state = names(probs), count = as.numeric(tree$nodes[[key]]$counts),
               prob = unname(probs), stringsAsFactors = FALSE)
  }))
}

fit_input <- function(input) {
  context_tree(input$sequences, max_depth = input$max_depth, min_count = input$min_count,
               smoothing = list("floor", ymin = input$ymin), alphabet = input$alphabet)
}

by_config <- function(tables) {
  do.call(rbind, Map(\(table, config) cbind(config = config, table, stringsAsFactors = FALSE), tables, names(tables)))
}

test_that("the contexts, their counts and their probabilities equal PST::pstree() over eight configurations", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()
  ours <- lapply(inputs, fit_input)
  theirs <- lapply(inputs, pst_fit_input)
  # Counts are integers: exact.
  expect_reference(by_config(lapply(ours, \(tree) subset(tt_state_table(tree), select = -prob))),
                   by_config(lapply(theirs, \(pst) subset(pst_state_table(pst), select = -prob))),
                   reference = "PST", covers = "context_tree()", quantity = "count of the next state after a context",
                   tolerance = 0, key = c("config", "pathway", "state"))
  expect_reference(by_config(lapply(ours, \(tree) subset(tt_state_table(tree), select = -count))),
                   by_config(lapply(theirs, \(pst) subset(pst_state_table(pst), select = -count))),
                   reference = "PST", covers = c("context_tree()", "query_pathway()"),
                   quantity = "probability of the next state after a context (floor smoother)", tolerance = 1e-10,
                   key = c("config", "pathway", "state"))
})

test_that("as.data.frame() of a tree lists PST's nodes with their depth and count", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()
  ours <- lapply(inputs, \(input) subset(as.data.frame(fit_input(input)), select = c(pathway, depth, count)))
  theirs <- lapply(inputs, \(input) pst_node_table(pst_fit_input(input)))
  # Rows must match both ways: the node set is the same, then each node's depth and count.
  expect_reference(by_config(ours), by_config(theirs), reference = "PST",
                   covers = c("context_tree()", "as.data.frame.transitiontrees()", "tree_pathways()"),
                   quantity = "depth and count of each node", tolerance = 0, key = c("config", "pathway"))
  expect_reference(vapply(inputs, \(input) n_nodes(fit_input(input)), integer(1)),
                   vapply(inputs, \(input) pst_n_nodes(pst_fit_input(input)), numeric(1)),
                   reference = "PST", covers = "n_nodes()", quantity = "number of nodes (PST: internal nodes + leaves)",
                   tolerance = 0)
})

test_that("a tree of the package's engagement data equals PST::pstree() of it", {
  skip_without_reference(c("PST", "TraMineR"))
  states <- sort(unique(stats::na.omit(unlist(engagement, use.names = FALSE))))
  tree <- context_tree(engagement, max_depth = 3, min_count = 10, smoothing = list("floor", ymin = 0.001))
  pst <- pst_fit(as.matrix(engagement), states, max_depth = 3, min_count = 10, ymin = 0.001)
  expect_reference(subset(tt_state_table(tree), select = -prob), subset(pst_state_table(pst), select = -prob),
                   reference = "PST", covers = "context_tree()", quantity = "engagement: count of the next state",
                   tolerance = 0, key = c("pathway", "state"))
  expect_reference(subset(tt_state_table(tree), select = -count), subset(pst_state_table(pst), select = -count),
                   reference = "PST", covers = "context_tree()", quantity = "engagement: probability of the next state",
                   tolerance = 1e-10, key = c("pathway", "state"))
})

test_that("weighted sequences give PST's weighted counts", {
  skip_without_reference(c("PST", "TraMineR"))
  input <- ref_inputs()[["c2"]]
  weights <- rep(c(1, 2, 1.5, 3), length.out = nrow(input$sequences))
  tree <- context_tree(input$sequences, max_depth = input$max_depth, min_count = 1,
                       smoothing = list("floor", ymin = 0), alphabet = input$alphabet, weights = weights)
  pst <- pst_fit(input$sequences, input$alphabet, input$max_depth, min_count = 1, ymin = 0, weights = weights)
  # Weighted counts are sums of the weights 1 to 3 in steps of 0.5: exact in binary.
  expect_reference(subset(tt_state_table(tree), select = -prob), subset(pst_state_table(pst), select = -prob),
                   reference = "PST", covers = "context_tree()", quantity = "weighted count of the next state",
                   tolerance = 1e-12, key = c("pathway", "state"))
  expect_reference(subset(tt_state_table(tree), select = -count), subset(pst_state_table(pst), select = -count),
                   reference = "PST", covers = "context_tree()", quantity = "weighted probability of the next state",
                   tolerance = 1e-12, key = c("pathway", "state"))
})

test_that("an order-1 tree without smoothing equals markovchain's MLE transition matrix", {
  skip_without_reference("markovchain")
  inputs <- lapply(list(c(2L, 60L, 20L, 11L), c(3L, 80L, 25L, 22L), c(4L, 120L, 30L, 33L)),
                   \(cfg) list(k = cfg[[1]], sequences = ref_chain(cfg[[2]], cfg[[3]], cfg[[1]], cfg[[4]])))
  per_input <- lapply(seq_along(inputs), \(i) {
    states <- LETTERS[seq_len(inputs[[i]]$k)]
    rows <- lapply(seq_len(nrow(inputs[[i]]$sequences)), \(r) inputs[[i]]$sequences[r, ])
    tree <- context_tree(inputs[[i]]$sequences, max_depth = 1, min_count = 1, smoothing = list("floor", ymin = 0),
                         alphabet = states)
    ours <- do.call(rbind, lapply(states, \(from) data.frame(input = i, from = from, to = states,
                                                              prob = unname(query_pathway(tree, from)))))
    estimate <- markovchain::markovchainFit(rows, method = "mle")$estimate
    matrix_mle <- methods::as(estimate, "matrix")[states, states]
    theirs <- data.frame(input = i, from = rep(states, each = length(states)), to = rep(states, length(states)),
                         prob = as.vector(t(matrix_mle)))
    list(ours = ours, theirs = theirs)
  })
  expect_reference(do.call(rbind, lapply(per_input, `[[`, "ours")), do.call(rbind, lapply(per_input, `[[`, "theirs")),
                   reference = "markovchain", covers = c("context_tree()", "query_pathway()"),
                   quantity = "order-1 transition probability", tolerance = 1e-12, key = c("input", "from", "to"))
})
