# query_pathway(), pathway_exists(), subtree(), state_frequencies(): reading a fitted tree. PST::query() returns
# the next-state distribution of a context, falling back to the longest suffix that is a node when the context
# itself is not; PST's node names say which contexts exist; TraMineR::seqstatf() counts the states.
# (The query comparison moved from local_testing_and_equivalence/test-equiv-PST.R.)

fit_input <- function(input) {
  context_tree(input$sequences, max_depth = input$max_depth, min_count = input$min_count,
               smoothing = list("floor", ymin = input$ymin), alphabet = input$alphabet)
}

# Contexts to ask about: every node, and each node lengthened by one older state (deeper than the tree or
# below min_count, so that both sides fall back to a suffix).
ask_contexts <- function(nodes, alphabet) {
  real <- setdiff(nodes, "(start)")
  longer <- as.vector(outer(alphabet, real, \(state, context) paste(state, context, sep = " -> ")))
  unique(c(real, longer))
}

test_that("query_pathway() equals PST::query() for nodes and for contexts that fall back to a suffix", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()
  sides <- lapply(names(inputs), \(config) {
    input <- inputs[[config]]
    tree <- fit_input(input)
    pst <- pst_fit_input(input)
    contexts <- ask_contexts(tree_pathways(tree)$pathway, input$alphabet)
    ours <- do.call(rbind, lapply(contexts, \(context) {
      probs <- query_pathway(tree, context)
      data.frame(config = config, pathway = context, state = names(probs), prob = unname(probs))
    }))
    theirs <- do.call(rbind, lapply(contexts, \(context) {
      probs <- pst_query(pst, context)
      data.frame(config = config, pathway = context, state = names(probs), prob = unname(probs))
    }))
    one_state <- data.frame(config = config, pathway = contexts,
                            prob = vapply(contexts, \(context) query_pathway(tree, context, next_state = input$alphabet[[1L]]),
                                          numeric(1), USE.NAMES = FALSE))
    their_state <- data.frame(config = config, pathway = contexts,
                              prob = vapply(contexts, \(context) pst_query(pst, context)[[input$alphabet[[1L]]]],
                                            numeric(1), USE.NAMES = FALSE))
    list(ours = ours, theirs = theirs, one_state = one_state, their_state = their_state)
  })
  stack <- \(part) do.call(rbind, lapply(sides, `[[`, part))
  expect_reference(stack("ours"), stack("theirs"), reference = "PST", covers = "query_pathway()",
                   quantity = "next-state distribution of a context", tolerance = 1e-10,
                   key = c("config", "pathway", "state"))
  expect_reference(stack("one_state"), stack("their_state"), reference = "PST", covers = "query_pathway()",
                   quantity = "probability of the first state (next_state =)", tolerance = 1e-10,
                   key = c("config", "pathway"))
})

test_that("pathway_exists() says yes exactly for PST's nodes", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()
  sides <- lapply(names(inputs), \(config) {
    input <- inputs[[config]]
    tree <- fit_input(input)
    known <- pst_node_table(pst_fit_input(input))$pathway
    contexts <- ask_contexts(tree_pathways(tree)$pathway, input$alphabet)
    list(ours = data.frame(config = config, pathway = contexts,
                           exists = vapply(contexts, \(context) pathway_exists(tree, context), logical(1), USE.NAMES = FALSE)),
         theirs = data.frame(config = config, pathway = contexts, exists = contexts %in% known))
  })
  expect_reference(do.call(rbind, lapply(sides, `[[`, "ours")), do.call(rbind, lapply(sides, `[[`, "theirs")),
                   reference = "PST", covers = "pathway_exists()", quantity = "context is a node of the tree",
                   tolerance = 0, key = c("config", "pathway"))
})

test_that("subtree() holds the context and every PST node that extends it into the past", {
  skip_without_reference(c("PST", "TraMineR"))
  input <- ref_inputs()[["c3"]]
  tree <- fit_input(input)
  pst <- pst_fit_input(input)
  pst_nodes <- pst_node_table(pst)
  sides <- lapply(c("A", "B -> C"), \(root) {
    ours <- subset(as.data.frame(subtree(tree, root)), select = c(pathway, depth, count, next_probability))
    # The nodes of the subtree: the context itself and the contexts that end with it.
    below <- subset(pst_nodes, pathway == root | endsWith(pathway, paste0(" -> ", root)))
    below$next_probability <- vapply(below$pathway, \(context) max(pst_query(pst, context)), numeric(1))
    list(ours = cbind(root = root, ours), theirs = cbind(root = root, below))
  })
  expect_reference(do.call(rbind, lapply(sides, `[[`, "ours")), do.call(rbind, lapply(sides, `[[`, "theirs")),
                   reference = "PST", covers = "subtree()", quantity = "nodes of the subtree: depth, count, modal probability",
                   tolerance = 1e-10, key = c("root", "pathway"))
})

test_that("state_frequencies() equals TraMineR::seqstatf() and PST's root counts", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()
  ours <- do.call(rbind, lapply(names(inputs), \(config) cbind(config = config, state_frequencies(fit_input(inputs[[config]])))))
  theirs <- do.call(rbind, lapply(names(inputs), \(config) {
    frequencies <- TraMineR::seqstatf(pst_sequences(inputs[[config]]$sequences, inputs[[config]]$alphabet))
    data.frame(config = config, state = rownames(frequencies), count = frequencies$Freq,
               proportion = frequencies$Percent / 100)
  }))
  expect_reference(ours, theirs, reference = "TraMineR", covers = "state_frequencies()",
                   quantity = "state count and proportion", tolerance = 1e-12, key = c("config", "state"))
  pst_root <- do.call(rbind, lapply(names(inputs), \(config) {
    root <- subset(pst_state_table(pst_fit_input(inputs[[config]])), pathway == "(start)")
    data.frame(config = config, state = root$state, count = root$count)
  }))
  expect_reference(subset(ours, select = c(config, state, count)), pst_root, reference = "PST",
                   covers = "state_frequencies()", quantity = "state count (PST root node)", tolerance = 0,
                   key = c("config", "state"))
})
