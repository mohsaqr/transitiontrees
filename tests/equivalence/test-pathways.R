# tree_pathways(), tree_dependence(), common_pathways(), divergent_pathways(), sharp_pathways(), mine_contexts(),
# summary() of a tree: the pathway tables. Each row describes one context of the tree: its count and modal next
# state (PST's node), how far its distribution departs from that of its immediate suffix (Kullback-Leibler
# divergence in bits, entropy::KL.plugin() of PST's two distributions), and the entropies (entropy::entropy.plugin()).
# The rankings are checked as the sorted values (the order of tied rows is not part of the result) and as the rows
# they return. PST::cmine() mines the contexts after which a state is likely.

fit_input <- function(input) {
  context_tree(input$sequences, max_depth = input$max_depth, min_count = input$min_count,
               smoothing = list("floor", ymin = input$ymin), alphabet = input$alphabet)
}

# ymin > 0: every probability is positive, so every divergence is finite.
smoothed_inputs <- function() Filter(\(input) input$ymin > 0, ref_inputs())

# PST's side of the pathway table: one row per node.
pst_pathways <- function(pst) {
  nodes <- pst_node_table(pst)
  probs <- pst_prob_list(pst)
  alphabet <- pst@alphabet
  parent <- ifelse(nodes$pathway == "(start)", NA_character_, ref_parent(nodes$pathway))
  modal <- \(p) alphabet[[which.max(p)]]
  nodes$likely_next <- vapply(nodes$pathway, \(context) modal(probs[[context]]), character(1))
  nodes$next_probability <- vapply(nodes$pathway, \(context) max(probs[[context]]), numeric(1))
  nodes$divergence <- vapply(seq_along(parent), \(i) if (is.na(parent[[i]])) NA_real_ else
    ref_kl(probs[[nodes$pathway[[i]]]], probs[[parent[[i]]]], unit = "log2"), numeric(1))
  nodes$changes_prediction <- vapply(seq_along(parent), \(i) if (is.na(parent[[i]])) NA else
    nodes$likely_next[[i]] != modal(probs[[parent[[i]]]]), logical(1))
  nodes$entropy <- vapply(nodes$pathway, \(context) ref_entropy(probs[[context]], unit = "log2"), numeric(1))
  nodes$entropy_before <- vapply(parent, \(context) if (is.na(context)) NA_real_ else ref_entropy(probs[[context]], unit = "log2"),
                                 numeric(1))
  nodes$entropy_drop <- nodes$entropy_before - nodes$entropy
  nodes$likely_before <- vapply(parent, \(context) if (is.na(context)) NA_character_ else modal(probs[[context]]), character(1))
  nodes
}

by_config <- function(make) {
  inputs <- smoothed_inputs()
  do.call(rbind, lapply(names(inputs), \(config) cbind(config = config, make(inputs[[config]]), stringsAsFactors = FALSE)))
}

# The values of a column in rank order: tie-invariant.
ranked <- function(table, column) data.frame(rank = seq_len(nrow(table)), value = table[[column]])

test_that("tree_pathways() and summary() give PST's counts and modal states, and divergences from entropy", {
  skip_without_reference(c("PST", "TraMineR", "entropy"))
  ours <- by_config(\(input) tree_pathways(fit_input(input)))
  theirs <- by_config(\(input) subset(pst_pathways(pst_fit_input(input)),
                                      select = c(pathway, depth, count, likely_next, next_probability, changes_prediction)))
  expect_reference(subset(ours, select = -divergence), theirs, reference = "PST",
                   covers = c("tree_pathways()", "summary.transitiontrees()"),
                   quantity = "context: depth, count, modal next state and its probability, change of modal state",
                   tolerance = 1e-10, key = c("config", "pathway"))
  expect_reference(subset(ours, select = c(config, pathway, divergence)),
                   by_config(\(input) subset(pst_pathways(pst_fit_input(input)), select = c(pathway, divergence))),
                   reference = "entropy", covers = "tree_pathways()",
                   quantity = "KL divergence from the immediate suffix (bits), of PST's distributions", tolerance = 1e-12,
                   key = c("config", "pathway"))
  # summary() of a tree carries the same table and the size of the tree (API gap: summary() returns a list
  # with no as.data.frame(); its fields are read with $).
  sizes <- by_config(\(input) {
    overview <- summary(fit_input(input))
    data.frame(n_states = overview$n_states, max_depth = overview$max_depth, n_nodes = overview$n_nodes)
  })
  pst_sizes <- by_config(\(input) {
    pst <- pst_fit_input(input)
    data.frame(n_states = length(pst@alphabet), max_depth = as.numeric(pst_summary(pst)@depth), n_nodes = pst_n_nodes(pst))
  })
  expect_reference(sizes, pst_sizes, reference = "PST", covers = "summary.transitiontrees()",
                   quantity = "states, depth and nodes of the tree", tolerance = 0, key = "config")
})

test_that("tree_dependence() gives each context's divergence and entropies from PST's distributions", {
  skip_without_reference(c("PST", "TraMineR", "entropy"))
  ours <- by_config(\(input) tree_dependence(fit_input(input)))
  theirs <- by_config(\(input) subset(pst_pathways(pst_fit_input(input)), pathway != "(start)",
                                      select = c(pathway, depth, count, divergence, entropy, entropy_before, entropy_drop,
                                                 likely_next, likely_before, changes_prediction)))
  expect_reference(ours, theirs, reference = "entropy", covers = "tree_dependence()",
                   quantity = "divergence, entropy (bits) and modal states of a context and its suffix", tolerance = 1e-12,
                   key = c("config", "pathway"))
  in_nats <- by_config(\(input) subset(tree_dependence(fit_input(input), base = exp(1)), select = c(pathway, divergence, entropy)))
  pst_nats <- by_config(\(input) {
    table <- pst_pathways(pst_fit_input(input))
    probs <- pst_prob_list(pst_fit_input(input))
    inner <- subset(table, pathway != "(start)")
    data.frame(pathway = inner$pathway,
               divergence = vapply(inner$pathway, \(context) ref_kl(probs[[context]], probs[[ref_parent(context)]], unit = "log"), numeric(1)),
               entropy = vapply(inner$pathway, \(context) ref_entropy(probs[[context]], unit = "log"), numeric(1)))
  })
  expect_reference(in_nats, pst_nats, reference = "entropy", covers = "tree_dependence()",
                   quantity = "divergence and entropy in nats (base = e)", tolerance = 1e-12, key = c("config", "pathway"))
})

test_that("common_, divergent_ and sharp_pathways() rank PST's contexts by count, divergence and modal probability", {
  skip_without_reference(c("PST", "TraMineR", "entropy"))
  inputs <- smoothed_inputs()
  sides <- lapply(names(inputs), \(config) {
    input <- inputs[[config]]
    tree <- fit_input(input)
    table <- pst_pathways(pst_fit_input(input))
    inner <- subset(table, pathway != "(start)")
    depth_two <- subset(table, depth == 2)
    flipping <- subset(inner, changes_prediction)
    list(ours = rbind(
           cbind(config = config, list = "common, top 10, count", ranked(common_pathways(tree, top = 10), "count")),
           cbind(config = config, list = "common at depth 2, top 5, count", ranked(common_pathways(tree, top = 5, depth = 2), "count")),
           cbind(config = config, list = "divergent, top 10, divergence", ranked(divergent_pathways(tree, top = 10), "divergence")),
           cbind(config = config, list = "divergent flips only, top 10, divergence",
                 ranked(divergent_pathways(tree, top = 10, flips_only = TRUE), "divergence")),
           cbind(config = config, list = "sharp, top 10, modal probability", ranked(sharp_pathways(tree, top = 10), "next_probability"))),
         theirs = rbind(
           cbind(config = config, list = "common, top 10, count", ranked(ref_head_by(table, "count", 10L), "count")),
           cbind(config = config, list = "common at depth 2, top 5, count", ranked(ref_head_by(depth_two, "count", 5L), "count")),
           cbind(config = config, list = "divergent, top 10, divergence", ranked(ref_head_by(inner, "divergence", 10L), "divergence")),
           cbind(config = config, list = "divergent flips only, top 10, divergence",
                 ranked(ref_head_by(flipping, "divergence", 10L), "divergence")),
           cbind(config = config, list = "sharp, top 10, modal probability", ranked(ref_head_by(table, "next_probability", 10L), "next_probability"))),
         rows_ours = rbind(cbind(config = config, list = "common", subset(common_pathways(tree, top = 10), select = c(pathway, count, next_probability, divergence))),
                           cbind(config = config, list = "divergent", subset(divergent_pathways(tree, top = 10), select = c(pathway, count, next_probability, divergence))),
                           cbind(config = config, list = "sharp", subset(sharp_pathways(tree, top = 10), select = c(pathway, count, next_probability, divergence)))),
         rows_theirs = rbind(cbind(config = config, list = "common", subset(table, pathway %in% common_pathways(tree, top = 10)$pathway, select = c(pathway, count, next_probability, divergence))),
                             cbind(config = config, list = "divergent", subset(table, pathway %in% divergent_pathways(tree, top = 10)$pathway, select = c(pathway, count, next_probability, divergence))),
                             cbind(config = config, list = "sharp", subset(table, pathway %in% sharp_pathways(tree, top = 10)$pathway, select = c(pathway, count, next_probability, divergence)))))
  })
  stack <- \(part) do.call(rbind, lapply(sides, `[[`, part))
  expect_reference(stack("ours"), stack("theirs"), reference = "entropy",
                   covers = c("common_pathways()", "divergent_pathways()", "sharp_pathways()"),
                   quantity = "ranked values (PST counts and modal probabilities, entropy divergences)", tolerance = 1e-12,
                   key = c("config", "list", "rank"))
  expect_reference(stack("rows_ours"), stack("rows_theirs"), reference = "entropy",
                   covers = c("common_pathways()", "divergent_pathways()", "sharp_pathways()"),
                   quantity = "returned contexts: count, modal probability, divergence", tolerance = 1e-12,
                   key = c("config", "list", "pathway"))
})

test_that("mine_contexts() returns the contexts PST::cmine() mines", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- smoothed_inputs()
  # Cutoffs that no probability of these fits equals, so that PST's '< pmax' and the package's '<= max_prob' agree.
  sides <- lapply(names(inputs), \(config) {
    input <- inputs[[config]]
    tree <- fit_input(input)
    pst <- pst_fit_input(input)
    nodes <- pst_node_table(pst)
    lapply(input$alphabet, \(state) {
      mined <- function(found) {
        contexts <- pst_to_arrow(vapply(found, \(entry) entry@context, character(1)))
        data.frame(pathway = contexts, count = nodes$count[match(contexts, nodes$pathway)],
                   prob = vapply(found, \(entry) sum(entry@.Data[1L, state]), numeric(1)))
      }
      high <- mined(suppressMessages(PST::cmine(pst, pmin = 0.4123, state = state)))
      low <- mined(suppressMessages(PST::cmine(pst, pmax = 0.2087, state = state)))
      list(ours = rbind(ref_tag(subset(mine_contexts(tree, state = state, min_prob = 0.4123), select = c(pathway, count, prob)),
                                config = config, state = state, cut = "min_prob 0.4123"),
                        ref_tag(subset(mine_contexts(tree, state = state, max_prob = 0.2087), select = c(pathway, count, prob)),
                                config = config, state = state, cut = "max_prob 0.2087")),
           theirs = rbind(ref_tag(high, config = config, state = state, cut = "min_prob 0.4123"),
                          ref_tag(low, config = config, state = state, cut = "max_prob 0.2087")))
    })
  })
  flat <- unlist(sides, recursive = FALSE)
  expect_reference(do.call(rbind, lapply(flat, `[[`, "ours")), do.call(rbind, lapply(flat, `[[`, "theirs")),
                   reference = "PST", covers = "mine_contexts()", quantity = "mined context: count and probability of the state",
                   tolerance = 1e-10, key = c("config", "state", "cut", "pathway"))
})
