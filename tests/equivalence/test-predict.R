# predict(), impute_sequences(), generate_sequences(), simulate(): what a tree says comes next. The next-state
# distribution after a history is PST::query() of that history (longest suffix in the tree); modal imputation
# is PST::impute(method = "pmax"). The random draws (imputation by probability, generated sequences) are checked
# by seeded Monte Carlo against the probabilities PST gives: the frequency of each outcome over N draws must be
# within 4 standard errors of PST's probability, for each of three seeds (tolerance fixed from N and PST's
# probabilities before the draws).

fit_input <- function(input) {
  context_tree(input$sequences, max_depth = input$max_depth, min_count = input$min_count,
               smoothing = list("floor", ymin = input$ymin), alphabet = input$alphabet)
}

# Histories of every length from 1 to 12: the first 12 sequences of an input, row i cut to i states.
histories <- function(sequences) {
  cut <- sequences[seq_len(12L), seq_len(12L)]
  cut[col(cut) > row(cut)] <- NA
  cut
}

# Four standard errors of a frequency over n draws, at the largest variance among the probabilities.
four_se <- function(probabilities, n) 4 * sqrt(max(probabilities * (1 - probabilities)) / n)

test_that("predict() gives PST::query() of each history, its matched context and its modal state", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()
  sides <- lapply(names(inputs), \(config) {
    input <- inputs[[config]]
    tree <- fit_input(input)
    pst <- pst_fit_input(input)
    cut <- histories(input$sequences)
    ours <- predict(tree, newdata = cut, tidy = TRUE)
    stories <- lapply(seq_len(nrow(cut)), \(i) as.character(stats::na.omit(cut[i, ])))
    contexts <- vapply(stories, paste, character(1), collapse = " -> ")
    probs <- do.call(rbind, lapply(contexts, \(context) pst_query(pst, context)))
    theirs <- data.frame(id = seq_along(contexts), context = vapply(contexts, \(context) pst_matched(pst, context), character(1)),
                         predicted = input$alphabet[apply(probs, 1L, which.max)], probability = apply(probs, 1L, max),
                         probs, check.names = FALSE, stringsAsFactors = FALSE)
    list(ours = cbind(config = config, subset(ours, select = c("id", "context", "predicted", "probability", input$alphabet))),
         theirs = cbind(config = config, theirs))
  })
  # The states differ between configurations (2 to 4): compare per state, long.
  long <- \(tables) do.call(rbind, lapply(tables, \(table) {
    rbind(data.frame(config = table$config, id = table$id, what = "context", value = table$context),
          data.frame(config = table$config, id = table$id, what = "predicted", value = table$predicted))
  }))
  numbers <- \(tables) do.call(rbind, lapply(tables, \(table) {
    states <- intersect(LETTERS, names(table))
    do.call(rbind, lapply(c("probability", states), \(column) data.frame(config = table$config, id = table$id,
                                                                       what = column, value = table[[column]])))
  }))
  ours <- lapply(sides, `[[`, "ours")
  theirs <- lapply(sides, `[[`, "theirs")
  expect_reference(long(ours), long(theirs), reference = "PST", covers = "predict.transitiontrees()",
                   quantity = "matched context and predicted state", tolerance = 0, key = c("config", "id", "what"))
  expect_reference(numbers(ours), numbers(theirs), reference = "PST", covers = "predict.transitiontrees()",
                   quantity = "next-state probability after a history", tolerance = 1e-10, key = c("config", "id", "what"))
})

# Sequences with gaps inside and at the start (not at the end: TraMineR deletes trailing missing values, and
# impute_sequences() leaves them, so neither side fills them).
with_gaps <- function(sequences) {
  gapped <- sequences[seq_len(30L), ]
  gaps <- cbind(seq_len(30L), (seq_len(30L) * 7L) %% (ncol(gapped) - 2L) + 1L)
  gapped[gaps] <- NA
  gapped[cbind(seq_len(30L), (seq_len(30L) * 3L) %% (ncol(gapped) - 3L) + 2L)] <- NA
  gapped[1L, 1L] <- NA
  gapped
}

test_that("impute_sequences(method = 'modal') fills every gap as PST::impute(method = 'pmax')", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()
  sides <- lapply(names(inputs), \(config) {
    input <- inputs[[config]]
    gapped <- with_gaps(input$sequences)
    missing_cells <- which(is.na(gapped), arr.ind = TRUE)
    filled <- impute_sequences(fit_input(input), newdata = gapped, method = "modal")
    pst <- pst_fit_input(input)
    theirs <- as.matrix(suppressMessages(PST::impute(pst, pst_sequences(gapped, input$alphabet), method = "pmax")))
    list(ours = data.frame(config = config, row = missing_cells[, 1L], column = missing_cells[, 2L], state = filled[missing_cells]),
         theirs = data.frame(config = config, row = missing_cells[, 1L], column = missing_cells[, 2L], state = theirs[missing_cells]))
  })
  expect_reference(do.call(rbind, lapply(sides, `[[`, "ours")), do.call(rbind, lapply(sides, `[[`, "theirs")),
                   reference = "PST", covers = "impute_sequences()", quantity = "state imputed into a gap (modal)",
                   tolerance = 0, key = c("config", "row", "column"))
})

test_that("impute_sequences(method = 'prob') draws a gap's state with PST's probabilities (seeded Monte Carlo)", {
  skip_without_reference(c("PST", "TraMineR"))
  input <- ref_inputs()[["c3"]]
  tree <- fit_input(input)
  pst <- pst_fit_input(input)
  draws <- 4000L
  template <- input$sequences[1L, seq_len(8L)]
  copies <- matrix(template, nrow = draws, ncol = 8L, byrow = TRUE)
  copies[, 6L] <- NA
  probabilities <- pst_query(pst, paste(template[seq_len(5L)], collapse = " -> "))
  frequencies <- do.call(rbind, lapply(c(11L, 22L, 33L), \(seed) {
    filled <- impute_sequences(tree, newdata = copies, method = "prob", seed = seed)
    data.frame(seed = seed, state = input$alphabet, share = as.vector(table(factor(filled[, 6L], input$alphabet))) / draws)
  }))
  expected <- data.frame(seed = rep(c(11L, 22L, 33L), each = length(input$alphabet)), state = input$alphabet,
                         share = rep(unname(probabilities), 3L))
  # Monte Carlo: 4 standard errors of a share over 4,000 draws.
  expect_reference(frequencies, expected, reference = "PST", covers = "impute_sequences()",
                   quantity = "share of each state drawn into a gap", tolerance = four_se(probabilities, draws),
                   key = c("seed", "state"))
})

# The probability of every sequence of three states under PST's tree, from PST::predict(output = "prob").
three_state_probabilities <- function(pst, alphabet) {
  grid <- as.matrix(expand.grid(T1 = alphabet, T2 = alphabet, T3 = alphabet, stringsAsFactors = FALSE))
  data.frame(sequence = apply(grid, 1L, paste, collapse = "-"),
             share = as.numeric(pst_predict(pst, pst_sequences(grid, alphabet), output = "prob")))
}

shares <- function(generated, alphabet) {
  grid <- as.matrix(expand.grid(T1 = alphabet, T2 = alphabet, T3 = alphabet, stringsAsFactors = FALSE))
  all_sequences <- apply(grid, 1L, paste, collapse = "-")
  drawn <- apply(as.matrix(generated), 1L, paste, collapse = "-")
  data.frame(sequence = all_sequences, share = as.vector(table(factor(drawn, all_sequences))) / length(drawn))
}

test_that("generate_sequences() and simulate() draw sequences with PST's sequence probabilities (seeded Monte Carlo)", {
  skip_without_reference(c("PST", "TraMineR"))
  input <- ref_inputs()[["c2"]]
  tree <- fit_input(input)
  expected_one <- three_state_probabilities(pst_fit_input(input), input$alphabet)
  draws <- 20000L
  seeds <- c(5L, 6L, 7L)
  expected <- do.call(rbind, lapply(seeds, \(seed) cbind(seed = seed, expected_one)))
  generated <- do.call(rbind, lapply(seeds, \(seed) {
    set.seed(seed)
    cbind(seed = seed, shares(generate_sequences(tree, n = draws, length = 3), input$alphabet))
  }))
  simulated <- do.call(rbind, lapply(seeds, \(seed) {
    cbind(seed = seed, shares(simulate(tree, nsim = draws, seed = seed, length = 3), input$alphabet))
  }))
  # Monte Carlo: 4 standard errors of a share over 20,000 draws.
  expect_reference(generated, expected, reference = "PST", covers = "generate_sequences()",
                   quantity = "share of each sequence of three states", tolerance = four_se(expected_one$share, draws),
                   key = c("seed", "sequence"))
  expect_reference(simulated, expected, reference = "PST", covers = "simulate.transitiontrees()",
                   quantity = "share of each simulated sequence of three states",
                   tolerance = four_se(expected_one$share, draws), key = c("seed", "sequence"))
})
