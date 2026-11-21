# logLik(), nobs(), score_positions(), score_sequences(), perplexity(), model_fit(), mine_sequences(): the
# predictive evaluation of a tree. PST computes the same quantities: logLik() and nobs() of the fit, and
# predict() the probability of every position (decomp = TRUE), of every sequence (output = "prob") and the mean
# log loss (output = "logloss"). The rest are formulas of PST's numbers, written out below.
# (logLik, nobs and per-position prediction moved from local_testing_and_equivalence/test-equiv-PST.R.)

fit_input <- function(input) {
  context_tree(input$sequences, max_depth = input$max_depth, min_count = input$min_count,
               smoothing = list("floor", ymin = input$ymin), alphabet = input$alphabet)
}

# Inputs with ymin > 0, so that every held-out position has a probability above zero.
smoothed_inputs <- function() Filter(\(input) input$ymin > 0, ref_inputs())

test_that("logLik() and its df, and nobs(), equal PST's", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()
  trees <- lapply(inputs, fit_input)
  psts <- lapply(inputs, pst_fit_input)
  # A sum of up to 3,000 logarithms: relative 1e-10.
  expect_reference(vapply(trees, \(tree) as.numeric(logLik(tree)), numeric(1)),
                   vapply(psts, \(pst) as.numeric(pst_logLik(pst)), numeric(1)),
                   reference = "PST", covers = "logLik.transitiontrees()", quantity = "in-sample log-likelihood",
                   tolerance = 1e-10, measure = "relative")
  expect_reference(vapply(trees, \(tree) as.numeric(attr(logLik(tree), "df")), numeric(1)),
                   vapply(psts, pst_freepar, numeric(1)),
                   reference = "PST", covers = "logLik.transitiontrees()", quantity = "free parameters (df)", tolerance = 0)
  expect_reference(vapply(trees, \(tree) as.numeric(nobs(tree)), numeric(1)), vapply(psts, pst_nobs, numeric(1)),
                   reference = "PST", covers = "nobs.transitiontrees()", quantity = "number of observations", tolerance = 0)
})

test_that("score_positions() gives PST::predict(decomp = TRUE) at every position, in and out of sample", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- smoothed_inputs()
  sides <- lapply(names(inputs), \(config) {
    input <- inputs[[config]]
    half <- seq_len(nrow(input$sequences) %/% 2L)
    train <- input$sequences[half, , drop = FALSE]
    test <- input$sequences[-half, , drop = FALSE]
    tree <- context_tree(train, max_depth = input$max_depth, min_count = input$min_count,
                         smoothing = list("floor", ymin = input$ymin), alphabet = input$alphabet)
    pst <- pst_fit(train, input$alphabet, input$max_depth, input$min_count, input$ymin)
    list(ours = rbind(cbind(config = config, sample = "training", subset(score_positions(tree, newdata = train),
                                                                         select = c(sequence_id, position, predicted_prob))),
                      cbind(config = config, sample = "held out", subset(score_positions(tree, newdata = test),
                                                                         select = c(sequence_id, position, predicted_prob)))),
         theirs = rbind(cbind(config = config, sample = "training", pst_positions(pst, train)),
                        cbind(config = config, sample = "held out", pst_positions(pst, test))),
         logLik = data.frame(config = config, logLik = as.numeric(logLik(tree, newdata = test))),
         pst_logLik = data.frame(config = config, logLik = sum(log(pst_positions(pst, test)$predicted_prob))))
  })
  stack <- \(part) do.call(rbind, lapply(sides, `[[`, part))
  expect_reference(stack("ours"), stack("theirs"), reference = "PST", covers = "score_positions()",
                   quantity = "predicted probability of the observed state", tolerance = 1e-12,
                   key = c("config", "sample", "sequence_id", "position"))
  expect_reference(stack("logLik"), stack("pst_logLik"), reference = "PST", covers = "logLik.transitiontrees()",
                   quantity = "held-out log-likelihood (sum of log PST::predict(decomp = TRUE))",
                   tolerance = 1e-10, measure = "relative", key = "config")
})

test_that("score_sequences() gives log PST::predict(output = 'prob') and 2^PST::predict(output = 'logloss')", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()
  sides <- lapply(names(inputs), \(config) {
    input <- inputs[[config]]
    tree <- fit_input(input)
    pst <- pst_fit_input(input)
    states <- pst_sequences(input$sequences, input$alphabet)
    probability <- as.numeric(pst_predict(pst, states, output = "prob"))
    logloss <- as.numeric(pst_predict(pst, states, output = "logloss"))
    decomposed <- pst_predict(pst, states, output = "prob", decomp = TRUE)
    list(ours = cbind(config = config, score_sequences(tree, newdata = input$sequences)),
         theirs = data.frame(config = config, sequence_id = seq_along(probability), n_scored = rowSums(!is.na(decomposed)),
                             log_lik = log(probability), perplexity = 2^logloss))
  })
  # A product of 25 probabilities, then a log: relative 1e-10.
  expect_reference(do.call(rbind, lapply(sides, `[[`, "ours")), do.call(rbind, lapply(sides, `[[`, "theirs")),
                   reference = "PST", covers = "score_sequences()", quantity = "per-sequence score",
                   tolerance = 1e-10, measure = "relative", key = c("config", "sequence_id"))
})

test_that("perplexity() is exp(-logLik / nobs) of PST's fit, in and out of sample", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- smoothed_inputs()
  ours <- vapply(inputs, \(input) c(training = perplexity(fit_input(input)),
                                    held_out = perplexity(fit_input(input), newdata = ref_chain(30, 20, length(input$alphabet), 99L))),
                 numeric(2))
  theirs <- vapply(inputs, \(input) {
    pst <- pst_fit_input(input)
    held_out <- pst_positions(pst, ref_chain(30, 20, length(input$alphabet), 99L))
    c(training = exp(-as.numeric(pst_logLik(pst)) / pst_nobs(pst)),
      held_out = exp(-mean(log(held_out$predicted_prob))))
  }, numeric(2))
  expect_reference(ours, theirs, reference = "independent: exp(-logLik / nobs) from PST's logLik(), nobs() and predict()",
                   covers = "perplexity()", quantity = "perplexity", tolerance = 1e-10, measure = "relative")
})

test_that("model_fit() reports PST's logLik, df and nobs, and AIC, BIC and perplexity from them", {
  skip_without_reference(c("PST", "TraMineR"))
  inputs <- ref_inputs()
  ours <- do.call(rbind, lapply(names(inputs), \(config) cbind(config = config, model_fit(fit_input(inputs[[config]])))))
  pst_rows <- do.call(rbind, lapply(names(inputs), \(config) {
    pst <- pst_fit_input(inputs[[config]])
    data.frame(config = config, logLik = as.numeric(pst_logLik(pst)), df = pst_freepar(pst), nobs = pst_nobs(pst))
  }))
  expect_reference(subset(ours, select = c(config, logLik, df, nobs)), pst_rows, reference = "PST",
                   covers = "model_fit()", quantity = "logLik, df and nobs", tolerance = 1e-10, measure = "relative",
                   key = "config")
  # stats' AIC() and BIC() of PST's logLik object (its df, and its nobs added as the attribute stats reads).
  information <- do.call(rbind, lapply(names(inputs), \(config) {
    pst <- pst_fit_input(inputs[[config]])
    likelihood <- pst_logLik(pst)
    attr(likelihood, "nobs") <- pst_nobs(pst)
    data.frame(config = config, AIC = stats::AIC(likelihood), BIC = stats::BIC(likelihood))
  }))
  expect_reference(subset(ours, select = c(config, AIC, BIC)), information, reference = "stats",
                   covers = "model_fit()", quantity = "AIC and BIC of PST's log-likelihood", tolerance = 1e-10,
                   measure = "relative", key = "config")
  expect_reference(subset(ours, select = c(config, perplexity)),
                   subset(transform(pst_rows, perplexity = exp(-logLik / nobs)), select = c(config, perplexity)),
                   reference = "independent: exp(-logLik / nobs) from PST's logLik() and nobs()", covers = "model_fit()",
                   quantity = "perplexity", tolerance = 1e-10, measure = "relative", key = "config")
})

test_that("mine_sequences() returns the sequences PST's log loss ranks most and least surprising", {
  skip_without_reference(c("PST", "TraMineR"))
  input <- ref_inputs()[["c3"]]
  tree <- fit_input(input)
  pst <- pst_fit_input(input)
  logloss <- as.numeric(pst_predict(pst, pst_sequences(input$sequences, input$alphabet), output = "logloss"))
  scores <- data.frame(sequence_id = seq_along(logloss), perplexity = 2^logloss)
  surprising <- ref_head_by(scores, "perplexity", n = 8L, decreasing = TRUE)
  expected <- ref_head_by(scores, "perplexity", n = 8L, decreasing = FALSE)
  expect_reference(subset(mine_sequences(tree, newdata = input$sequences, n = 8, which = "surprising"),
                          select = c(sequence_id, perplexity)),
                   surprising, reference = "PST", covers = "mine_sequences()",
                   quantity = "the 8 most surprising sequences: perplexity (2^log loss)", tolerance = 1e-10,
                   measure = "relative", key = "sequence_id")
  expect_reference(subset(mine_sequences(tree, newdata = input$sequences, n = 8, which = "expected"),
                          select = c(sequence_id, perplexity)),
                   expected, reference = "PST", covers = "mine_sequences()",
                   quantity = "the 8 least surprising sequences: perplexity (2^log loss)", tolerance = 1e-10,
                   measure = "relative", key = "sequence_id")
})
