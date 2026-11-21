# prepare_input(): a long event log to one sequence per row. tna::prepare_data() (Saqr et al.) does the same
# reshaping, including sessions split where the gap between two events of an actor exceeds time_threshold. Both
# are driven over every combination of actor (none, one column, two columns), time (none, timestamp), order
# (none, each order column), action (each state column) and threshold, on the package's two long datasets. The
# contract is the same set of sequences: tna sorts its session ids by their native type, prepare_input() as text,
# so the rows can come in a different order, which does not change a fitted tree.
# (Moved from local_testing_and_equivalence/test-equiv-tna-prepare.R.)

# tna::prepare_data() takes bare column names: the call is built from symbols (two actor columns as c(a, b)).
tna_prepare <- function(data, actor, time, action, order, time_threshold) {
  arguments <- list(data = data, action = as.name(action), time_threshold = time_threshold)
  if (!is.null(actor)) arguments$actor <- if (length(actor) > 1L) as.call(c(as.name("c"), lapply(actor, as.name))) else as.name(actor)
  if (!is.null(time)) arguments$time <- as.name(time)
  if (!is.null(order)) arguments$order <- as.name(order)
  suppressMessages(do.call(tna::prepare_data, arguments))$sequence_data
}

# A wide table as its multiset of rows, compactly: each row's fingerprint sum_j code(state_j) * sqrt(prime_j)
# (code 0 for a missing value; the square roots of distinct primes are linearly independent over the rationals,
# so two different rows of codes do not share a fingerprint), sorted. Equal sorted fingerprints mean the same
# rows, whatever their order. Beside it, the table's width and number of rows.
primes <- Filter(\(n) all(n %% seq_len(floor(sqrt(n)))[-1L] != 0), 2:3000)
row_multiset <- function(wide, alphabet) {
  cells <- unname(as.matrix(as.data.frame(wide)))
  codes <- matrix(match(cells, alphabet, nomatch = 0L), nrow(cells))
  fingerprints <- sort(as.vector(codes %*% sqrt(primes[seq_len(ncol(codes))])))
  list(shape = data.frame(what = c("columns", "rows"), value = c(ncol(cells), nrow(cells))),
       rows = data.frame(rank = seq_along(fingerprints), fingerprint = fingerprints))
}

label <- function(x) if (is.null(x)) "none" else paste(x, collapse = "+")

combinations <- function(actors, times, orders, actions, thresholds) {
  grid <- expand.grid(actor = seq_along(actors), time = seq_along(times), order = seq_along(orders),
                      action = seq_along(actions), threshold = seq_along(thresholds))
  # The threshold only matters with a time column: without one, only the first threshold.
  grid <- subset(grid, !(vapply(time, \(i) is.null(times[[i]]), logical(1)) & threshold > 1L))
  lapply(seq_len(nrow(grid)), \(i) list(actor = actors[[grid$actor[[i]]]], time = times[[grid$time[[i]]]],
                                        order = orders[[grid$order[[i]]]], action = actions[[grid$action[[i]]]],
                                        threshold = thresholds[[grid$threshold[[i]]]]))
}

test_that("prepare_input() gives the sequences of tna::prepare_data() over every combination of arguments", {
  skip_without_reference("tna")
  specs <- list(
    group_regulation_long = list(data = group_regulation_long,
                                 combinations = combinations(list(NULL, "Actor", c("Actor", "Group"), c("Actor", "Course")),
                                                             list(NULL, "Time"), list(NULL), list("Action"),
                                                             c(60, 300, 900, 1800, 3600, 7200, 86400, 1e9))),
    ai_long = list(data = ai_long,
                   combinations = combinations(list(NULL, "project", "session_id", c("project", "session_id")),
                                               list(NULL, "timestamp"), list(NULL, "code_order", "order_in_session"),
                                               list("code", "cluster"), c(60, 300, 900, 1800, 3600, 7200, 1e9))))
  sides <- unlist(lapply(names(specs), \(dataset) lapply(specs[[dataset]]$combinations, \(combination) {
    name <- sprintf("%s: actor=%s time=%s order=%s action=%s threshold=%g", dataset, label(combination$actor),
                    label(combination$time), label(combination$order), combination$action, combination$threshold)
    ours <- prepare_input(specs[[dataset]]$data, actor = combination$actor, time = combination$time,
                          action = combination$action, order = combination$order, time_threshold = combination$threshold)
    theirs <- tna_prepare(specs[[dataset]]$data, combination$actor, combination$time, combination$action,
                          combination$order, combination$threshold)
    alphabet <- sort(unique(as.character(specs[[dataset]]$data[[combination$action]])))
    list(ours = row_multiset(ours, alphabet), theirs = row_multiset(theirs, alphabet), name = name)
  })), recursive = FALSE)
  stack <- \(side, part) do.call(rbind, lapply(sides, \(entry) ref_tag(entry[[side]][[part]], combination = entry$name)))
  expect_reference(stack("ours", "shape"), stack("theirs", "shape"), reference = "tna", covers = "prepare_input()",
                   quantity = "columns and rows of the sequence table", tolerance = 0, key = c("combination", "what"))
  expect_reference(stack("ours", "rows"), stack("theirs", "rows"), reference = "tna", covers = "prepare_input()",
                   quantity = "sorted row fingerprints (the set of sequences)", tolerance = 0, key = c("combination", "rank"))
})
