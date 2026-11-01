#' Student engagement trajectories
#'
#' An example set of categorical learning-engagement trajectories used in
#' the \pkg{transitiontrees} examples and in the \dQuote{Getting started},
#' \dQuote{Advanced analysis} and \dQuote{Visualization} vignettes.
#' Each row is one learner, each column a time-step, and each cell the
#' engagement state at that step. Trailing \code{NA}s mark the end of a
#' trajectory. This wide character matrix is exactly the shape
#' \code{\link{context_tree}()} consumes.
#'
#' @format A character matrix with 138 rows (learners) and 15 columns
#'   (time-steps). Three states: \code{"Active"}, \code{"Average"},
#'   \code{"Disengaged"}.
#'
#' @source Bundled example dataset.
#'
#' @examples
#' data(trajectories)
#' dim(trajectories)
#' tree <- context_tree(trajectories)
#' tree
"trajectories"

#' Collaborative-regulation events (long format)
#'
#' A long, one-row-per-event log of collaborative regulation moves, with
#' timestamps. Used to demonstrate long-format loading: reshape with
#' \code{\link{prepare_input}()} (or name \code{action =} in
#' \code{\link{context_tree}()}) to split each actor's events into
#' time-gap sessions. Bundled example dataset.
#'
#' @format A \code{data.frame} with 27533 rows and 6 columns:
#'   \describe{
#'     \item{Actor}{integer; the learner.}
#'     \item{Achiever}{character; an achievement-level covariate.}
#'     \item{Group}{numeric; the collaboration group.}
#'     \item{Course}{character; the course.}
#'     \item{Time}{POSIXct; the event timestamp.}
#'     \item{Action}{character; the regulation move (the state).}
#'   }
#' @source Bundled example dataset.
#' @examples
#' data(group_regulation_long)
#' context_tree(group_regulation_long, actor = "Actor", time = "Time",
#'              action = "Action", max_depth = 3L)
"group_regulation_long"

#' AI-collaboration messages (long format)
#'
#' A long, one-row-per-message log from an AI-assisted collaboration
#' study, with Unix timestamps and an explicit session id. Used to
#' demonstrate long-format loading with Unix time and sessions. Bundled
#' example dataset.
#'
#' @format A \code{data.frame} with 8551 rows and 9 columns:
#'   \describe{
#'     \item{message_id}{integer; the message the code was assigned to.}
#'     \item{project}{character; the project (the actor).}
#'     \item{session_id}{character; the session id.}
#'     \item{timestamp}{integer; Unix time in seconds.}
#'     \item{session_date}{character; the session date (YYYY-MM-DD).}
#'     \item{code}{character; the move (eight codes, the fine state).}
#'     \item{cluster}{character; the coarse state the code belongs to.}
#'     \item{code_order}{integer; order of the code within its message.}
#'     \item{order_in_session}{integer; order within the session.}
#'   }
#' @source Bundled example dataset.
#' @examples
#' data(ai_long)
#' context_tree(ai_long, actor = "project", time = "timestamp",
#'              action = "code", max_depth = 2L)
"ai_long"

#' Student engagement state sequences (wide format)
#'
#' A wide set of student engagement-state sequences, one learner per row
#' and one time-step per column. Used to demonstrate fitting from wide
#' sequence data and from the network objects (\pkg{tna} / \pkg{Nestimate})
#' built on top of it. Bundled example dataset.
#'
#' @format A \code{data.frame} with 1000 rows (learners) and 25 columns
#'   (time-steps). States \code{"Active"}, \code{"Average"},
#'   \code{"Disengaged"}; \code{NA} marks missing / dropped-out positions.
#' @source Bundled example dataset.
#' @examples
#' data(engagement)
#' context_tree(engagement, max_depth = 2L)
"engagement"
