#' Student engagement trajectories
#'
#' @description
#' \code{trajectories} is a coded dataset of the engagement states of 138
#' students at 15 successive time points, one row per student and one column per
#' time point. Each state is \code{"Active"}, \code{"Average"} or
#' \code{"Disengaged"}. Missing values occur only as a run at the end of a row
#' and mark that the student's trajectory ended (drop-out); two rows are missing
#' throughout and are ignored when a tree is fitted. The matrix has the wide
#' layout that \code{\link{context_tree}()} fits directly, and is small enough
#' for fast examples.
#'
#' @format A character matrix with 138 rows (students) and 15 columns
#'   (time points, named \code{"1"} to \code{"15"}), holding
#'   \code{"Active"}, \code{"Average"}, \code{"Disengaged"} or \code{NA}.
#'   118 rows are complete; 1870 states are observed in all.
#'
#' @source Bundled example dataset. The same matrix is distributed with the
#'   \pkg{Nestimate} package as \code{trajectories}.
#'
#' @seealso \code{\link{engagement}} for a larger engagement dataset;
#'   \code{\link{context_tree}()}.
#'
#' @examples
#' tree <- context_tree(trajectories, max_depth = 3)
#' tree
#' state_frequencies(tree)
"trajectories"

#' Regulation of collaborative learning (long event log)
#'
#' @description
#' \code{group_regulation_long} is a coded dataset of the regulation actions of
#' 2,000 students working in 200 small groups across three courses, one row per
#' action with a timestamp. Each action is one of nine regulation moves:
#' \code{adapt}, \code{cohesion}, \code{consensus}, \code{coregulate},
#' \code{discuss}, \code{emotion}, \code{monitor}, \code{plan} and
#' \code{synthesis}. Each student is classified as a high or low achiever (1,000
#' each), which makes the data suited to comparing the regulation pathways of
#' two groups.
#'
#' The data are in long format. Name the columns in
#' \code{\link{context_tree}()} (\code{actor = "Actor"},
#' \code{time = "Time"}, \code{action = "Action"}) or reshape them with
#' \code{\link{prepare_input}()}; with the default 15-minute session gap,
#' each student's actions form one sequence of 2 to 26 actions.
#'
#' @format A \code{data.frame} with 27,533 rows (actions) and 6 columns:
#'   \describe{
#'     \item{Actor}{Integer student identifier, 1 to 2,000.}
#'     \item{Achiever}{Achievement level, \code{"High"} or \code{"Low"},
#'       constant within a student.}
#'     \item{Group}{Numeric collaboration group, 1 to 200; each student
#'       belongs to one group.}
#'     \item{Course}{Course, \code{"A"}, \code{"B"} or \code{"C"}.}
#'     \item{Time}{\code{POSIXct} timestamp of the action.}
#'     \item{Action}{The regulation move, the state of the sequence.}
#'   }
#'
#' @source Identical to \code{group_regulation_long} in the \pkg{tna}
#'   package (there a tibble), where it is described as generated
#'   synthetically from that package's \code{group_regulation} data.
#'
#' @seealso \code{\link{compare_groups}()}, \code{\link{prepare_input}()}.
#'
#' @examples
#' groups <- context_tree(group_regulation_long, actor = "Actor",
#'                        time = "Time", action = "Action",
#'                        group = "Achiever", max_depth = 2)
#' summary(groups)
"group_regulation_long"

#' AI turns in human-AI programming sessions (long event log)
#'
#' @description
#' \code{ai_long} is a coded dataset of the turns of an AI coding assistant in
#' 428 human-AI pair-programming sessions across 34 projects, held between 5
#' February and 8 March 2026, one row per code. A turn (message) can carry
#' several codes, in the order given by \code{code_order}. Each code is one of
#' eight assistant actions (\code{Ask}, \code{Delegate}, \code{Execute},
#' \code{Explain}, \code{Investigate}, \code{Plan}, \code{Repair},
#' \code{Report}), grouped into three clusters (\code{Action},
#' \code{Communication}, \code{Repair}).
#'
#' Timestamps have day resolution (midnight, UTC+2, of the session date), so
#' the order of events within a day comes from \code{order_in_session},
#' which counts the human turns of the session as well and therefore has
#' gaps. Use \code{actor = "session_id"} with
#' \code{order = "order_in_session"} for one sequence per session, or
#' \code{time = "timestamp"} to split sequences at changes of day.
#'
#' @format A \code{data.frame} with 8,551 rows (codes) and 9 columns:
#'   \describe{
#'     \item{message_id}{Integer identifier of the turn; the codes of one
#'       turn share it.}
#'     \item{project}{Project identifier, \code{"Project_1"} to
#'       \code{"Project_34"}.}
#'     \item{session_id}{Session identifier, a 12-character hash.}
#'     \item{timestamp}{Integer Unix time in seconds (day resolution).}
#'     \item{session_date}{Date of the row, \code{"YYYY-MM-DD"}.}
#'     \item{code}{The assistant action, the state of the sequence.}
#'     \item{cluster}{The cluster of the code.}
#'     \item{code_order}{Order of the code within its turn, from 1.}
#'     \item{order_in_session}{Position of the row in its session, counting
#'       human and AI rows together.}
#'   }
#'
#' @source Bundled example dataset. The same data are distributed with the
#'   \pkg{Nestimate} package as \code{ai_long}.
#'
#' @seealso \code{\link{prepare_input}()}, \code{\link{context_tree}()}.
#'
#' @examples
#' tree <- context_tree(ai_long, actor = "session_id",
#'                      order = "order_in_session", action = "code",
#'                      max_depth = 2)
#' state_frequencies(tree)
#' common_pathways(tree, top = 5, depth = 2)
"ai_long"

#' Student engagement sequences (wide format)
#'
#' @description
#' \code{engagement} is a synthetic dataset of the engagement states of 1,000
#' students over a study program, generated to resemble the engagement
#' trajectories analyzed by Saqr, \enc{López-Pernas}{Lopez-Pernas}, Helske and
#' Hrastinski (2023). It has one row per student and one column per time point
#' (23 to 25 observed states per student). Each state is \code{"Active"},
#' \code{"Average"} or \code{"Disengaged"}; missing values occur only at the end
#' of a row.
#'
#' @format A \code{data.frame} with 1,000 rows (students) and 25 character
#'   columns \code{T1} to \code{T25}, holding \code{"Active"},
#'   \code{"Average"}, \code{"Disengaged"} or \code{NA}.
#'
#' @source The \code{engagement} sequence object of the \pkg{tna} package,
#'   converted to character with its void and missing codes set to
#'   \code{NA} (identical values).
#'
#' @references
#' Saqr, M., \enc{López-Pernas}{Lopez-Pernas}, S., Helske, S., &
#' Hrastinski, S. (2023). The longitudinal association between engagement
#' and achievement varies by time, students' profiles, and achievement
#' state: A full program study. \emph{Computers & Education}, 199, 104787.
#' \doi{10.1016/j.compedu.2023.104787}
#'
#' @seealso \code{\link{trajectories}} for a smaller engagement dataset;
#'   \code{\link{context_tree}()}.
#'
#' @examples
#' tree <- context_tree(engagement, max_depth = 3, min_count = 20)
#' tree
#' divergent_pathways(tree, top = 5)
"engagement"
