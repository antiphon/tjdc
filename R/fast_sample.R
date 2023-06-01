#' Sample Element Fast
#'
#'
tj_fast_sample <- function(prob) {
  sum(runif(1) > cumsum(prob)) + 1
}
