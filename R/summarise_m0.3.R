#' Summarise fit of M0.3 as a tibble
#'
#' @param x object fitted by this package's tj_fit_m0.3
#' @param burnin if mcmc traces included (`keep_hist` was not false), recalculate summaries using the iteration vector (burnin*nsteps):nsteps
#' @param thin_steps take only every `thin_steps` from the iteration vector
#'
#' @import stars dplyr
#' @export
tj_summarise_m0.3 <- function(x, ..., cells, burnin, thin_steps) {

  if(missing(cells)) cells <- x$cell_info$cell # do alla
  # in case we re-summarise:
  h <- x$keep_hist
  #
  ss <- function(x)       c(mean = mean(x), var = var(x), q=quantile(x, c(0.025, .1, 0.5, .9, 0.975)))
  #
  if( (length(h) == 1 && h) &
     !missing(burnin) & !missing(thin_steps)) {
    # we got all, re summarise
    o                    <- seq(burnin * x$niter, x$niter, by = thin_steps )
    K                    <- length(x$timesteps)
    post_k               <- apply(x$hist_k[o,cells,drop=FALSE], 2, \(x) { f <- tabulate(x, K); f/sum(f)})
    post_z               <- colMeans(x$hist_z[o,cells,drop=FALSE])
    post_theta_m         <- apply(x$hist_theta[o,cells,, drop=FALSE], c(2, 3), mean)
    post_theta_sd        <- apply(x$hist_theta[o,cells,, drop=FALSE], c(2, 3), sd )
    post_sigma2          <- ss(x$hist_sigma2[o])
  }
  else {
    post_k               <- x$k[,cells, drop=FALSE]# |> t()
    post_z               <- x$z[cells]
    post_theta_m         <- x$theta_m[,cells, drop=FALSE]  |> t()
    post_theta_sd      <- x$theta_Suppertri[c(1,4,6), cells, drop=FALSE] |> t()
    post_sigma2          <- x$sigma2
  }
  #
  #browser()
  # Posterior map
  pred <- x$cell_info |>
    filter(cell %in% cells) |>
    mutate(
      theta_mean = post_theta_m |> as_tibble(.name_repair = "unique") |> setNames(c("a","b","d")),
      theta_sd  = post_theta_sd |> as_tibble(.name_repair = "unique") |> setNames(c("a","b","d")),
      pred_jump_prob = post_z,
      pred_jump_mode = apply(post_k, 2, which.max),
      pred_jump_prob_at_mode = apply(post_k, 2, max)
           )
  attr(pred, "sigma2") <- post_sigma2
  attr(pred, "timesteps") <- x$timesteps
  pred
}



#' Prediction the series per cell
#'
#' @param x object fitted by this package's tj_fit_m0.3
#' @param cells which cells to predict? Default: all
#' @param s output of summarise_m0.3
#'
#' @import dplyr
#' @export
tj_predict_m0.3 <- function(x, cells, ..., s) {
  #
  if(missing(s)) s <- summarise_m0.3(x, cells = cells, ...)
  # Modes
  if(missing(cells)) cells <- s$cell # do all
  if(!all(cells %in% s$cell)) stop("cells and summary do not match properly.")

  timesteps <- attr(s, "timesteps")

  wra <- function(k) {
    sk  <- s |> filter(cell == k)
    j   <- sk$pred_jump_mode
    kv <- rep(1, length(timesteps))
    kv[1:j] <- 0
    be  <- sk$theta_mean |> as.matrix()
    y <- c(cbind(1, timesteps, kv) %*% t(be))
    tibble(time = timesteps, value=y)
  }

  s |>
    filter(cell %in% cells) |>
    select(cell) |>
    rowwise() |>
    summarise(cell = cell, wra(cell))
}


#' Get the history
#'
#' Compile the MCMC trace for one cell
#'
#' @param x fit_m0.3 with `keep_hist=TRUE`
#' @param cell which cell
#' @details
#'
#' Combine with posterior-package.
#'
#' @import dplyr
#' @export

tj_trace_m0.3 <- function(x, cell, ...) {
  h <- x$keep_hist
  if(is.null(x$hist_k)) stop("Model not estimated with keep_hist = TRUE")
  th <- x$hist_theta[,cell,] |>
    as_tibble() |>
    setNames(c("a", "b","d"))
  tj <- tibble(k = x$hist_k[,cell], z = x$hist_z[,cell])
  ts <- tibble(sigma2 = c(x$hist_sigma2) )

  bind_cols(th, tj, ts)
}

