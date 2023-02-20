#' Simultaneous linear regression with correlated changepoint estimation
#'
#' For each pixel in a spatiotemporal raster x time datacube,estimate a linear regression model
#' with 0/1 jumps such that the jump event is correlated within neighbouring pixels.
#'
#' @param x input stars object
#' @param niter number of gibbs sampler interations
#' @param prior_k prior jump probabilities. Either a vector (will be normalised) or a single value for the probability of jump somewhere.
#' @param dat prepared data frame (for development, please dont use).
#' @param prior_theta prior (m=vec, S=matrx) for intercept and trend parameter, same for all cells
#' @param prior_delta prior(m=num, s2=num>0, a=num, b=num>a) truncated normal prior for eac jump, trucate to [a,b]
#' @param keep_hist If TRUE, keep the mcmc trajectories of all variables. If vector of cell indices, keep history for only those (to save space). Give at least 2 cell-indices or logical, otherwise unexpected behaviour. (default: FALSE).
#' @param neighbourhood_type How is neighbourhood defined (default: "square")
#' @param neighbourhood_range How large neighbourhood (default: 1)
#' @param ... ignored
#'
#' @details We assume the datacube holds for each pixel (x-y) a series of values (`attrvar`) in time (`timevar`).
#'
#' The model estimate for each series the best linear regression line such that there might be
#' a jump of size delta after one of the time points 1...T-1. The special estimate "jump after T"
#' means no jump is predicted.
#'
#'
#' @import looptimer dplyr
#' @export


tj_fit_m0.6 <- function(x,
                        dat = NULL, # make sure you know what you are doing here!
                        timevar = "z",
                        attrvar = "values",
                        niter = 1000,
                        gamma = 0,
                        prior_theta  = list(m = c(a=0, b=0),
                                            S = diag(1e5*c(1, 1))),
                        prior_delta = list(m=0, s2 = 1e5, a=-Inf, b=Inf), # truncate delta
                        prior_sigma2 = c(shape = 2, rate = 1), # inv-gamma
                        prior_k = 0.5,
                        verbose = FALSE,
                        ctrl = list(burnin = 0.5, thin_steps = 1),
                        keep_hist = TRUE,
                        cells_to_ignore = NULL,
                        ignore_cells_with_na  =  TRUE,
                        neighbourhood_type = "square",
                        neighbourhood_range = 1,
                        ... # ignored
) {

  t0 <- Sys.time()
  cat2 <- if(verbose) function(...) message(appendLF=FALSE, ...) else function(...) NULL
  if(verbose) timer1 <- looptimer(n = niter, endline = "     \r", prefix = "[m0.6 mc]", printevery = ceiling(niter/200)  )
  #
  # Prepare input:
  # stars to data frame
  if(is.null(dat)) dat <- tj_stars_to_data(x, timevar, attrvar)
  # check that we have expected names in the dat
  gdims <- attr(dat, "grid")
  if( !all(c("cell", "time", "value") %in% names(dat))  | is.null(gdims))
    stop("'dat' ill formatted. need at least cell, time, value-columns, and `gdims` attribute. See `tj_stars_to_data`. ")
  #
  #
  #
  #
  # assume complete rectangular grid.
  nr <- gdims['nrow']
  nc <- gdims['ncol']
  cells <- unique(dat$cell)
  if(length(keep_hist) > 1) if(!all(keep_hist %in% cells)) stop("Input error: keep_hist")
  #
  n <- nr * nc # the complete raster
  nlist <-  if(gamma != 0)  lapply(1:n, tj_cell_neighbours,
                                   nr = nr, nc = nc,
                                   type = neighbourhood_type,
                                   range = neighbourhood_range) else vector("list", n)

  # We might have weights.
  neighbours_have_weights <- FALSE
  if(!is.null(attr(nlist[[1]], "weight"))){
    nlist_weights <- lapply(nlist, attr, "weight")
    neighbours_have_weights <- TRUE
  }

  #browser()

  # Keep this info
  timesteps <- sort( unique(dat$time) )
  cell_info <- dat[dat$time == timesteps[1], c("cell", "c_x", "c_y")]


  ####################################################################################
  # Check for missing values and ignore them?
  if(ignore_cells_with_na) {
    # check which have na and flag them
    bad <- dat |>
      dplyr::arrange(cell, time) |>
      dplyr::group_by(cell) |>
      dplyr::summarise(n_na = sum(is.na(value))) |>
      dplyr::filter(n_na > 0)
    cells_to_ignore <- union(cells_to_ignore, bad$cell)
  }
  # These cells will be in play.
  cells_to_consider <- setdiff(cells, cells_to_ignore)
  # remove the unwanted from neighbourhoods
  if(!is.null(cells_to_ignore)) {
    #
    #nlist <- lapply(nlist, setdiff, cells_to_ignore)
    nlist_keep <- lapply(nlist, \(n) !(n%in% cells_to_ignore))
    nlist <- lapply(seq_along(nlist), \(i) nlist[[i]][nlist_keep[[i]]] )
    if(neighbours_have_weights)
      nlist_weights <- lapply(seq_along(nlist_weights), \(i) nlist_weights[[i]][nlist_keep[[i]]] )
  }
  #
  # Keep only relevant, and arrange by cell and then time.
  dat <- dat |>
    filter(cell %in% cells_to_consider) |>
    arrange(cell, time)
  ####################################################################################
  # Check for hermits. These will need to be known to avoid subsetting G with numeric(0)
  nnsizes <- sapply(nlist, length)
  #
  #
  # format value-data into a matrix:
  K <- length(timesteps)  # assuming all same length
  Y <- matrix(ncol = n, nrow = K)
  #browser()
  # This is very crusial step here! Make sure dat is in the right order.
  Y[, cells_to_consider] <- dat$value
  N <- nrow(dat)
  ##########################################
  # priors
  if(length(prior_k) == 1) prior_k <- c( rep(1, K-1), (K-1) * (1-prior_k)/prior_k )
  prior_k <- rep(prior_k, K)[1:K]
  prior_k <- prior_k / sum(prior_k)

  ##########################################
  # Precompute for linear stuff
  X <- cbind(1, timesteps)
  X1_list    <- lapply(1:K, \(k) cbind(X, rep(0:1, c(k,K-k))) )
  pX         <- ncol(X)+1
  tX1_list   <- lapply(X1_list, t)
  tX1X1_list <- lapply(X1_list, \(x1) t(x1)%*%x1)
  tX <- t(X)
  tXX <- tX%*%X
  #
  # delta, jumpsize, is the last parameter in theta
  theta_priorS  <- prior_theta$S
  theta_priorSi <- prior_theta$S |> solve()
  theta_priorm  <- prior_theta$m
  theta_priorQm <- theta_priorSi %*% theta_priorm
  #
  # @todo Externalize this.
  rtnorm <- function( mu, sigma, a = prior_delta$a, b = prior_delta$b, n = 1) {
    #truncnorm::rtruncnorm(n = n, a = a, b = b, mean = mu, sd = sigma)
    qnorm( pnorm( (a-mu)/sigma ) + runif(n) * (pnorm((b-mu)/sigma)-pnorm((a-mu)/sigma)) ) *sigma + mu
  }
  theta  <- mvtnorm::rmvnorm(n, prior_theta$m, prior_theta$S)
  delta  <- rtnorm(mu = prior_delta$m,
                   sigma = sqrt(prior_delta$s2), n = n)
  theta  <- cbind(theta, delta)
  #
  # The interaction matrix
  G <- diag(gamma, K)
  # initialise chain: Random jump times
  kvec   <- sample(K, n, replace=TRUE, prob = prior_k + 0.0001)
  #
  # Strategy: Update delta, and others separately all but delta.

  sigma2  <- 1/rgamma(1, prior_sigma2[1], prior_sigma2[2] )

  # fixed over the pixel-wise update
  #theta_Snew_list <- lapply(1:K, \(ki)
  #                          solve(tX1X1_list[[ki]]/sigma2 + theta_priorSi))
  #chol_theta_Snew_list <- lapply(theta_Snew_list, chol)
  theta_Snew <- solve(tXX / sigma2 + theta_priorSi)
  chol_theta_Snew <- chol(theta_Snew)

  # follow something
  history_k       <- matrix(NA, nrow = niter, ncol = n)
  history_z       <- matrix(NA, nrow = niter, ncol = n) # probability of jump anywhere
  history_theta   <- array( dim = c(niter, n, ncol(theta)))
  history_sigma2  <- matrix(NA, nrow = niter, ncol = 1)
  #
  # verbosity freq
  B <- ceiling( niter / 100)

  # Then we iterate
  for(it in 1:niter) {
    resid2_sum <- 0 # for sigma update

    for(i in sample(cells_to_consider) ) {
      #
      yi   <- Y[,i]
      ei   <- c(yi - X %*% theta[i, -3]) # assuming 3d
      ei2  <- ei^2
      eid  <- ei - theta[i, 3]
      eid2 <- eid^2

      # update k: square sums of before-after jump
      ss_upto  <- cumsum(ei2[-K])
      ss_after <- sum(eid2[-1]) - cumsum(c(0, eid2[-c(1,K)]) )
      SS     <- c(ss_upto + ss_after, sum(ei2)) # and without jump
      # Then the Potts smoothing
      if(nnsizes[[i]]>0){
        pott <- if(neighbours_have_weights) rowSums(nlist_weights[[i]] * G[, kvec[ nlist[[i]] ], drop=FALSE  ])
                else rowSums(G[, kvec[ nlist[[i]] ], drop=FALSE  ]) # n.of kth type neighbours times weight.
      }
      else pott <- 0
      # Then the probability to be softmaxed
      lp   <- -1/(2*sigma2) * SS  + pott + log( prior_k )
      pv   <- exp( lp-max(lp) ) # shift to avoid singularities, softmax is invariant
      pv   <- pv / sum(pv)
      #if(any(is.na(pv))) browser()
      knew <- sample(K, 1, prob = pv )
      kvec[i] <- knew
      # track jump-anywhere probability
      history_z[it,i] <- sum(pv[-K])
      #
      # Update linear parameters
      # Sample full covariance Normal using chol.
      # The residuals, given delta
      yiprime            <- yi
      yiprime[-(1:knew)] <- yiprime[-(1:knew)] - theta[i,3]
      mnew            <- theta_Snew %*% ( tX%*%yiprime /sigma2  + theta_priorQm)
      thetanew        <- (rnorm(nrow(theta_Snew)) %*% chol_theta_Snew)[1,] + mnew
      #
      ei              <- c(yi - X %*% thetanew) # new residuals without jump
      #browser()
      # then delta, given theta, truncated normal
      ds2new          <- 1/((K-knew)/sigma2 + 1/prior_delta$s2 )
      dmnew           <- ds2new * (sum(ei[-(1:knew)])/sigma2 +
                                     prior_delta$m/prior_delta$s2)
      deltanew        <-  rtnorm(dmnew, sqrt(ds2new) )
      #if(knew < 11) browser()
      #
      thetanew        <- c(thetanew, deltanew)
      theta[i,]       <- thetanew
      #
      # new residual square sum, for updating sigma
      ei[-(1:knew)]  <- ei[-(1:knew)] - theta[i, 3]
      resid2_sum     <- resid2_sum + sum(ei^2)
      #
      #if(it == 300 & i == 333) browser()
    } # eo pixel wise iterations
    # Update sigma2
    a_new  <- prior_sigma2[1] + (n*K)/2
    b_new  <- prior_sigma2[2] + resid2_sum/2
    sigma2 <- 1/rgamma(1, a_new, b_new)
    #theta_Snew_list <- lapply(1:K, \(ki)
    #                          solve(tX1X1_list[[ki]]/sigma2 + theta_priorSi))
    #chol_theta_Snew_list <- lapply(theta_Snew_list, chol)
    theta_Snew  <- solve(tXX/sigma2 + theta_priorSi)
    chol_theta_Snew <- chol(theta_Snew)

    #
    #if(verbose && (it %% B == 0)) cat2(sprintf("m0.6mc: %3.0f%%   \r", 100*it/niter))
    if(verbose) print( timer1 <- looptimer(timer1) )
    history_k[it,]         <- kvec
    history_theta[it,,]    <- theta
    history_sigma2[it,]    <- sigma2
  }
  cat2("\n")
  #browser()
  #
  # Compute best estimates
  o                    <- seq(ctrl$burnin * niter , niter, by = ctrl$thin_steps )
  post_k               <- apply(history_k[o,,drop=FALSE], 2, \(x) { f <- tabulate(x, K); f/sum(f)})
  post_z               <- colMeans(history_z[o,,drop=FALSE])
  post_theta_m         <- apply(history_theta[o,,], 2, colMeans)
  post_theta_Suppertri <- apply(history_theta[o,,], 2, \(x){ S <- cov(x); c(S[upper.tri(S, diag = TRUE)])} )
  ss <- function(x)       c(mean = mean(x), var = var(x), q=quantile(x, c(0.025, .1, 0.5, .9, 0.975)))
  post_sigma2          <- ss(history_sigma2[o])
  #
  post_k[,cells_to_ignore ] <- NA
  post_z[cells_to_ignore ]  <- NA
  # post_thetam_m[cells_to_ignore ] <- NA

  # remove trace objects if not needed
  if(length(keep_hist) > 1) {
    history_k <- history_k[, keep_hist, drop = FALSE]
    history_z <- history_z[, keep_hist, drop = FALSE]
    history_theta <- history_theta[, keep_hist, , drop = FALSE]
  }
  else if(!keep_hist) { # drop all
    history_k <- NULL
    history_z <- NULL
    history_theta <- NULL
    #history_sigma2 <- NULL # always keep sigma2 history
  }
  cell_info <- cell_info |>
    mutate( ignored = cell %in% cells_to_ignore )

  out <- list(took = Sys.time() - t0,
              k = post_k,
              z = post_z,
              theta_m  = post_theta_m,
              theta_Suppertri = post_theta_Suppertri,
              sigma2 = post_sigma2,
              hist_k = history_k,
              hist_z = history_z,
              hist_theta = history_theta,
              hist_sigma2 = history_sigma2,
              keep_hist = keep_hist,
              niter = niter,
              cell_info  = cell_info,
              timesteps = timesteps,
              ctrl = ctrl,
              gdims = gdims,
              cells_ignored = cells_to_ignore,
              neighbourhood = list(type = neighbourhood_type,
                                   range = neighbourhood_range))

  class(out) <- c(is(out), "tj_fit_mc")
  out
}


######################################################################################
#' Fit m0.6 alternative parametrisation
#'
#' Wrapper for m0.6 to be called with data and a single list-item of run-parameters
#'
#' @seealso [tj_cfg_m0.6] for creating required list.
#'
#' @export
tj_fit_m0.6s <- function(..., cfg = tj_cfg_m0.6()) {
  do.call( tj_fit_m0.6, c(list(...), cfg) )
}


#' Compile paramaters for m0.6
#'
#' @export
tj_cfg_m0.6 <- function(...) {
  cfg0 <- list(
    niter = 1000,
    gamma = 0,
    prior_theta  = list(m = c(a=0, b=0), S = diag(1e5*c(1, 1))), # 2d-Normal
    prior_delta = list(m = 0, s2 = 1e5, a = -Inf, b = Inf), # truncated normal
    prior_sigma2 = c(shape = 2, rate = 1), # inv-gamma
    prior_k = 0.5,
    verbose = FALSE,
    ctrl = list(burnin = 0.5, thin_steps = 1),
    keep_hist = TRUE,
    cells_to_ignore = NULL,
    ignore_cells_with_na  =  TRUE,
    neighbourhood_type = "square",
    neighbourhood_range = 1
  )
  #browser()
  new <- list(...)
  for(e in names(new)) {
    cfg0[[e]] <- new[[e]]
  }
  cfg0
}
