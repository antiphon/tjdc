% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gaussian_regression_with_mcmc-v3.R
\name{tj_fit_m0.3}
\alias{tj_fit_m0.3}
\title{Simultaneous linear regression with correlated changepoint estimation}
\usage{
tj_fit_m0.3(
  x,
  dat = NULL,
  timevar = "z",
  attrvar = "values",
  niter = 1000,
  gamma = 0,
  prior_theta = list(m = c(a = 0, b = 0, d = 0), S = diag(1e+05 * c(1, 1, 1))),
  prior_sigma2 = c(shape = 2, rate = 1),
  prior_k = 0.5,
  verbose = FALSE,
  ctrl = list(burnin = 0.5, thin_steps = 1),
  keep_hist = TRUE,
  cells_to_ignore = NULL,
  ignore_cells_with_na = TRUE,
  truncate_jump_at_mean = 0,
  ...
)
}
\arguments{
\item{x}{input stars object}

\item{dat}{prepared data frame (for development, please dont use).}

\item{niter}{number of gibbs sampler interations}

\item{prior_k}{prior jump probabilities. Either a vector (will be normalised) or a single value for the probability of jump somewhere.}

\item{keep_hist}{If TRUE, keep the mcmc trajectories of all variables. If vector of cell indices, keep history for only those (to save space). Give at least 2 cell-indices or logical, otherwise unexpected behaviour. (default: FALSE).}

\item{truncate_jump_at_mean}{use a truncated normal prior for the jump? Will truncate at prior mean ad-hoc, sign says which side to keep.}

\item{...}{ignored}
}
\description{
For each pixel in a spatiotemporal raster x time datacube, estimate a linear regression model
with 0/1 jumps such that the jump event is correlated within neighbouring pixels.
}
\details{
We assume the datacube holds for each pixel (x-y) a series of values (`attrvar`) in time (`timevar`).

The model estimate for each series the best linear regression line such that there might be
a jump of size delta after one of the time points 1...T-1. The special estimate "jump after T"
means no jump is predicted.
}
