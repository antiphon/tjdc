#' Contextual Mann-Kendall Trend Test
#'
#' Mann-Kendall trend test for 2D rasterStacks, assuming layers
#' represent timepoints.
#'
#' @param x stars-object, with x,y and time dimensions, and at least one attribute.
#' @param ... ignored
#' @param neighbourhood 2:queen, 1:rook, 0:none (classical M-K test) Only 2 and 0 implemented.
#' @param calc_slope Calculate Theil-Sen slope estimate for each series? Default:FALSE
#' @param timevar Name of the time dimension. Default "z"
#' @param attrvar Name of the value attribute. Default "values
#' @details  Assume that the values of each rasterStack location (cell) over time (layers) are a time-series. This function calculates the Mann-Kendall trend test for each series, and returns the relevant statistics as a new rasterStack.
#'
#'  The **Contextual** version (default; see ref.) simply averages the
#'  Mann-Kendall-test's S-statistic at
#'  each location over neighbours, very much like running 'focal(S, w=matrix(1/9,3,3))'
#'  on the non-contextual/cellwise statistic raster.
#'  However, this function also calculates the adjusted variances (and hence the p-values).
#'
#' NOTE assumes equal interval timeseries if 'time' is missing. Matters only for the slope calculation.
#'
#'
#' @return A stars-object with following attributes: 'S' for the Mann-Kendall Statistic (classical or smoothed);
#'  's2' for the variance of the statistic; 'p' for p-value of trend detected;
#' (optional) 'slope' for the Theil-Sen slope estimate.
#'
#' @references
#' Neeti, N. and Eastman J.R. (2011) A Contextual Mann-Kendall Approach for the Assesment of Trend Significance in Image Time Series, \emph{Transactions in GIS}
#' @useDynLib tjdc
#' @import stars dplyr tidyr
#' @export

tj_contextual_mann_kendall <- function(x, ..., neighbourhood = 2, calc_slope = FALSE,
                                       timevar = "z",
                                       attrvar = "values") {

  if(!neighbourhood %in% c(0,1,2)) stop("Only '0', '1' and '2' neighbourhoods implemented.")
  dims <- dim(x)
  if(!all( c("x", "y") %in% names(dims) )) stop("Dimensions 'x' and 'y' not found.")
  #if( !canProcessInMemory(x) ) stop("'canProcessInMemory' return FALSE")
  t0 <- Sys.time()

  # center each series? do it in tibble pipe
  #x <- st_apply(x[attrvar], c("x", "y"), \(v) v - mean(v, na.rm=TRUE), .fname = timevar)

  # Extract values.
  Xt <- as_tibble(x) |>
    select(all_of(c("x", "y", timevar, attrvar))) |>
    transmute(x, y, time = get(!!timevar), value = get(!!attrvar) ) |>
    # center each series
    group_by(x, y) |>
    mutate(value = value - mean(value, na.rm=TRUE)) |>
    ungroup() |>
    pivot_wider(values_from = value, names_from = time)

  X <- Xt |>
    # Make a matrix, cells x time
    select(-x, -y) |>
    as.matrix()

  # check time vector
  time_values <- st_get_dimension_values(x, timevar)
  if(length(time_values) != ncol(X) || !is.numeric(time_values)) stop(paste0("'time' should be a numeric increasing vector of length ", ncol(X)))

  # compute
  res <- c_contextual_mann_kendall( x = X,
                                    nrow = dims['y'],
                                    time = time_values,
                                    neighbourhood = neighbourhood,
                                    calc_slope = calc_slope)

  # gather results
  V <- res$S_and_s2
  Sm <- V[,1]
  s2 <- V[,2]
  D <- (Sm > 0) * 1
  Z <- (Sm+(1-2*D))/sqrt(s2)
  p <- 2 * (1 - abs(pnorm(abs(Z))))

  out  <- slice(x, !!timevar, 1) |>
    mutate(values = NULL,
               S  = Sm,
               s2 = s2,
               p  = p)

  if(calc_slope)  out$slope <- res$slope

  attr(out, "timing") <-Sys.time() - t0
  out
}

