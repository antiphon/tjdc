#' Divide bounding box into subsets, return polygons
#'
#' Divide both sides by a factor 2^q
#'
#' @import sf
#' @export
tj_divide_bbox <- function(bb,
                           power_of_2 = 2,
                           q = power_of_2,
                           expand = 0) {

  poly3 <- bb |> st_as_sfc()
  # Buffer?
  ex <- expand
  bb3 <- bb
  # shift
  poly30 <- poly3 - bb3[c(1,2)]
  factr  <- 1/(2^q/(1+2*ex))
  polytmpl <- poly30 * factr - c(bb3[3]-bb3[1], bb3[4]-bb3[2])*ex * factr # corner in
  # shifts
  n <- 2^q+1
  crns <- expand.grid(seq(bb3[1], bb3[3], l = n)[-n], seq(bb3[2], bb3[4], l = n)[-n])
  polys <- apply(crns, 1, \(s) polytmpl + s) |> do.call(what = c) |> st_as_sf( crs = st_crs(poly3)) |>
    mutate(id = sprintf("%04i", 1:n()) )
  polys
}


#' Tile stars raster
#'
#' Divide both sides by a factor 2^q
#' @param x stars raster
#' @param n two-vector giving number of splits per x,y
#' @param buffer two-vector of buffering the sub-rasters in each direction
#' @param cells return cell numbers of x (TRUE) or st-polygons?
#'
#' @return A set of st-polygons, cells, and row-col-ranges
#' @import sf stars
#' @export
tj_divide_raster <- function(x,
                             n = c(2,2),
                             buffer = c(0,0)) {
  rc   <- dim(x)[c('y', 'x')]
  nx <- n[1]
  ny <- n[2]
  dx <- ceiling(rc[2]/nx)
  dy <- ceiling(rc[1]/ny)
  rclist <- vector("list", length = nx * ny)
  wlist  <- vector("list", length = nx * ny)
  cellist <- vector("list", length = nx * ny)
  n <- 1
  for (i in seq_len(nx) - 1) {
    for (j in seq_len(ny) - 1) {
      x0 <- 1 +       i * dx
      x1 <- min(rc[2], 1 + (i + 1) * dx )
      y0 <- 1 +       j * dy
      y1 <- min(rc[1], 1 + (j + 1) * dy )
      rclist[[n]] <- c(y0, y1, x0, x1) # original picture rc without buffer
      xr <- c(max(1,     x0 - buffer[1]):min(rc[2], x1 + buffer[1]))
      yr <- c(max(1,     y0 - buffer[2]):min(rc[1], y1 + buffer[2]))
      wlist[[n]]  <- st_bbox( x[,xr, yr, 1, drop=TRUE] ) |> st_as_sfc()
      xxx <- expand.grid(xr, yr)
      cellist[[n]] <- xxx[,1] + (xxx[,2]-1) * nx
      n <- n + 1L
    }
  }

  list(cell = cellist, poly = wlist, rc = rclist)
}


