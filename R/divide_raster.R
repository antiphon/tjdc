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
  warning("tj_divide_box has been deprecated. Use tj_divide_raster instead.")
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
#' Divide the domain of raster into regular tiling, with optional overlap.
#'
#' @param x stars raster
#' @param n two-vector giving number of tiles per in x and y, in that order.
#' @param buffer two-vector of buffering the sub-rasters in each direction, in raster pixels.
#' @param open.upper Remove overlap? (TRUE)
#'
#' @details
#' Will operate on the row-col indices of the stars-object, with row top to bottom along y-axis, and
#' col left to right on x axis.
#'
#' The tiles will go top-bottom, left to right order in space. The buffer is on each tile, so
#' for buffer (t,t) will lead to neighbouring tiles having 2t cell-width row and/or col in common. To get
#' 2t+1-cell overlaps set `open.upper=FALSE`, but note then the minimum is 1-cell overlap.
#'
#'
#'
#' @return A list with st-polygons, cells, and row-col-ranges, and row-col-ranges without a buffer.
#'
#'
#' The row-col ranges are closed by default,
#' so even with `buffer = c(0,0)` they will overlap at the ends. This can be avoided by setting `open.upper=TRUE` whence the
#' range will not include the highest row/col, apart from the tile at the bottom row /right columns.
#'
#' Cell-sets are computed from the row-col ranges.
#'
#' The polygons are given by the st_bbox of `x` subindexed with the [expand.grid]'d vectors between the row-col ranges.
#'
#'
#' @import sf stars dplyr
#' @export
tj_divide_raster <- function(x,
                             n = c(2,2),
                             buffer = c(0,0),
                             open.upper = TRUE) {

  if(length(buffer) != 2) buffer <- rep(buffer, 2)[1:2]
  rc   <- dim(x)[c('y', 'x')] |> unname()
  nx <- n[1]
  ny <- n[2]


  dx <- ceiling(rc[2]/nx)
  dy <- ceiling(rc[1]/ny)
  rclist <- vector("list", length = nx * ny)
  rclist0 <- vector("list", length = nx * ny) # no buffer
  wlist  <- vector("list", length = nx * ny)
  cellist <- vector("list", length = nx * ny)
  n <- 1
  for (j in seq_len(ny) - 1) {
    for (i in seq_len(nx) - 1) {
      x0 <- 1 +       i * dx
      x1 <- min(rc[2], 1 + (i + 1) * dx )
      y0 <- 1 +       j * dy
      y1 <- min(rc[1], 1 + (j + 1) * dy )
      if(open.upper) {
        if(j < ny-1){ # not last column
          y1 <- y1 - 1
        }
        if(i < nx-1) { # not last row
          x1 <- x1 - 1
        }
      }
      # Buffer
      xr <- c(max(1,     x0 - buffer[1]), min(rc[2], x1 + buffer[1]))
      yr <- c(max(1,     y0 - buffer[2]), min(rc[1], y1 + buffer[2]))

      # stars bounding box
      w  <- st_bbox( x[,xr[1]:xr[2], yr[1]:yr[2], 1, drop=TRUE] ) |> st_as_sfc()
      # cells
      xxx <- expand.grid(xr[1]:xr[2], yr[1]:yr[2])
      cl <- xxx[,1] + (xxx[,2]-1) * rc[2]
      # store
      rclist0[[n]] <- c(rmin = y0, rmax = y1, cmin = x0, cmax = x1 ) # original picture rc without buffer
      rclist[[n]]  <- c(rmin = yr[1], rmax = yr[2], cmin = xr[1], cmax = xr[2])
      wlist[[n]]   <- w
      cellist[[n]] <- cl
      n <- n + 1L
    }
  }
  wlist <- lapply(wlist, st_as_sf) |> bind_rows() |> mutate(id = 1:n())
  rclist <- rclist |> bind_rows() |> mutate(id = 1:n())
  rclist0 <- rclist0 |> bind_rows() |> mutate(id = 1:n())

  out <- list(cell = cellist, poly = wlist, rc = rclist,
       rc_no_buffer = rclist0, buffer = buffer, n = c(nx, ny), open.upper=open.upper)
  class(out) <- c( is(out), "tj_tiling" )
  out
}


