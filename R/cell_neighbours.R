#' Raster Cell Neighbourhoods
#'
#' @param i focal cell
#' @param nr number of rows in the xy-grid
#' @param nc number of columns in the xy-grid
#' @param ... ignored
#' @param type "square" (default), "circle", "sqexp"
#' @param range Range of neighbourhood, default is 1 for 1-step-neighbourhoods
#'
#' @details
#' The default is `type='square'` and `range=1`, which is 8 adjacent pixels. Increasing `range=2` includes the neighbours' neighbours (24 pixels).
#'
#' To return neighbourhood with corners removed, use `type='circle'` which drops the pixels from square-neighbourhood with distance larger than `range`.
#'
#' To weight the further away neighbours, use `type='sqexp'`. It is similar to circle, but with weight decreasing exponentially in squared distance,
#' scaled to be 1 for distance 1.
#'
#' Non-toroidal neighborhood, so at edges there will be less neighbours.
#'
#' Will assume cell-numbering
#' cell = col + (row-1) * ncol.
#' @returns
#' A vector of indices. If `type='sqexp'`, will return the weights in the attribute 'weight'.
#'
#' @export
tj_cell_neighbours <- function(i, nr, nc, ..., type = "square", range = 1) {
  #
  # the radius in pixels
  sqexp <- (type == "sqexp")
  if(sqexp) range <- range + 1
  M <- if(type=="square") max(1, floor(range)) else ceiling(range)
  #
  # template for row-column indices
  i0  <- ceiling((M*2+1)^2/2)
  rcn <- expand.grid(-M:M, -M:M, KEEP.OUT.ATTRS = FALSE)[-i0, ]
  # clear irrelevant
  if(type == "circle" | sqexp) {
    # clip to radius
    rcn <- rcn[(rcn[,1]^2+rcn[,2]^2) <= (range^2), , drop =FALSE]
  }
  if(sqexp) {
    #
    D2 <- (rcn[,1]^2+rcn[,2]^2)
    # set 1 at distance 1
    W  <- exp(-(D2-1)/range^2)
    # clip corners?
  }
  #  browser()
  # focal rc
  rc <- tj_i2rc(i, nr, nc)
  #browser()
  o <- lapply(1:nrow(rc), \(j) {
    # shift template
    rcj <- cbind(rcn[,1] + rc[j,1], rcn[,2] + rc[j,2])
    # crop to domain
    good <- (rcj[,1] > 0) & (rcj[,1] <= nr) & (rcj[,2]>0) & (rcj[,2] <= nc)
    out <- tj_rc2i(rcj[good,], nr, nc)
    if(sqexp) attr(out, "weight") <- W[good]
    out
  })
  if(length(o) == 1) o[[1]]
  else o
}

#' Cell to Row-Column
#'
#' @param i cell
#' @param nr rowss
#' @param nc colums
#'
tj_i2rc <- function(i, nr, nc) {
  r  <- ceiling( i/nc )
  cbind(r = r, c = i - (r-1)*nc)
}

#' Row-Col to Cell
#'
#' @param i cell
#' @param nr rowss
#' @param nc colums
#'
tj_rc2i <- function(rc, nr, nc) {
  rc <- rbind(rc)
  (rc[,1]-1) * nc + rc[,2]
}



#' Raster Cell Neighbourhoods
#'
#' @param i focal cell
#' @param nr number of rows in the xy-grid
#' @param nc number of columns in the xy-grid
#' @param ... ignored
#' @details
#' Non-toroidal neighborhood, so at edges there will be less neighbours.
#'
#' Will assume cell-numbering
#' cell = col + (row-1) * ncol.
#' @returns
#' A vector of indices. If `type='sqexp'`, will return the weights in the attribute 'weight'.
#'
#' @export
tj_cell_neighbours_old <- function(i, nr, nc, ...) {
  rc <- tj_i2rc(i, nr, nc)


  # template for row-column indices:
  rcn <- cbind(c(-1, 0, 1, -1, 1, -1, 0, 1),
               c(-1, -1, -1, 0, 0, 1, 1, 1))


  o <- lapply(1:nrow(rc), \(j) {
    rcj <- cbind(rcn[,1] + rc[j,1], rcn[,2] + rc[j,2])
    rcj <- rcj[ rcj[,1] > 0 & rcj[,1] <= nr & rcj[,2]>0 & rcj[,2] <= nc, ]
    oj <- tj_rc2i( rcj, nr, nc)
    oj[oj > 0 & oj <= (nr*nc)]
  })
  if(length(o)) o[[1]]
  else o
}
