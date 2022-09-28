#' Raster Cell Neighbourhoods
#'
#' @param i focal cell
#' @param nr number of rows in the xy-grid
#' @param nc number of columns in the xy-grid
#' @param ... ignored
#'
#' @details
#' Non-toroidal, 8 nearest aka queen neighbours.
#'
#' Will assume cell-numbering
#' cell = col + (row-1) * ncol.
#'
#' @export
tj_cell_neighbours <- function(i, nr, nc, ...) {

  rc2i <- function(rc, nr, nc) {
    rc <- rbind(rc)
    (rc[,1]-1) * nc + rc[,2]
  }
  i2rc <- function(i, nr, nc) {
    r  <- ceiling( i/nc )
    cbind(r = r, c = i - (r-1)*nc)
  }
  rc <- i2rc(i, nr, nc)
  rcn <- cbind(c(-1, 0, 1, -1, 1, -1, 0, 1),
               c(-1, -1, -1, 0, 0, 1, 1, 1))
  o <- lapply(1:nrow(rc), \(j) {
    rcj <- cbind(rcn[,1] + rc[j,1], rcn[,2] + rc[j,2])
    rcj <- rcj[ rcj[,1] > 0 & rcj[,1] <= nr & rcj[,2]>0 & rcj[,2] <= nc, ]
    oj <- rc2i( rcj, nr, nc)
    oj[oj > 0 & oj <= (nr*nc)]
  })
  if(length(o)) o[[1]]
  else o
}
