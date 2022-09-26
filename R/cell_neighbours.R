#' Raster Cell Neighbourhoods
#'
#' @param i focal cell
#' @param nr number of rows in the xy-raster
#' @param nc number of columns in the xy-raster
#'
#' @details
#' Non-toroidal.
#'
#' NOTE we assume that the cell indices go first row then col, i.e. cell=1,2,3... are (row,col) (1,1), (1,2), (1,3)...
#' AND  we assume that rows go in "bottom to top" fashion in persieved y-dimension. Make sure data is like this! Stars::as_tibble
#' might have negative dimension deltas.
#'
#'
tj_cell_neighbours <- function(i, nr, nc) {

  rc2i <- function(rc, nr, nc) {
    rc <- rbind(rc)
    (rc[,1]-1) * nc + rc[,2]
  }
  i2rc <- function(i, nr, nc) {
    r  <- ceiling( i/nc )
    cbind(r=r, c = i - (r-1)*nc)
  }
  # target (c,r)
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

