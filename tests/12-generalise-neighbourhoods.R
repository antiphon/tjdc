# Dev the generalised neighbourhood calculator

# cell = col + (row-1) * ncol.


library(devtools)
load_all()

if(!exists("nc")) {
  nc <- 35
  nr <- 20
  xy <- expand.grid(c=1:nc, r=1:nr)#[,2:1]
}

##############################################
# The neighbourhood calculator generalisation
#
# check speed
if(0) {
  i <- nc+36+5+35*2
  mb <- microbenchmark::microbenchmark(
    old = n1 <- tj_cell_neighbours_old(i, nr, nc),
    new = n2 <- tj_cell_neighbours(i, nr, nc, range = 2, type = "sqexp"),
  times = 20)

  print(mb)
  plot(xy, asp = 1, col = "gray50")
  points(xy[i,,drop=FALSE], col = 3, pch = 19)
  points( tj_i2rc(n2, nr, nc)[,2:1], col = 2, pch = 3)
}

# check type and range

if(0) {
  iv <- c(10 * nc + c(5, 15, 30))
  n  <- list()
  for(i in seq_along(iv))
    n[[i]] <- tj_cell_neighbours(iv[i], nr, nc, range = r <- 1 * i, type = t <- "sqexp")
  n2 <- c(unlist(n))
  w  <- unlist( lapply(n, \(v) attr(v, "weight")) )
  if(is.null(w)) w <- 1

  par(bg = "black")
  plot(xy[-n2,], asp = 1, pch = 19, col = "gray50")
  points(xy[iv,,drop=FALSE], col = 3, pch = 19)
  points( tj_i2rc(n2, nr, nc)[,2:1], col = 2, pch = 15, cex = w*1.5)
}

### Check normalisation to 8
if(1) {
  iv <- c(10 * nc + c(5, 15, 30))
  n  <- list()
  for(i in seq_along(iv))
    n[[i]] <- tj_cell_neighbours(iv[i], nr, nc, range = r <- 1 * i, type = t <- "sqexp", scale_to_8 = TRUE )
  n2 <- c(unlist(n))
  wl <- lapply(n, \(v) attr(v, "weight"))
  sapply(wl, sum)
  w  <- unlist( wl )
  if(is.null(w)) w <- 1

  par(bg = "black")
  plot(xy[-n2,], asp = 1, pch = 19, col = "gray50")
  points(xy[iv,,drop=FALSE], col = 3, pch = 19)
  points( tj_i2rc(n2, nr, nc)[,2:1], col = 2, pch = 15, cex = w*1.5)
}





