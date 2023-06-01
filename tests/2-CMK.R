# Test the mann-kendall calls,
library(devtools)

load_all()

o <- load("data/test_stacks2.rda") |> get()
x0 <- o$trend

x <- x0 |> st_as_stars()
names(x) <- "values"
x <- stars::st_set_dimensions(x, "band", 1:20)

######
# Check if the pointwise non-contextuals are equal

if(0){
fun <-function(v,...) if(sum(is.na(v))>1) return(NA) else tj_mann_kendall(v, est_beta = FALSE)$p

t0 <- system.time(   a <- st_apply(x, c(1,2), fun, .fname = "p")$p )
t2 <- system.time(   b <- (s <- tj_contextual_mann_kendall(x, neighbourhood = 0))$p )

print(all.equal(a,b))
print(rbind(t0,t2))
}



# Check forward indexing, now something is different...
if(0) {
  d  <- dim(x)
  nc <- d['x'] # how stars does it
  nr <- d['y']

  xy <- as_tibble(x)

  # some examples
  iv <- c(1,
          4+1 + 9 * nc, # fourth of 10th row
          nc*nr) # last cell
  plot(xy[,1:2], asp = 1)
  for(i in iv){
    a <- forward_neighbour_cells_queen_col_row(i - 1, nr, nc) + 1 # zero indexing
    points(xy[i, 1:2, drop=F], pch = 3, col = ci <- which(i==iv))
    points(xy[a, 1:2], pch = 19, col = ci)
  }
}




## Make sure matrix cells are in right order for the the previously developed code
if(0) {
  o <- load("data/test_stacks2.rda") |> get()
  x0 <- o$trend
  x <- x0 |> st_as_stars()
  x <- stars::st_set_dimensions(x, "band", 1:20)
  names(x0) <- 1:20
  # old
  o <- load("data/test_stacks2.rda") |> get()
  x0 <- o$trend
  X0v <- raster::as.data.frame(x0, xy = TRUE)
  X0 <- raster::values(x0)
  x  <- st_as_stars(x0)

  # new
  Xv <- X <- as_tibble(x)
  X <- Xv |>
    # Make a matrix, cells x time
    pivot_wider(values_from = layer.1, names_from = band) |>
    select(-x, -y) |>
    as.matrix()
  all.equal(X0, X, check.attributes = FALSE)
  # calc wise
  res0 <- ConMK:::c_contextual_mann_kendall( X0,
                                    nrow(x0),
                                    time = 1:20,
                                    neighbourhood = 2,
                                    calc_slope = 0)
  res <- c_contextual_mann_kendall( X,
                                    ncol(x),
                                    time = 1:20,
                                    neighbourhood = 2,
                                    calc_slope = 0)
  all.equal(res0$S_and_s2, res$S_and_s2, check.attributes = FALSE)

}




# Then check if the smoother works
if(0) {
  o <- load("data/test_stacks2.rda") |> get()
  x0 <- o$trend
  x <- stars::st_as_stars(x0) |> stars::st_set_dimensions("band", values = 1:20)
  #x <- o$artrend
  # manipulate time
  #
  t0 <- system.time(   a <- (s0 <- ConMK::contextual_mann_kendall(x0, neighbourhood = 2))$p[] )
  t2 <- system.time(   b <- (s1 <- tj_contextual_mann_kendall(x, neighbourhood = 2, attrvar = "layer.1", timevar = "band"))$p |> c())
  #
  print( all.equal(a,b) )
  p0 <- s0$p |> st_as_stars()
  p1 <- s1['p']
  plot(p1-p0)
  (c(p0, p1)  < .05)  |> st_redimension() |> plot()
  print(rbind(t0, t2))
   # ok
}


# Then check if the smoother works, 2
if(1) {
  o <- load("data/testcubes_trend.rda") |> get()
  x <- o$trend
  x0 <- o$trend |> stars:::st_as_raster("raster")
  #
  t0 <- system.time(   a <- (s0 <- ConMK::contextual_mann_kendall(x0, neighbourhood = 2))$p[] )
  t2 <- system.time(   b <- (s1 <- tj_contextual_mann_kendall(x, neighbourhood = 2))$p |> c())
  #
  print( all.equal(a,b) )
  p0 <- s0$p |> st_as_stars()
  # some bug, "not identical"
  p1 <- p0 |> mutate(p = s1['p']$p |> c())
  plot(p1-p0)
  (c(p0, p1)  < .05)  |> st_redimension() |> plot()
  print(rbind(t0, t2))
  # ok
}


