# Test m0.3 fitter

library(devtools)
load_all()

o <- load("data/testcubes_jump.rda") |> get()
x0 <- o$noise
x  <- x0
# as dataframe
d <- tj_stars_to_data(x) |>
  left_join(tj_stars_to_data(o$jump) |> filter(x == x[11]) |> transmute(cell, hasjump = y != 0))


# Example data frame
datx <- d |> filter(hasjump & cell %in% c(326:363))



a <- tj_fit_m0.3(x, timevar = "z",
                 attrvar = "values",
                 gamma = 0.0,
                 prior_k = 0.9,
                 verbose = TRUE,
                 niter = 1000,
                 keep_hist = TRUE,
                 prior_theta = list(m = c(0, 0, 0),
                                     S = diag( c(1e4, 1e4, 1e2) )),
                 truncate_jump_at_mean = -1,
                 dat = datx)

#
b <- a$cell_info |> mutate(prob = a$z[a$cell_info$cell])
#
dp <- datx  |> left_join(b)

library(ggplot2)
#dp |>
#  ggplot() + geom_raster(aes(c_x, c_y, fill = prob)) + coord_fixed()



wra <- function(k) {
  be <- a$theta_m[,k]
  j  <- which.max(a$k[,k])
  xv <- unique(datx$x)
  kv <- rep(1, length(xv))
  kv[1:j] <- 0
  tibble(x = xv, y=c(cbind(1, xv, kv) %*% be) )
}
pred <- b |> rowwise() |> summarise(wra(cell), cell = cell)

#plot(a$hist_sigma2)
print( dp |> ggplot() +
  geom_line(aes(x, y, group= cell, alpha = prob)) +
  geom_line(aes(x, y), col = 2, data = pred) +
    scale_alpha_continuous(limits = 0:1) +
  facet_wrap(~cell) + theme_bw() )

k <- 297
#plot(a$hist_theta[,k,3])#
