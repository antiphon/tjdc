# Test m0.5 fitter

library(devtools)
load_all()

o <- load("data/testcubes_jump.rda") |> get()
x0 <- o$noise
x  <- x0

datx <- tj_stars_to_data(x)
cell_info <- datx |> filter(x == x[1]) |> mutate(isna = is.na(y), x=NULL, y = NULL)

subs <- tj_divide_bbox(st_bbox(x), expand = 0.1)
xy <- st_intersects(x, subs, as_points=F, sparse=FALSE)   #|> filter(z == 1) |> arrange(y, x)
subsets <- apply(xy, 2, which)
#




af <- tj_fit_m0.3_dac(x,
                 timevar = "z",
                 attrvar = "values",
                 gamma = 0.1,
                 subsets = subsets,
                 prior_k = 0.5,
                 verbose = TRUE,
                 niter = 30,
                 ncores = 8,
                 keep_hist = TRUE,
                 prior_theta = list(m = c(0, 0, 0),
                                    S = diag( c(1e4, 1e4, 1e2) )),
                 truncate_jump_at_mean = 0,
                 just_fit = TRUE
                 #dat = datx
                 )

a <- tj_stitcher_m0.5_v1.1(af, subsets)

#
b <- cell_info |> mutate(prob = c(a$z) )
#
dp <- datx  |> left_join(b)

library(ggplot2)
dp |>
 ggplot() + geom_raster(aes(c_x, c_y, fill = prob)) + coord_fixed()



wra <- function(k) {
  be <- a$theta_m[,k]
  j  <- which.max(a$k[,k])
  xv <- unique(datx$x)
  kv <- rep(1, length(xv))
  kv[1:j] <- 0
  tibble(x = xv, y=c(cbind(1, xv, kv) %*% be) )
}
pred <- b |> filter(!isna) |> rowwise() |> summarise(wra(cell), cell = cell)

#plot(a$hist_sigma2)
print( dp |> ggplot() +
         geom_line(aes(x, y, group= cell, alpha = prob)) +
         geom_line(aes(x, y), col = 2, data = pred) +
         scale_alpha_continuous(limits = 0:1) +
         facet_wrap(~cell) + theme_bw() )

k <- 297
plot(a$hist_theta[,k,3])#


