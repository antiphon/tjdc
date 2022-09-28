# Test m0.5 fitter

library(devtools)
load_all()

o <- load("data/testcubes_jump.rda") |> get()
x0 <- o$noise
x  <- x0

datx <- tj_stars_to_data(x)
cell_info <- datx |> group_by(cell, col, row, c_x, c_y) |> summarise(isna = is.na(sum(value)) ) |> ungroup()

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
                 niter = 2000,
                 ncores = 8,
                 keep_history = TRUE,
                 prior_theta = list(m = c(0, 0, 0),
                                    S = diag( c(1e4, 1e4, 1e2) )),
                 truncate_jump_at_mean = 0,
                 just_fit = TRUE
                 )

a <- tj_stitcher_m0.5_v1.1(af, subsets)

#
b <- cell_info |> mutate(prob = c(a$z) )
#
dp <- datx  |> left_join(b)
#
library(ggplot2)
b |>
 ggplot() + geom_raster(aes(c_x, c_y, fill = prob )) +
  scale_fill_viridis_c(option = "turbo") +
  geom_sf(data = subs, aes(col = id), fill = NA) +
  coord_sf(expand=F)



wra <- function(k) {
  be <- a$theta_m[,k]
  j  <- which.max(a$k[,k])
  xv <- unique(datx$time)
  kv <- rep(1, length(xv))
  kv[1:j] <- 0
  tibble(time = xv, value=c(cbind(1, xv, kv) %*% be) )
}

pred <- b |> filter(!isna) |> rowwise() |> summarise(wra(cell), cell = cell)

set.seed(1)
excel <- sample(cell_info$cell, 20)

predx <- pred |> filter(cell %in% excel)

dpx <- dp |>filter(cell %in% excel)
#plot(a$hist_sigma2)
print( dpx |>
         ggplot() +
         geom_line(aes(time, value, group= cell)) +
         geom_line(data = predx, aes(time, value), col = "green") +
         scale_alpha_continuous(limits = 0:1) +
         facet_wrap(~cell) + theme_bw() )




####### Check particular cell:
k <- 105
list_of_fits <- af
cell_map <- lapply(seq_along(list_of_fits), \(si) list_of_fits[[si]]$cell_mapping |>
                     mutate(subset = si)) |>
  bind_rows()
fi <- cell_map |> filter(oldcell == k) |> slice()
f1 <- af[[fi$subset[1]]]
th <- f1$hist_theta[,fi$cell[1],] |> data.frame() |> setNames(c("a","b","d"))
th$z <-  f1$hist_z[,fi$cell[1]]
bayesplot::mcmc_trace(th)


