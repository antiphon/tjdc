# Test m0.5 fitter

library(devtools)
load_all()

o <- load("data/testcubes_jump.rda") |> get()
x0 <- o$noise
r  <- x0

datx <- tj_stars_to_data(r)
cell_info <- datx |> group_by(cell, col, row, c_x, c_y) |> summarise(isna = is.na(sum(value)) ) |> ungroup()

subs <- tj_divide_bbox(st_bbox(r), expand = 0.1)
xy <- st_intersects(r, subs, as_points=F, sparse=FALSE)   #|> filter(z == 1) |> arrange(y, x)
subsets <- apply(xy, 2, which)
#




af <- tj_fit_m0.3_dac(r,
                 timevar = "z",
                 attrvar = "values",
                 gamma = 0.1,
                 prior_k = 0.5,
                 subsets = subsets,
                 verbose = TRUE,
                 niter = 1000,
                 ncores = 8,
                 keep_hist = TRUE,
                 just_fit = TRUE
                 )

a <- tj_stitcher_m0.5_v1.1(af, subsets)
#
##############################################
# Check helpers

# Summarise
s <- summarise_m0.5(a)


library(ggplot2)
s |>
 ggplot() + geom_raster(aes(c_x, c_y, fill = pred_jump_prob )) +
  scale_fill_viridis_c(option = "turbo") + theme_void() +
  #geom_sf(data = subs, aes(col = id), fill = NA) +
  coord_sf(expand=F)


set.seed(1)
excel <- (cell_info |> filter(!isna)|>sample_n(20))$cell

predx <- predict_m0.3(s = s, cells = excel)

dpx <- datx |>filter(cell %in% excel)
#plot(a$hist_sigma2)
print( dpx |>
         ggplot() +
         geom_line(aes(time, value, group= cell)) +
         geom_line(data = predx, aes(time, value), col = "green") +
         scale_alpha_continuous(limits = 0:1) +
         facet_wrap(~cell) + theme_bw() )


####### Check particular cell's history:
cell <- 605

th <- trace_m0.5(af, cell)
dr <- posterior::as_draws(th)
bayesplot::mcmc_trace(dr)


