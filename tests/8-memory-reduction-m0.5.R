# New architecture to save memory on large mosaics
#
# Reserve memory only for the tracked cells in m0.3.

# Test m0.5 fitter

library(devtools)
load_all()

o <- load("data/testcubes_jump.rda") |> get()
r <- o$noise
# add some NA's

r$values[]# <- 1


subs <- tj_divide_raster(r, buffer = 1, n = c(3,2))
#

cfg <- tj_cfg_m0.3(                      verbose = TRUE,
                                         niter = 10,
                                         ncores = 1,
                                         keep_hist = TRUE,
                                         just_fit = TRUE
)

set.seed(1)

af <- tj_fit_m0.3_dac(r,
                      timevar = "z",
                      attrvar = "values",
                      subsets = subs$cell,
                      verbose = TRUE,
                      niter = 100,
                      ncores = 1,
                      keep_hist = 1,
                      just_fit = TRUE
)

object.size(af)



####################################################################################################
cfg <- tj_cfg_m0.3(                      verbose   = !TRUE,
                                         niter     = 100,
                                         ncores    = 1,
                                         keep_hist = TRUE,
                                         just_fit  = TRUE
)

tmpf <- tempdir()
load_all()

set.seed(1)
bf <- tj_fit_m0.5(r, timevar = "z", attrvar = "values", tiling = subs, cfg = cfg)


####################################################################################################
s <- tj_summarise_m0.3(bf[[1]])
plot(s |> select(c_x, c_y, pred_jump_prob) |> st_as_stars(dims =c("c_x", "c_y")))

