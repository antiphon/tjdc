# Test m0.x fitter

library(devtools)
load_all()

o <- load("data/testcubes_jump.rda") |> get()
x0 <- o$noise
r  <- x0



tiling <- tj_divide_raster(r, buffer = 2)


# m0.3:
r3 <- tj_fit_m0.x(r, timevar = "z", tiling = tiling, model_variant = "m0.3", verbose=TRUE, verbose2 = TRUE, ncores = 1)

r6 <- tj_fit_m0.x(r, timevar = "z", tiling = tiling, model_variant = "m0.6", verbose=TRUE, ncores = 4)
