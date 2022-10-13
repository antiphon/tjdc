# Test m0.x Fitter

library(devtools)
load_all()

o <- load("data/testcubes_jump.rda") |> get()
x0 <- o$noise
r  <- x0



tiling <- tj_divide_raster(r, buffer = 2)

ffl <- list()
for(m in c("m0.3", "m0.6")) {
  ffl[[m]] <- tj_fit_m0.x(r, timevar = "z",
                  tiling = tiling, gamma = .2,
                  keep_hist = c(1,30),
                  niter = 1000,
                  prior_sigma2 = c(10, 1),
                  prior_delta = list(m = 0, s2 = 1e5, a = -Inf, b = 0),
                  model_variant = m, verbose=TRUE,
                  verbose2 = TRUE, ncofl = 4)
  print(ffl[[m]][[1]]$took)
}
library(ggplot2)

s <- lapply(names(ffl), \(m) tj_summarise_m0.3(ffl[[m]][[i <- 1]]) |> mutate(mod = m)) |> bind_rows()
ggplot(s) + geom_raster(aes(c_x, c_y, fill = pred_jump_prob)) + coord_fixed(ex=F) + facet_wrap(~mod)



e <- lapply(names(ffl), \(m) tj_trace_m0.3(ffl[[m]][[i <- 1]], 1))

library(posterior)
mc <- as_draws(e)
bayesplot::mcmc_trace(mc)

