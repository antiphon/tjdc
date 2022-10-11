# Test m0.6 fitter

library(devtools)
load_all()

o <- load("data/testcubes_jump.rda") |> get()
x0 <- o$noise
x  <- x0
# as dataframe
jumpd <- tj_stars_to_data(o$jump, timevar = "z") |> filter(time == 11)
d <- tj_stars_to_data(x) |>
  left_join(jumpd |> transmute(cell, hasjump = value != 0))


# Example data frame
datx <- d |> filter(hasjump & cell %in% c(326:363))



da <- tj_fit_m0.6(x, timevar = "z",
                  attrvar = "values",
                  gamma = 0.2,
                  prior_k = 0.9,
                  verbose = TRUE,
                  niter = 1000,
                  keep_hist = TRUE,
                  prior_delta = list(m = 0,
                                     s2 = 10, a = -Inf, b = 0)
)
# or
da <- tj_fit_m0.6s(x, timevar="z",
                   cfg = tj_cfg_m0.6(gamma = 0.2, prior_k = 0.9, niter = 1000, keep_hist = TRUE, verbose=TRUE))




a <- da
## Trace
k <- 177
z <- posterior::as_draws(tj_trace_m0.3(a, cell = 2))
bayesplot::mcmc_trace(z)
#bayesplot::mcmc_dens(z)
posterior::summarise_draws(z)


#
library(ggplot2)

a <- da
s <- tj_summarise_m0.3(da, thin_steps = 1, burnin = .7)

exc <- sample(da$cell_info$cell, 5)

ds <- s |> left_join(d) |> filter(cell %in% exc)

pred <- tj_predict_m0.3(s = s, cells = exc)

print( ds |> ggplot() +
         geom_line(aes(time, value, group= cell, alpha = pred_jump_prob)) +
         geom_line(data = pred, aes(time, value), col = 2) +
         scale_alpha_continuous(limits = 0:1) +
         facet_wrap(~cell) + theme_bw() )



