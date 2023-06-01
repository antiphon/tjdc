# Scale the neighoburhood to sum to 8 for sacling invariance.
#

library(devtools)
load_all()

if(!exists("x")) {
o <- load("data/testcubes_jump.rda") |> get()
x0 <- o$artrend
x  <- x0 |> slice("x", 10:25) |> slice("y", 5:15)
jdat <- tj_stars_to_data(o$jump[,,,10] )
xdat <- tj_stars_to_data(x) |>
  left_join(jdat |> transmute(c_x, c_y, jump = value != 0))
}
mos <-  tj_divide_raster(x, n = c(2,2), buffer = c(2,2))


# Example data frame
mfl <- tibble()
tims <- list()
for(typ in c("square", "circle", "sqexp")[1] )
  for(nr in 2)
    for(scaleit in c(T,F)){
    tim <- system.time(
      mf <- tj_fit_m0.x(x,
                        timevar = "z",
                        attrvar = "values",
                        gamma = 0.2,
                        tiling = mos,
                        prior_k = 0.5,
                        verbose = TRUE,
                        niter = 100,
                        neighbourhood_type = typ,
                        neighbourhood_range = nr,
                        neighbourhood_scale_to_8 = scaleit,
                        keep_all_in_memory = TRUE
      )
    )
    tims[[paste(typ, nr, scaleit)]] <- tim
    me <- tj_stitcher_m0.x(mf, mos)
    z  <- tjdc:::tj_summarise_m0.5(me)
    dp <- z |>
      select(c_x, c_y, pred_jump_prob) |>
      mutate(ntype = typ, nrange = nr, scale8 = scaleit)
    mfl <- bind_rows(mfl, dp)
  }
############################### Compare
mflx <- left_join(mfl, xdat |> filter(time == 1))

library(ggplot2)
mflx |> ggplot() +
  geom_raster(aes(c_x, c_y, fill = pred_jump_prob > 0.5)) + coord_fixed(exp=F) +
  facet_grid(nrange ~ ntype + scale8)

mflx |> group_by(ntype, nrange) |>
  summarise(fp = mean(pred_jump_prob[!jump] > 0.5),
       tp = mean(pred_jump_prob[jump] > 0.5))
# timings
as_tibble(tims) |> t()



# eof
