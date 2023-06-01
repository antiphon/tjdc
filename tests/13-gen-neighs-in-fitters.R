# Include generalised neighbourhood into model fitters
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
# as dataframe

# Example data frame
mfl <- tibble()
tims <- list()
for(typ in c("square", "circle", "sqexp")[c(1,3)] )
  for(nr in 1:2){
    tim <- system.time(
      mf <- tj_fit_m0.6(x,
                        timevar = "z",
                        attrvar = "values",
                        gamma = 0.2,
                        prior_k = 0.5,
                        verbose = TRUE,
                        niter = 1000,
                        neighbourhood_type = typ,
                        neighbourhood_range = nr,
                        neighbourhood_scale_to_8 = TRUE
      )
    )
    tims[[paste(typ, nr)]] <- tim
    z  <- tjdc:::tj_summarise_m0.3(mf)
    dp <- z |>
      select(c_x, c_y, pred_jump_prob) |>
      mutate(ntype = typ, nrange = nr)
    mfl <- bind_rows(mfl, dp)
  }
############################### Compare
mflx <- left_join(mfl, xdat |> filter(time == 1))

library(ggplot2)
mflx |> ggplot() +
  geom_raster(aes(c_x, c_y, fill = pred_jump_prob > 0.5)) + coord_fixed(exp=F) +
  facet_grid(nrange ~ ntype)

mflx |> group_by(ntype, nrange) |>
  summarise(fp = mean(pred_jump_prob[!jump] > 0.5),
       tp = mean(pred_jump_prob[jump] > 0.5))
# timings
as_tibble(tims) |> t()



# eof
