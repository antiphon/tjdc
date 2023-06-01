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
mf <- function(f) f(x,
                    timevar = "z",
                    attrvar = "values",
                    gamma = 0.2,
                    prior_k = 0.5,
                    verbose = !TRUE,
                    niter = 300,
                    neighbourhood_type = typ,
                    neighbourhood_range = nr,
                    neighbourhood_scale_to_8 = FALSE
)

tims <- rbind(

  v1 = system.time( mf(tj_fit_m0.6_v1)),
  v1.1 = system.time( mf(tj_fit_m0.6) )



)

############################### Compare
print(tims)



# eof
