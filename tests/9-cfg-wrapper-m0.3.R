# Test the cfg-wrapper for m0.3

library(devtools)
load_all()

o <- load("data/testcubes_jump.rda") |> get()
r <- o$noise

set.seed(1)
da <- tj_fit_m0.3(r, timevar = "z",
                  attrvar = "values",
                  gamma = 0.0,
                  prior_k = 0.9,
                  verbose = TRUE,
                  niter = 10,
)

load_all()
set.seed(1)
db <- tj_fit_m0.3s(x=r, timevar="z", attrvar="values",
                   cfg = tj_cfg_m0.3(niter = 10,
                                     gamma = 0.0, prior_k = .9, verbose = TRUE))

all.equal(da, db)
