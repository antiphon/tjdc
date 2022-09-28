# Check tiling

# Debug, something with the mosaic not being correct

library(devtools)
load_all()

o <- load("data/testcubes_jump.rda") |> get()
x <- o$noise



##############################################################################
#
s <- tj_divide_raster(x, n = c(5,2), buffer = 3)
