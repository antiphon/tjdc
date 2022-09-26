# Check cell neighbours works with teh as-data

library(devtools)
load_all()

o <- load("data/testcubes_jump.rda") |> get()
x0 <- o$noise

d <- tj_stars_to_data(x0)

i <- 80

ni <- tj_cell_neighbours(i , nr = max(d$row), nc = max(d$col))

par(mfrow=c(2,1))

plot( d[d$x == 1, c("c_x", "c_y") ] , asp = 1)
points(d[ d$cell %in% c(ni) & d$x == 1, c("c_x", "c_y")], col = 2, pch = 19)
points(d[ d$cell %in% c(i) & d$x == 1, c("c_x", "c_y")], col = 1, pch = 19)

