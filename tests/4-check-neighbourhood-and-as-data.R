# Check cell neighbours works with teh as-data

library(devtools)
load_all()

o <- load("data/testcubes_jump.rda") |> get()
x0 <- o$noise

d <- tj_stars_to_data(x0)

gd <- attr(d, "grid")




par(mfrow=c(2,1))
plot( d[d$time == 1, c("c_x", "c_y") ] , asp = 1)

iv <- c(topleft = 1, topright = gd['ncol'], bottomleft = gd['ncol']*(gd['nrow']-1)+1, bottomright = prod(gd),
        somewhere = 315)

for(i in iv){
  ni <- tj_cell_neighbours(i , nr = max(d$row), nc = max(d$col))
  points(d[ d$cell %in% c(ni) & d$time == 1, c("c_x", "c_y")], col = 1, pch = 19)
  points(d[ d$cell %in% c(i) & d$time == 1, c("c_x", "c_y")], col =  which(i==iv)+1, pch = 19)

}
