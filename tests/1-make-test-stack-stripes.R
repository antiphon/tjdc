# Generate test stack for jump analysis
# artificial stripes, and also make time-reverse of it
library(stars)
library(dplyr)

# Stacks stripes.*: small
nr <- 6
nt <- 11  # cut to length
nc <- 10 * nt + 2 # add to horiz. ends 1pixel buffer
Nx <- nc * nr
N <- nc * nr * nt
bb <- st_bbox(c(xmin = -1, xmax = nc-1, ymin = 0, ymax=nr, zmin=1, zmax=nt))
#

s0 <- stars::st_as_stars(bb, nx = nc, ny = nr, nz = nt, n = N)
st_crs(s0)  <- st_crs(3067)
r0 <- s0[,,,1, drop=TRUE]
xy <- st_coordinates(r0)[,-3] |> as.matrix()

# Jump times.
jump  <- r0 * 0 + 1
jump$values <- nt
for(i in 1:nt) jump$values[((i-1)*10+ifelse(i>1, 2, 1)):(i * 10+1),] <-  i
# time reversed
jump2 <- nt - jump + 1

##
# draw clear lines
# with just noise
##
sig <- 10
s2 <- s1  <- s0
alpha <- 100
deltaf <- -0.3
set.seed(1)
for(i in 1:nt) {
  v <- rnorm(nc*nr, 0, sig) + alpha
  s1$values[,, i]           <-  v + deltaf * alpha * (jump$values < i)
  s2$values[,, nt-i+1]      <-  v + deltaf * alpha * (jump2$values < (nt-i+1) )
}
#
#### store
easystripes <- list(strip = s1, strip_rev = s2, jump = jump, jump_rev = jump2)
save(easystripes, file="data/easystripes.rda", compress = TRUE)





### check
if(0) {

  x1 <- tjdc::tj_stars_to_data(s1)
  x2 <- tjdc::tj_stars_to_data(s2)
  bind_rows(x1 |> mutate(t="f"),
            x2 |> mutate(t ="b") )|> filter(cell %in% (1:4) ) |>
  ggplot() +
    geom_point(aes(time, value, col = t)) + facet_wrap(~cell)

}
if(0){
  plot(jump2)
}
