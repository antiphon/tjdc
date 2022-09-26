# Generate test stack for jump analysis

library(stars)
library(dplyr)

# Stacks 1.*: small
nc <- 35
nr <- 20
nt <- 11  # cut to length
Nx <- nc * nr
N <- nc * nr * nt
bb <- st_bbox(c(xmin = 0, xmax = nc/nr,ymin = 0, ymax=1, zmin=0, zmax=20))
#

s0 <- stars::st_as_stars(bb, nx = nc, ny = nr, nz = nt, n = N)
st_crs(s0)  <- st_crs(3067)
r0 <- s0[,,,1, drop=TRUE]
xy <- st_coordinates(r0)[,-3] |> as.matrix()

# Mask for trend
D <- colSums((t(xy) - c(.5,.5))^2)
them <- which(D <= 0.2^2)
# Mask for jump
JD <- colSums((t(xy) - c(.7,.5))^2)
Jthem <- which(JD <= 0.2^2)
# Jump serieses
jump  <- s0
jump$values <- 0
jthem_rc <- stars:::colrow_from_xy(xy[Jthem,], jump)
for(i in 1:nrow(jthem_rc)) jump$values[jthem_rc[i,1], jthem_rc[i,2],] <- c(rep(0,floor(nt/2) ), rep(-3, ceiling(nt/2) ))
# Mask for NA
c0 <- c(.7, 0)
themNA <- which(colSums((c0-t(xy))^2) < 0.2^2 )
mask <- r0 |> mutate(values = 1, values = replace(values, themNA, NA))

##
s <- 0.1
# just noise
v1n <- rnorm(N, 0, s)
s0n <- s0 |> mutate(values = v1n)
s1n <- s0n * mask + jump


# noise with strong AR(1) noise
v1 <- v1n
A <- Matrix::bandSparse(nt, k = 0:1, diagonals = list(rep(1, nt), rep(0.5, nt)), symmetric=TRUE)
U <- Matrix::chol(A)

s1a <- st_apply(s0n, c("x","y"), \(v) as.numeric(U%*%v), .fname="z") |> aperm(c(2,3,1))
s1a <- s1a * mask + jump

# noise with trendy area in the middle
D <- colSums((t(xy) - c(.4,.5))^2)
them <- which(D <= 0.2^2)
hastrend <- mask * (r0 |> mutate(values = replace(values, them, 1)))

s1t <- s1n
them_rc <- stars:::colrow_from_xy(xy[them,], s1t)
for(i in 1:nrow(them_rc)) s1t$values[them_rc[i,1], them_rc[i,2],] <- s1t$values[them_rc[i,1], them_rc[i,2],] + 1:nt*.3/nt



# AR + trend
s1at <- s1a
for(i in 1:nrow(them_rc)) s1at$values[them_rc[i,1], them_rc[i,2],] <- s1at$values[them_rc[i,1], them_rc[i,2],] + 1:nt*.3/nt
#plot(s1at[k][1,])

###
out <- list(noise=s1n, ar=s1a, trend=s1t, artrend=s1at, hastrend = hastrend, jump = jump)

#### store
test_stacks3 <- out
save(test_stacks3, file="data/testcubes_jump.rda", compress = TRUE)



