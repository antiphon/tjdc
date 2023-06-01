# Generate test stack for trend analysis

library(stars)
library(dplyr)

# Stacks 1.*: small
nc <- 60
nr <- 40
nt <- 20  # cut to length
Nx <- nc * nr
N <- nc * nr * nt
bb <- st_bbox(c(xmin = 0, xmax = nc/nr,ymin = 0, ymax=1, zmin=0, zmax=20))
#
#r0 <- raster(nrows = nr, ncols = nc, crs = CRS("+init=epsg:3067"), xmn=0, xmx=1,
#             ymn = 0, ymx=nr/nc)

s0 <- stars::st_as_stars(bb, nx = nc, ny = nr, nz = nt, n = N)
st_crs(s0)  <- st_crs(3067)
r0 <- s0[,,,1, drop=TRUE]

##
s <- 0.1
# just noise
v1n <- rnorm(N, 0, s)
s0n <- s0 |> mutate(values = v1n) #mask(setValues(s0, v1n), mask)


# Mask for NA
xy <- st_coordinates(r0)[,-3] |> as.matrix()
R2 <- 0.2^2
c0 <- c(.7, 0)
themNA <- which(colSums((c0-t(xy))^2) < R2 )
mask <- r0 |> mutate(values = 1, values = replace(values, themNA, NA))
s1n <- s0n * mask

#
# noise with strong AR(1) noise
v1 <- v1n
A <- Matrix::bandSparse(nt, k = 0:1, diagonals = list(rep(1, nt), rep(0.5, nt)), symmetric=TRUE)
U <- Matrix::chol(A)

s1a <- st_apply(s0n, c("x","y"), \(v) as.numeric(U%*%v), .fname="z") |> aperm(c(2,3,1))
s1a <- s1a * mask

# noise with trendy area in the middle
D <- colSums((t(xy) - c(.5,.5))^2)
them <- which(D <= 0.11^2)
hastrend <- mask * (r0 |> mutate(values = replace(values, them, 1)))

s1t <- s1n
them_rc <- stars:::colrow_from_xy(xy[them,], s1t)
for(i in 1:nrow(them_rc)) s1t$values[them_rc[i,1], them_rc[i,2],] <- s1t$values[them_rc[i,1], them_rc[i,2],] + 1:nt*.3/nt



# AR + trend
s1at <- s1a
for(i in 1:nrow(them_rc)) s1at$values[them_rc[i,1], them_rc[i,2],] <- s1at$values[them_rc[i,1], them_rc[i,2],] + 1:nt*.3/nt
#plot(s1at[k][1,])

###
out <- list(noise=s1n, ar=s1a, trend=s1t, artrend=s1at, hastrend = hastrend)

#### store
test_stacks2 <- out
save(test_stacks2, file="data/testcubes_trend.rda", compress = TRUE)



