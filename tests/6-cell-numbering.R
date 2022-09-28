# Cell number mapping between raster, terra, and stars


devtools::load_all()

library(stars)

# Start with something


library(dplyr)
m <- expand.grid(x=1:9, y = 1:3) |> mutate(val = x^2 + y, z = 0)

#

rs <- stars::st_as_stars(m, dims = c("x", "y", "z"), xy = c("x", "y"))
rs1 <- rs[,,1,drop=T]
# terra: tp - lr, spatial top
library(terra)
td <- terra::rast(rs1) |>
  as.data.frame(cell = TRUE, xy = TRUE) |>
  as_tibble()

# stars: tp - lr
ts <- as_tibble(rs)

# raster: tp - lr
library(raster)
tr <- raster(rast(rs1)) |> as.data.frame(cell = TRUE, xy = TRUE) |>
  as_tibble()

# mine: tp -lr now
tt <- tj_stars_to_data(rs, attrvar = "val", timevar = "z")

# get some, and check same coord and val in others
them <- sample(tt$cell, 5)
ex   <- tt |> filter(cell %in% them)

