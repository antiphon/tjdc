# Debug, something with the mosaic not being correct

library(devtools)
load_all()

o <- load("data/testcubes_jump.rda") |> get()
r <- o$noise
r <- st_set_dimensions(r, "x", values = st_get_dimension_values(r, "x")+1000)
r <- st_set_dimensions(r, "y", values = st_get_dimension_values(r, "y")+10000)
subs <- tj_divide_bbox(st_bbox(r), expand = 0.1)
xy <- st_intersects(r, subs, as_points=F, sparse=FALSE)   #|> filter(z == 1) |> arrange(y, x)
subsets <- apply(xy, 2, which)
#




af <- tj_fit_m0.3_dac(r,
                 timevar = "z",
                 attrvar = "values",
                 subsets = subsets,
                 niter = 10,
                 ncores = 8,
                 keep_hist = TRUE,
                 just_fit = TRUE
                 )

# Summarise

f <- af[[1]]
s <- summarise_m0.3(f)

library(ggplot2)
library(sf)
s |>
 ggplot() +
  geom_sf(data = st_bbox(r) |> st_as_sfc()) +
  geom_raster(aes(c_x, c_y, fill = pred_jump_prob )) +
  geom_sf(data = subs[9,], aes(col = id), alpha = .3) +
  geom_sf(data = subs, aes(col = id), alpha = .3) +
  coord_sf(expand=F) + theme_void()

