# Check tiling

# Debug, something with the mosaic not being correct

library(devtools)

o <- load("data/testcubes_jump.rda") |> get()
x <- o$noise
d <- tj_stars_to_data(x, "z", "values")


##############################################################################
#
load_all()
s <- tj_divide_raster(x, n = c(5,2), buffer = c(0,0) + 2, open.upper = TRUE)

p <- s$poly
re <- s$rc

c <- lapply(seq_along(s$cell), \(si) tibble( cell = s$cell[[si]], s = si)) |> bind_rows()
ci <- c |> group_by(cell) |> summarise(id = paste0(s, collapse=","))


##############
#
library(ggplot2)
library(sf)
# some tiles only
i <-   c(2,3,4,8) #1:nrow(p)
j <- sapply(i, grepl,ci$id) |> apply(1,any)
# Check coordinate space:

p1 <- ggplot() +
  geom_tile(data = d, aes(c_x, c_y), alpha = 1, fill = 1) +
  geom_point(data = d, aes(c_x, c_y), col = "white", alpha  =.2, size = .5) +
  ggnewscale::new_scale("fill") +
  geom_sf(data = p[i,], aes(fill = factor(id) ), col = NA, alpha = .7) +
  geom_sf_text(data = p[i,], aes(label = id ))


# Check row-col space

p2 <- ggplot(d) +
  geom_tile(aes(col, row), fill =  1) +
  geom_point(data = d, aes(col, row), col = "white", alpha  =.2, size = .5) +
  ggnewscale::new_scale("fill") +
  geom_rect(data = re[i,], aes(xmin = cmin, xmax = cmax, ymin = rmin, ymax=rmax, fill = factor(id) ), alpha = .6, col = NA)+
  coord_fixed() + scale_y_reverse()

## Cells
p3 <- ggplot(d |> left_join(ci) |> filter(cell %in% ci$cell[j]) ) +
  geom_tile(data = d, aes(col, row), alpha = 1, fill = 1) +
  geom_point(data = d, aes(col, row), col = "white", alpha  =.2, size = .5) +
  ggnewscale::new_scale("fill") +
  geom_point(aes(col, row, col = id), alpha = 1) +
  coord_fixed() + scale_y_reverse()

library(patchwork)

#print(p1/p2/p3 )
print(p1)

