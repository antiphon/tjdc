#' Raster to Polygons
#' Take stars-raster and create a polygons of those pixels with >0 values using st_union
#' @param r stars object
#'
#' @return sf-object of polygons
#'
#' @import stars sf dplyr
#' @export

tj_stars_to_polygon_union <- function( r, which = 1) {
  if(length(names(r)) > 1 || length(dim(r)) > 2) stop("r should have only one layer and 2 dimensions.")
  p <- st_as_sf(r > 0)
  names(p)[1] <- "val"
  pt <- p |> filter( val == TRUE)
  pts <- pt |> st_union()
  # good enough actually
  pts
}

