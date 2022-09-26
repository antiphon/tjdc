#' Process input stars to legacy data format
#'
#' @param x stars object
#' @import dplyr
#' @import stars
#' @export

tj_stars_to_data <- function(x, timevar = "z", attrvar = "values") {
#
  dims <- st_dimensions(x)

  if(!all(c("x", "y", timevar) %in% names(dims))) stop("Dimension x,y not found.")

  dat <- as_tibble(x) |>
    transmute(c_x = x,
              col = as.factor(x) |> as.integer(),
              c_y = y,
              row = as.factor(y) |> as.integer(),
              cell = as.factor( sprintf("%010i_%010i",
                                #        max(row) + sign(dims$y$delta) * row,
                                #        max(col) + sign(dims$x$delta) * col)
                                row, col)
                                ) |> as.integer() ,
              x   = get(timevar),
              y   = get(attrvar),
              )
  # att grid info
  attr(dat, "grid") <- c(nrow = max(dat$row), ncol = max(dat$col))
  dat  |> arrange(cell, x) # x = time in below algorithm!
}
