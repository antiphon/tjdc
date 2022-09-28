#' Process input stars to legacy data format
#'
#' Convert stars object into a tibble with variable names recogniced by the fitters.
#'
#' @param x stars object
#' @param timevar which dimension is time
#' @param attrvar which attribute to keep
#'
#' @details
#' This is essentially a wrapper for the as_tibble-method of stars-objects. The difference
#' is that only one attribute will be kept, chosen with `attrvar`. We assume that
#' spatial dimensions have names 'x' and 'y', and time dimension is defined by `timevar`.
#'
#' The function adds a cell number in the same way `terra::as.data.frame` does:
#' Spatially top-to-bottom, then left to right. That means cell=1 is for the cell with
#' coordinate (min(x), max(y)), and last cell is (max(x), min(y)).
#' @return
#' A tibble with variables `c_x`, `c_y` for spatial coordinates and `time` for temporal coordinates,
#' `col,row` for  the top-bottom - left-right indices in the data matrix,
#' `cell` as described above, and `value` for the values of the chosen attribute.
#'
#' The output also has an attribute grid=c(nrow, ncol).
#'
#' @import dplyr
#' @import stars
#' @export

tj_stars_to_data <- function(x, timevar = "z", attrvar = "values") {
#
  dims <- dim(x)
  if(!all(c("x", "y", timevar) %in% names(dims)))
    stop(sprintf("Dimension (x,y,%s) not found.", timevar))

  dat <- as_tibble(x) |>
    transmute(c_x     = x,
              c_y     = y,
              col     = as.factor(x) |> as.integer(),
              row     = as.factor(max(y)  -  y) |> as.integer(),
              cell    = col + (row-1) * dims['x'],
              time    = get(timevar),
              value   = get(attrvar),
              )
  # att grid info
  attr(dat, "grid") <- c(nrow = max(dat$row), ncol = max(dat$col))
  dat
}
