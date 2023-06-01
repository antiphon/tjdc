
#' Summarise fit of M0.5 stiched results
#'
#' @param x object fitted by this package's tj_fit_m0.5 and then stitched
#'
#' @import stars dplyr
#' @export
tj_summarise_m0.5 <- function(x, ..., cells) {

  if(missing(cells)) cells <- x$cells # do alla

  # for each cell only once
  cinfo <- x$cell_info |>
    group_by(cell, c_x, c_y) |>
    summarise(in_subsets = n()) |>
    ungroup()
  # just gather
  pred <- cinfo |>
    filter(cell %in% cells) |>
    mutate(
      theta_mean = t(x$theta_m[,cell]) |> as_tibble() |>setNames(c("a","b","d")),
      #theta_sd  = t(x$) |> as_tibble() |> setNames(c("a","b","d")),
      pred_jump_prob = x$z[cell],
      pred_jump_mode = apply(x$k[,cell,drop=FALSE], 2, which.max),
      pred_jump_prob_at_mode = apply(x$k[,cell,drop=FALSE], 2, max)
    )
  attr(pred, "sigma2") <- x$sigma2
  attr(pred, "timesteps") <- x$timesteps
  pred
}



#' Retrieve trace for particular cell from list of mosaic fits
#'
#' @param x list of fits from the tj_fit_m0.5
#' @param cell which cell, in the original raster
#'
#' @export

tj_trace_m0.5 <- function(x, cell, add_sub_info = FALSE) {
  # find it and get it:
  them <-  which ( sapply(x, \(f) cell %in% f$cell_mapping$oldcell ) )

  tr <- lapply(them, \(si) {
    fi <- x[[si]]
    ncell <- fi$cell_mapping$cell[ match( cell, fi$cell_mapping$oldcell ) ]
    tc <- trace_m0.3(x[[si]], ncell)
    if(add_sub_info) tc <- tc |> mutate(subset=si, subset_cell = ncell, cell = cell)
    tc
  })
  if(length(them) == 1) tr[[them]] else tr
}


#' Retrieve trace for particular cell from list of mosaic fits
#'
#' @param x list of fits from the tj_fit_m0.3_dac
#' @param cell which cell, in the original raster
#'
#' @export

tj_trace_m0.3_dac <- function(x, cell, add_sub_info = FALSE) {
  # find it and get it:
  them <-  which ( sapply(x, \(f) cell %in% f$cell_mapping$oldcell ) )

  tr <- lapply(them, \(si) {
    fi <- x[[si]]
    ncell <- fi$cell_mapping$cell[ match( cell, fi$cell_mapping$oldcell ) ]
    tc <- trace_m0.3(x[[si]], ncell)
    if(add_sub_info) tc <- tc |> mutate(subset=si, subset_cell = ncell, cell = cell)
    tc
  })
  if(length(them) == 1) tr[[them]] else tr
}
