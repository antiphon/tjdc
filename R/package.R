#' Trends and Jumps in Thin DataCubes
#'
#' Pixel level time-series analysis with the help of neighboring pixels.
#'
#' @details
#'
#' The `tjdc` package implements methods for fast trend and changepoint detection in stacked raster data, or data cubes.
#'
#' Current version is based on the `stars`-package data management, but this might change in the future.
#'
#' The main idea of the package is to provide user-friendly wrappers to split-calc-merge type of operations
#' on high spatial dimension, low temporal dimension data-cubes. Motivation and research focus has been
#' on datacubes for the Finnish National Forest Inventory with spatial resolution of 16m x 16m but with mere 10-15 time steps. Such
#' series are below the lower threshold for sophisticated change point and trend analysis techniques,
#' but assuming spatial correlation between nearby pixels we might be able to improve the power.
#'
#' The statistical functionality is roughly in two groups.
#'
#' @section Trend analysis:
#' This is functionality carried over from the `ConMK` package. The main ideas revolves around
#' the contextual Mann-Kendall's test [tj_contextual_mann_kendall()]. Nothing new in this package
#' compared to the `ConMK` other than switching the retired `raster` package to `stars`.
#'
#' @section Change-point analysis:
#' For short data sequences with high noise, detecting reliably even at most one change-point is challenging.
#' A simple linear trend model with a possible jump is estimated so that the jump events are correlated
#' between neighbours.
#'
#'
#'
"_PACKAGE"
