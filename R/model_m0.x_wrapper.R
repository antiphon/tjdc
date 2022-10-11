#' Wrapper to fit models of the 0.x-family to a mosaic
#'
#' Split datacube into sub-rasters, fit model 0.x's, and stitch together. Divide-and-conquer.
#'
#' @param x stars datacube
#' @param tiling given by [tj_divide_raster()]
#' @param cores number of cores to use in multicore computation using `doParallel`
#' @param cfg m0.x parameters, see [tj_cfg_m0.3(), tj_cfg_m0.6()]
#' @param ... parameters passed on to [tj_fit_m0.3(), tj_fit_m0.6()] (in case cfg not used)

#'
#' @details
#' Wrap fitting of m0.x in subsets for multicore computation. The tiling splits the stars-object
#' into subrasters, then [tj_fit_m0.3(), tj_fit_m0.6()] is called on each subraster.
#'
#' @return List of fits (keep_all_in_memory=TRUE) or list of files where the fits are stored (fitpath provided).
#' Note: The fits have cell_info element where the 'truecell' column provides the mapping to the original 'x'.
#'
#' @seealso [tj_cfg_m0.3(), tj_cfg_m0.6()]
#'
#' @import looptimer doParallel foreach stars dplyr
#' @export
tj_fit_m0.x <- function(x,
                        timevar     = "year",
                        attrvar     = "values",
                        tiling,
                        cfg,
                        ...,
                        ncores = 1 ,
                        verbose = TRUE, verbose2 = verbose, # extra
                        fitpath = NULL,
                        fitid = "youridhere",
                        keep_all_in_memory = TRUE, # keep all fits in memory? Only works under fitpath given
                        model_variant = "m0.3",
                        recalc = FALSE # incase find already fit pieces
) {
  t0 <- Sys.time()
  dbg <- verbose2
  if(!is(tiling, "tj_tiling")) stop("'tiling': Not from tj_divide_raster")
  if(!is(x, "stars")) stop("'x': Only stars supported.")

  if(is.null(fitpath) & !keep_all_in_memory) stop("mixed signals: no fitpath but also dont want to keep in memory.")

  if(!model_variant %in% c("m0.3", "m0.6")) stop("model_variant needs to be on of: 'm0.3', 'm0.6'.")

  if(missing(cfg)) cfg <- get(sprintf("tj_cfg_%s", model_variant))(...)

  tj_fit_fun <- get(sprintf("tj_fit_%ss", model_variant))

  ### Main call function for each tile.
  fit_one <- function(i) {
    # subdivide by row-column
    rc <- tiling$rc[i, ] |> unlist()
    ri <- x[, rc['cmin']:rc['cmax'], rc['rmin']:rc['rmax']]
    # go
    res <- tj_fit_fun( x= ri, timevar = timevar, attrvar = attrvar, cfg = cfg)
    # add cell mapping
    res$cell_info |> mutate(truecell = tiling$cell[[i]])
    res
  }

  # prepapre wrapping
  ntiles       <- prod(tiling$n)
  tiles_needed <- 1:ntiles
  # File management:
  if(!is.null(fitpath)) {
    # filename template for this run
    froot          <- sprintf("%s/tjdc_%s_%s_%%04i.rds", model_variant, fitpath, fitid)
    fit_file_names <- sprintf(froot, tiles_needed)
    fit_not_found  <- !file.exists(fit_file_names)
    if( any(!fit_not_found) & !recalc ) {
      tiles_needed <- tiles_needed[ fit_not_found ]
      message(sprintf("Found existing fits for tiles: %s. Will not recalculate (recalc=FALSE).",
                      paste0(which(!fit_not_found), collapse = ",")))
    }
  }

  ## Start cluster
  if(ncores>1) doParallel::registerDoParallel(ncores)
  # Cluster needs
  exports <- c(
    "x",
    "cfg",
    "tiling",
    "tj_fit_fun",
    "mutate")


  ## Loop over
  blocks <- looptimer::split_block(tiles_needed, ncores)

  ## Start clock
  timer0 <- looptimer( n = length(blocks) ,
                       fg = 2,
                       endline = "     \r",
                       prefix = sprintf("[m0.x %ic]", ncores))
  r5 <- list() # results
  ##### Main looper ######
  for(bl in blocks) {

    ### fit them
    if(ncores>1) {
      rx <- foreach(x = bl, .export = exports) %dopar% fit_one(x)
    } else rx <- lapply(bl, fit_one)
    # keep
    r5[bl] <- rx


    ## Store fits? Drop afterwards?
    if(!is.null(fitpath)) {
      for(bli in bl) {
        saveRDS(r5[[bli]], fit_file_names[bli])
        if(dbg) message(sprintf("\nWrote %s.", fit_file_names[bli]), appendLF = TRUE)
      }
      if(!keep_all_in_memory) r5[bl] <- NULL
    }
    ##
    ## Done main loop
    if(verbose) print(timer0 <- looptimer(timer0))
  }#####

  # Clear up
  if(ncores>1)doParallel::stopImplicitCluster()
  if(verbose) message(summary(timer0))

  if(!keep_all_in_memory) return(fit_file_names)
  # if we actually want to keep all in memory, load missing
  if(!is.null(fitpath)) {
    for(bli in which(!fit_not_found) ) r5[[bli]] <- readRDS(fit_file_names[bli])
    message(sprintf("[keep_all_in_memory=TRUE] Restored %i fits.", sum(!fit_not_found)))
  }
  # done and done.
  return(r5)
}


