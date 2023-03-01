#' Wrapper to fit model 3 to a mosaic
#'
#' (OBSOLETE) Split datacube into sub-rasters, fit model 0.3's, and stitch together. Divide-and-conquer.
#'
#' @param x stars datacube
#' @param tiling given by [tj_divide_raster()]
#' @param cores number of cores to use in multicore computation using `doParallel`
#' @param cfg m0.3 parameters, see [tj_cfg_m0.3()]
#' @param ... parameters passed on to [tj_fit_m0.3()] (in case cfg not used)

#'
#' @details
#' Wrap fitting of m0.3 in subsets for multicore computation. The tiling splits the stars-object
#' into subrasters, then [tj_fit_m0.3()] is called on each subraster.
#'
#' This has been replaced by the more general tj_fit_m0.x.
#'
#' @return List of fits (keep_all_in_memory=TRUE) or list of files where the fits are stored (fitpath provided).
#' Note: The fits have cell_info element where the 'truecell' column provides the mapping to the original 'x'.
#'
#'
#'
#' @seealso [tj_cfg_m0.3()] [tj_fit_m0.x]
#'
#' @import looptimer doParallel foreach stars dplyr
tj_fit_m0.5 <- function(x,
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
                        recalc = FALSE # incase find already fit pieces
) {
  t0 <- Sys.time()
  dbg <- verbose2
  if(!is(tiling, "tj_tiling")) stop("'tiling': Not from tj_divide_raster")
  if(!is(x, "stars")) stop("'x': Only stars supported.")

  if(is.null(fitpath) & !keep_all_in_memory) stop("mixed signals: no fitpath but also dont want to keep in memory.")

  if(missing(cfg)) cfg <- tj_cfg_m0.3(...)

  ### Main call function for each tile.
  fit_one <- function(i) {
    # subdivide by row-column
    rc <- tiling$rc[i, ] |> unlist()
    ri <- x[, rc['cmin']:rc['cmax'], rc['rmin']:rc['rmax']]
    # go
    res <- tj_fit_m0.3s( x= ri, timevar = timevar, attrvar = attrvar, cfg = cfg)
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
    froot          <- sprintf("%s/tjdc_%s_%%04i.rds", fitpath, fitid)
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
    "tj_fit_m0.3s", "mutate")


  ## Loop over
  blocks <- looptimer::split_block(tiles_needed, ncores)

  ## Start clock
  timer0 <- looptimer( n = length(blocks) ,
                       fg = 2,
                       endline = "     \r",
                       prefix = sprintf("[m0.5 %ic]", ncores))
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





#' Wrapper to fit model 3 to a mosaic, old version
#'
#' Split datacube into sub-rasters, fit model 0.3's, and stitch together. Divide-and-conquer.
#'
#' @param x stars datacube
#' @param subsets list of cell indices
#' @param ... ignored
#'
#' @details
#' Wrap fitting of m0.3 in subsets for multicore computation.
#'
#' Obsolete. Use tj_fit_m0.x.
#'
#'
#' @import looptimer doParallel foreach

tj_fit_m0.3_dac <- function(x,
                            timevar = "z",
                            attrvar = "values",
                            dat = NULL, # alternative to x, timevar, attrvar
                            prior_theta  = list(m = c(a=0, b=0, d=0),
                                                S = diag(c(1, 1, 1)*1e5)),
                            prior_sigma2 = c(shape = 2, rate = 1), # inv-gamma
                            prior_k = 0.5,
                            gamma = 0,
                            verbose = TRUE, dbg=FALSE,
                            method = "mc",
                            niter = 100,
                            subsets = NULL,
                            cells_to_ignore = NULL,
                            truncate_jump_at_mean = 0,
                            ...,
                            ctrl = list(eps = 0.1, burnin = .5, thin_steps = 1),
                            stitcher = tj_stitcher_m0.5_v1.1,
                            ncores = 1 ,
                            fitpath = NULL,
                            fitid = "youridhere",
                            just_fit = !is.null(fitpath), # do not parse? Only makes sense with
                            keep_all_in_memory = TRUE, # keep all fits in memory? Only works under fitpath given
                            recalc = FALSE, # incase find already fit pieces
                            keep_hist = FALSE # keep MCMC history? Will increase file size (niter+1) times.

) {
  #
  if(is.null(subsets) || !is.list(subsets)) stop (" 'subsets', a list of cells per subset, not provided.")
  #
  if(!keep_all_in_memory & is.null(fitpath) & !just_fit) stop("Parameters make no sense (keep in memory & only fit).")

  if(!is.null(fitpath)) if(!dir.exists(fitpath)) stop(sprintf("path %s does not exist.", fitpath))

  #
  if(is.null(dat)) dat <- tj_stars_to_data(x, timevar, attrvar)
  # check that we have expected names in the dat
  gdims <- attr(dat, "grid")
  if( !all(c("cell", "time", "value") %in% names(dat))  | is.null(gdims))
    stop("'dat' ill formatted. need at least cell, time, value-columns, and `gdims` attribute. See `tj_stars_to_data`. ")
  #
  t0 <- Sys.time()
  cells <- 1:prod(gdims)
  nsubsets <- length(subsets)
  #
  rnc <- gdims['ncol']
  rnr <- gdims['nrow']

  # Ball up parameters
  .mcpars <<- list(prior_sigma2 = prior_sigma2,
                   prior_k = prior_k,
                   prior_theta = prior_theta,
                   gamma = gamma,
                   niter = niter,
                   rnr = rnr, ctrl = ctrl,
                   rnc = rnc,
                   cells_to_ignore = cells_to_ignore,
                   keep_hist = keep_hist,
                   truncate_jump_at_mean = truncate_jump_at_mean)

  # Choose model to fit:
  .fitter <<- tj_fit_m0.3
  .subsets <<- subsets
  #
  ###############################
  # Fitter function, call for each sub-data
  #
  #
  fun1 <- function(subset_i) {
        # chip off the data and re-parameterise grid
        j <- subsets[[subset_i]]
        di <- .datx |>
          subset(cell %in% j) |>
          mutate(
            oldcell  = cell,
            row   = factor(row) |> as.integer(),
            col   = factor(col) |> as.integer(),
            cell  = col + (row-1) * max(col))
        #
        # solve the cells to ignore parameter
        # make sure to drop the dimensions down to save memory
        attr(di, "gdims") <- c(nrow = max(di$row), ncol = max(di$col))
        #
        #browser()
        fit1 <- .fitter(dat         = di,
                        niter       = .mcpars$niter,
                        prior_theta = .mcpars$prior_theta,
                        prior_k     = .mcpars$prior_k,
                        gamma       = .mcpars$gamma,
                        prior_sigma2 = .mcpars$prior_sigma2,
                        ignore_cells_with_na  =  TRUE,
                        verbose     = FALSE,
                        keep_hist   = .mcpars$keep_hist,
                        truncate_jump_at_mean = .mcpars$truncate_jump_at_mean,
                        ctrl = .mcpars$ctrl)
        # Add old cell index
        fit1$cell_mapping <- di |> filter(time == time[1]) |> select(cell, oldcell)
        # for now
        fit1
  }
  # go in bloks  for parallel
  subsets_needed <- 1:nsubsets
  # Check if intermediate files exist
  if(!is.null(fitpath)) {
    # filename template for this run
    froot          <- sprintf("%s/tjdc_%s_%%04i.rds", fitpath, fitid)
    fit_file_names <- sprintf(froot, subsets_needed)
    fit_not_found  <- !file.exists(fit_file_names)
    if( any(!fit_not_found) & !recalc ) {
      subsets_needed <- subsets_needed[ fit_not_found ]
      message(sprintf("Found fits for subsets: %s. Will not recalculate (recalc=FALSE).",
                      paste0(which(!fit_not_found), collapse = ",")))
    }

  }

  bloks <- looptimer::split_block(subsets_needed, ncores)
  if(ncores>1) doParallel::registerDoParallel(ncores)


  timer0 <- looptimer( n = length(bloks) , fg = 3, endline = "     \r", prefix = sprintf("[m0.4 %s %ic]", method, ncores))

  r5 <- list() # results
  # Main looper

  for(bl in bloks) {
    # this data suffices to be passed on to workers
    jj <- subsets[bl] |>
      unlist() |>
  #    c(bl) |> legacy from m0.4
      unique()
    .datx <<- subset(dat, cell %in% jj)

    exports <- c(#"tj_cell_neighbours",
                 ".datx", ".mcpars",
                 ".subsets",
                 "transmute",
                 ".fitter", "looptimer")

    if(ncores>1) {
      rx <- foreach(x = bl, .export = exports) %dopar% fun1(x)
    } else rx <- lapply(bl, fun1)
    if(verbose) print(timer0 <- looptimer(timer0))
    r5[bl] <- rx

    #browser()
    # Store fits? Drop afterwards?
    if(!is.null(fitpath)) {
      for(bli in bl) {
        if(dbg) message(sprintf("Writing %s...", fit_file_names[bli]), appendLF = FALSE)
        saveRDS(r5[[bli]], fit_file_names[bli])
        if(dbg) message(" ok.")
      }
      if(!keep_all_in_memory) r5[bl] <- NULL
    }
    #browser()
  }
  if(ncores>1)doParallel::stopImplicitCluster()
  if(verbose) message(summary(timer0))

  # In case we don't stitch
  if(just_fit) {
    if(!is.null(fitpath)) return(fit_file_names)
    else return(r5)
  }  # done
  # Load missing
  if(!is.null(fitpath)) {
    for(bli in which(!fit_not_found) ) r5[[bli]] <- readRDS(fit_file_names[bli])
    message(sprintf("Restored %i fits.", sum(!fit_not_found)))
  }
  # This is the magical part.
  # gather
  res <- stitcher(r5, subsets)
  #
  res$took <-  Sys.time() - t0
  #
  res
}



###############################################
#' Stitch a mosaic model fit
#'
#' Average overlapping regions' pixel estimates. Old format.
#'
#'
#' @export

tj_stitcher_m0.5_v1.1 <- function(list_of_fits, subsets) {
  # figure out which have repeats
  mess   <- unlist(subsets)
  ncells <- max(mess)
  x      <- tabulate(mess, ncells)
  est_cells  <- which(x > 0) # which cells do we actually have in subsets
  are_repeated_txt      <- which(x > 1) |> as.character()
  are_repeated          <- are_repeated_txt  |> as.integer()
  are_repeated_idx      <- match( are_repeated, est_cells )
  not_repeated_sets     <- lapply(subsets, setdiff, are_repeated)
  not_repeated_sets_idx <- lapply(not_repeated_sets, match, est_cells)
  #
  # Do global variables separate, here the pixelwises
  vars <- c("k", "z", "theta_m")
  posts <- list()
  #browser()
  for(var in vars) {
    # Go first for the non-repeated
    V1 <- lapply(seq_along(subsets), \(si) {
      ss_idx          <- match( not_repeated_sets[[si]], subsets[[si]]  )
      o               <- rbind(list_of_fits[[si]][[var]])[,ss_idx]
      attr(o, "cidx") <- not_repeated_sets[[si]]
      o
    })
    # then for the repeated:
    Vx <- lapply(seq_along(subsets), \(si) {
      idx    <- which(are_repeated %in% subsets[[si]])
      ss_idx <- match(are_repeated[idx], subsets[[si]]  )
      o      <- rbind(list_of_fits[[si]][[var]])[,ss_idx, drop=FALSE]
      attr(o, "cidx") <- are_repeated[idx]
      o
    })
    # gather
    V <- matrix(NA, nrow = nrow( rbind( list_of_fits[[1]][[var]] ) ), ncol = ncells)
    for(si in seq_along(subsets)) {
      V[, attr( V1[[si]], "cidx" )] <- V1[[si]]
      # And then: pointwise means
      if(length(are_repeated)) { ###
        r_c <- attr( Vx[[si]], "cidx" )
        p1  <- V[, r_c] # old
        p2  <- Vx[[si]]
        # for averaging
        nc <- x[r_c]
        p2 <- t( t(p2)/nc )
        #
        p12 <- list(p1, p2) |> simplify2array()
        ndi <- length(dim(p12))

        # Average, but only those that are not completely NA
        m_m <- apply(p12, 1:(ndi-1), \(v) if(all(is.na(v))) NA else sum(v, na.rm=TRUE))
        # store
        V[, attr( Vx[[si]], "cidx" )] <- m_m
      }
    }
    # and mean

    posts[[var]] <- V
  }
  # Gather cell information, remapping to original cell indices.
  cinfo <- lapply(seq_along(list_of_fits), \(i) {
    f <- list_of_fits[[i]]
    f$cell_info |> left_join(f$cell_mapping, by = "cell") |> mutate(subset_cell = cell,
                                                                    cell = oldcell,
                                                                    subset = i,
                                                                    oldcell = NULL)
  }) |>
    bind_rows()
  #
  posts$sigma2 <- sapply(list_of_fits, getElement, "sigma2") |> t() |> colMeans()
  posts$took   <- mean(lapply(list_of_fits, getElement, "took") |>
                         sapply(as.numeric, units = "secs")) |> as.difftime(units = "secs") # arbitrary
  posts$hist_sigma2 <- sapply(list_of_fits, getElement, "hist_sigma2")
  posts$cells  <- est_cells
  posts$cells_in_many <- are_repeated
  posts$cell_info <- cinfo
  posts$timesteps <- list_of_fits[[1]]$timesteps

  posts
}



