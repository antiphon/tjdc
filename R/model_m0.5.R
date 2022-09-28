#' Wrapper to fit model 3 to a mosaic
#'
#' Split datacube into sub-rasters, fit model 0.3's, and stitch together. Divide-and-conquer.
#'
#' @param x stars datacube
#' @param subsets st-polygons defining the subsets.
#' @param ... ignored
#'
#' @details
#' Wrap fitting of m0.3 in subsets for multicore computation.
#'
#'
#' @import looptimer doParallel foreach
#' @export

tj_fit_m0.3_dac <- function(x,
                            timevar = "z",
                            attrvar = "values",
                            dat = NULL, # alternative to x, timevar, attrvar
                            prior_theta  = list(m = c(a=0, b=0, d=0),
                                                S = diag(c(1, 1, 1)*1e5)),
                            prior_sigma2 = c(shape = 2, rate = 1), # inv-gamma
                            prior_k,
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
        cma <- di |> filter(x == x[1])
        #cells_to_ignore <- cma$cell[ match(.mcpars$cells_to_ignore, cma$oldcell) |> na.omit() ]
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
#' Average overlapping regions' pixel estimates. Version 2, use included cell mappings.
#'
#' @import dplyr
#' @export

tj_stitcher_m0.5_v1.2 <- function(list_of_fits, subsets) {
  stop("Not finished.")
  # figure out which have repeats
  cell_map <- lapply(seq_along(list_of_fits), \(si) list_of_fits[[si]]$cell_mapping |>
                       mutate(subset = si)) |>
    bind_rows()
  #
  cell_stats <- cell_map |> group_by(oldcell) |> count()
  allcells    <- unlist(cell_map$oldcell)
  ncells      <- max(allcells)
  cellcount   <- tabulate(allcells, ncells)
  est_cells  <- which(cellcount > 0) # which cells do we actually have in subsets
  are_repeated_txt      <- which(cellcount > 1) |> as.character()
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
  #
  posts$sigma2        <- sapply(list_of_fits, getElement, "sigma2") |> t() |> colMeans()
  posts$took          <- mean(lapply(list_of_fits, getElement, "took") |>
                         sapply(as.numeric, units = "secs")) |> as.difftime(units = "secs") # arbitrary
  posts$hist_sigma2   <- sapply(list_of_fits, getElement, "hist_sigma2")
  posts$cells         <- est_cells
  posts$cells_in_many <- are_repeated
  posts
}



###############################################
#' Stitch a mosaic model fit
#'
#' Average overlapping regions' pixel estimates.
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
  #
  posts$sigma2 <- sapply(list_of_fits, getElement, "sigma2") |> t() |> colMeans()
  posts$took   <- mean(lapply(list_of_fits, getElement, "took") |>
                         sapply(as.numeric, units = "secs")) |> as.difftime(units = "secs") # arbitrary
  posts$hist_sigma2 <- sapply(list_of_fits, getElement, "hist_sigma2")
  posts$cells  <- est_cells
  posts$cells_in_many <- are_repeated
  posts
}



################################################################
#' Divide bounding box into subsets, return polygons
#'
#' Divide both sides by a factor 2^q
#'
#' @import sf
#' @export
tj_divide_bbox <- function(bb,
                        power_of_2 = 2,
                        q = power_of_2,
                        expand = 0 ) {

  poly3 <- bb |> st_as_sfc()
  # Buffer?
  ex <- expand
  bb3 <- bb
  # shift
  poly30 <- poly3 - bb3[c(1,2)]
  factr  <- 1/(2^q/(1+2*ex))
  polytmpl <- poly30 * factr - c(bb3[3]-bb3[1], bb3[4]-bb3[2])*ex * factr # corner in
  # shifts
  n <- 2^q+1
  crns <- expand.grid(seq(bb3[1], bb3[3], l = n)[-n], seq(bb3[2], bb3[4], l = n)[-n])
  polys <- apply(crns, 1, \(s) polytmpl + s) |> do.call(what = c) |> st_as_sf( crs = st_crs(poly3)) |>
    mutate(id = sprintf("%04i", 1:n()) )
  polys
}



