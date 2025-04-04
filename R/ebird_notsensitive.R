#' @export
get_species <- function(x){
  stopifnot(is.character(x))
  r <- ebirdst::ebirdst_runs
  x <- tolower(trimws(x))
  code <- match(x, tolower(r$species_code))
  sci <- match(x, tolower(r$scientific_name))
  com <- match(x, tolower(r$common_name))
  if(is.na(code)){
    return(r$species_code[dplyr::coalesce(code, sci, com)]  )
  } else{
    return(r$common_name[dplyr::coalesce(code, sci, com)])
  }
}

create_basemap <- function(basemap, bbox = NULL, zeros = NULL) {
  stopifnot(all(c("land", "lakes") %in% names(basemap)))
  if (is.null(bbox)) {
    bbox <- sf::st_bbox(basemap$land)
  } else if (is.numeric(bbox)) {
    bbox <- sf::st_bbox(bbox, crs = sf::st_crs(basemap$land))
  } else if (!inherits(bbox, "bbox")) {
    stop("Invalid bbox object")
  }
  if (!is.null(zeros)) {
    if (inherits(zeros, c("sf", "sfc", "SpatRaster"))) {
      zeros <- list(col_zero = zeros)
    } else {
      stopifnot(length(zeros) == 2, names(zeros) == c("assumed", "predicted"))
      names(zeros) <- c("col_d_zero_assumed", "col_d_zero_pred")
    }
    for (i in seq_along(zeros)) {
      if (is.null(zeros[[i]])) {
        zeros[[i]] <- NULL
        next
      }
      stopifnot(inherits(zeros[[i]], c("sf", "sfc", "SpatRaster")))
    }
  }

  mt <- ebirdstwf::map_theme

  # square bounding box
  e <- square_extent(bbox)
  bbox <- extent_to_sf(e, crs = sf::st_crs(basemap$land))
  plot(bbox, col = NA, border = NA, axes = FALSE, bty = "n", reset = FALSE)

  # land
  plot(sf::st_geometry(basemap$land),
       col = mt$col_land, border = NA,
       axes = FALSE, bty = "n", reset = FALSE, add = TRUE)

  # zeros
  for (nm in names(zeros)) {
    if (inherits(zeros[[nm]], c("sf", "sfc"))) {
      plot(sf::st_geometry(zeros[[nm]]),
           col = mt[[nm]], border = NA,
           axes = FALSE, bty = "n", reset = FALSE, add = TRUE)
    } else if (inherits(zeros[[nm]], "SpatRaster")) {
      terra::plot(zeros[[nm]], col = mt[[nm]],
                  maxcell = terra::ncell(zeros[[nm]]),
                  legend = FALSE, axes = FALSE, bty = "n",
                  reset = FALSE, add = TRUE)
    }
  }

  # lakes
  plot(sf::st_geometry(basemap$lakes),
       col = mt$col_lakes, border = mt$col_land_border,
       lwd = mt$lwd_land_border,
       axes = FALSE, bty = "n", reset = FALSE, add = TRUE)

  # land border
  plot(sf::st_geometry(basemap$land),
       col = NA, border = mt$col_land_border,
       lwd = mt$lwd_land_border,
       axes = FALSE, bty = "n", reset = FALSE, add = TRUE)

  invisible()
}

#' Extract `params[[@param_name]]` from multiple interactive stixels, and union them into a complete set
#' @export
union_stixel_params <- function(dir_interactive_stx, stxl_ids=NULL, param_name){
  files_available <- grep('fit-predict-stixel',
                          list.files(dir_interactive_stx, full.names = TRUE),
                          value=TRUE)
  if(!is.null(stxl_ids)){
    interactive_stixel_files <- grep(pattern = paste(stxl_ids, collapse="|"),
                                     x = files_available,
                                     value = TRUE)
  } else{
    # Load all stixels in the directory
    interactive_stixel_files <- files_available
  }

  if(param_name == "srd"){
    return(do.call(dplyr::bind_rows, lapply(interactive_stixel_files, function(interactive_stixel){
      message("Reading: ", basename(interactive_stixel), " for stixel params")
      return(readRDS(interactive_stixel)$srd)
    })))
  }
  unique(unlist(lapply(interactive_stixel_files, function(interactive_stixel){
    message("Reading: ", basename(interactive_stixel), " for stixel params")
    ps <- readRDS(interactive_stixel)$params
    return(ps[[param_name]])
  })))
}


#' Run a quantile random forest and a regression random forest for each species in stx_dat.
#' Then compute out-of-bag predictions back to the data using both models. In the
#' case of the regression RF, compute ppois(obs, pred) to get P(count <= observed count).
#' @param dat data representing one stixel
#' @param dir_interactive_stx directory from which to retrieve stixel parameters,
#' specifically `PREDICTOR_LIST`, `ALWAYS_SPLIT_COVARS`, and `FACTOR_COVARS`
#' @param model_outdir where to write models (for all species, for all stixel IDs)
#' @export
compute_species_rfs <- function(stx_dat, spp_colnames, #dir_interactive_stx,
                               predictors,
                               always_splits,
                               factors,
                               model_outdir,
                               overwrite){
  stopifnot(dir.exists(file.path(model_outdir, 'data')))
  always_splits <- intersect(always_splits, predictors)
  factors <- intersect(factors, predictors)
  stopifnot(dir.exists(model_outdir))
  stopifnot(all(predictors %in% names(stx_dat)))
  stopifnot(all(factors %in% names(stx_dat)))
  stopifnot(all(always_splits %in% names(stx_dat)))

  stxl_id <- unique(stx_dat$stixel_id)
  stopifnot(length(stxl_id) == 1)


  # Change counts slashes to zeroes:
  sp <- stx_dat %>% dplyr::select(tidyselect::all_of(spp_colnames))
  sp[sp == -2] <- 0

  # Remove species without sufficient observations in this stixel:
  to_remove <- sapply(sp, function(x){
    # at least 10 pos obs with counts over all time
    length(x[x > 0]) <= 10
  })
  to_remove <- to_remove[to_remove] %>% names()
  stx_dat <- cbind(stx_dat %>% dplyr::select(-tidyselect::all_of(spp_colnames)),
                   sp %>% dplyr::select(-tidyselect::all_of(to_remove)))
  rm(sp)

  ####
  # Modify data classes
  ####
  classes <- sapply(stx_dat %>%
                    dplyr::select(-tidyselect::any_of(spp_colnames)),
                    class)

  for(varname in names(which(classes == "logical"))) {
    stx_dat[[varname]] <- as.integer(stx_dat[[varname]])
  }

  ####
  # Ensure proper formatting of factors
  ####
  for(i in 1:length(factors)){
    stx_dat[[factors[i]]] <- as.factor(stx_dat[[factors[i]]])
  }

  ####
  # Species-independent runs:
  ####
  for(species_name in setdiff(spp_colnames, to_remove)){
    message("\n**", species_name, "**\n")
    outfile <- file.path(model_outdir, sprintf('data/%s_datpred_%s.rds', species_name, stxl_id))
    if(file.exists(outfile) & overwrite == FALSE){
      message("\nFile exists.\nSkipping ", outfile)
      next()
    } else{
      ####
      # Separate data containing real counts and data containing X's only:
      ####
      sp_dat_counts <- stx_dat[stx_dat[[species_name]] > 0,] %>%
        dplyr::select(tidyselect::all_of(
          c('checklist_id', 'observer_id', "year", species_name, predictors)))

      sp_dat_x <- stx_dat[stx_dat[[species_name]] == -1, ] %>%
        dplyr::select(tidyselect::all_of(
          c('checklist_id', 'observer_id', "year", species_name, predictors)))


      ####
      # Create formula:
      ####
      form <- paste0(sprintf('%s ~ ', species_name),
                     paste0(predictors,
                            collapse=' + '))


      ####
      # Compute count models:
      ####
      message("Fitting quantile RF model")
      modq <- ebirdstwf:::run_rf(f=form,
                                 dat=sp_dat_counts,
                                 filename_rf=file.path(model_outdir,
                                                       sprintf('%s_rfquantile_%s.rds', species_name, stxl_id)),
                                 quantreg=TRUE,
                                 # need this for OOB quantile prediction:
                                 keep.inbag=TRUE,
                                 importance='impurity',
                                 respect.unordered.factors = "order",
                                 always.split.variables = always_splits)
      quantile_predictions <- predict(modq, type='quantiles', quantiles = c(0.25, 0.50, 0.75))$predictions
      colnames(quantile_predictions) <- sprintf('%s_%s',
                                                gsub("quantile= ", 'pred',
                                                     colnames(quantile_predictions), fixed=TRUE),
                                                species_name)
      sp_dat_counts <- cbind(sp_dat_counts, as.data.frame(quantile_predictions))

      message("Fitting regression RF model")
      mod <- ebirdstwf:::run_rf(f=form,
                                # filter to pos obs for this particular species:
                                # Note that Xs are coded as 1s
                                dat=sp_dat_counts,
                                filename_rf=file.path(model_outdir, sprintf('%s_rf_%s.rds', species_name, stxl_id)),
                                importance='impurity',
                                respect.unordered.factors = "order",
                                always.split.variables = always_splits)

      # Regression RF predictions:
      sp_dat_counts[[sprintf('pred_%s', species_name)]] <- mod$predictions

      # Compute p(X < x) using regression RF predictions
      sp_dat_counts[[sprintf('obs_%s', species_name)]] <- sp_dat_counts[[species_name]]
      sp_dat_counts[[sprintf('ppois_%s', species_name)]] <- ppois(q=sp_dat_counts[[sprintf('obs_%s', species_name)]],
                                                                  lambda=sp_dat_counts[[sprintf('pred_%s', species_name)]])
      ####
      # Use both quantile and regression models to make predictions to X's
      ####
      if(nrow(sp_dat_x) > 0) {
        sp_dat_x[[sprintf('obs_%s', species_name)]] <- sp_dat_x[[species_name]] # all -1
        # Regression RF predictions (out-of-sample):
        sp_dat_x[[sprintf('pred_%s', species_name)]] <- predict(mod, sp_dat_x)$predictions
        # Quantile RF predictions:
        qps <- predict(modq, sp_dat_x, type='quantiles', quantiles = c(0.25, 0.50, 0.75))$predictions
        colnames(qps) <- sprintf('%s_%s',
                                 gsub("quantile= ", 'pred', colnames(qps), fixed=TRUE),
                                 species_name)
        sp_dat_x <- bind_cols(sp_dat_x, as.data.frame(qps))

        # To get a p-value (from regression RF prediction), treat X's as 1's
        sp_dat_x[[sprintf('ppois_%s', species_name)]] <- ppois(q=1,
                                                               lambda=sp_dat_x[[sprintf('pred_%s', species_name)]])


      }
      # bind the X's with the positive counts:
      sp_dat <- bind_rows(sp_dat_counts, sp_dat_x) %>%
        select(any_of(c('checklist_id', 'observer_id', sprintf('obs_%s', species_name),
                        sprintf('ppois_%s', species_name))),
               contains('pred'))
      saveRDS(sp_dat, outfile)


    }

}



}
#' Compute some function on each stixel in the data. If `stixel_ids`==NULL, will operate on each unique `stixel_id` present in the data
#' Returns function results in a list
# stixelwise <- function(dat, fun, stixel_ids=NULL, ...){
#   stopifnot('stixel_id' %in% names(dat))
#   if(is.null(stixel_ids)){
#     stixel_ids <- unique(dat$stixel_id)
#   }
#   lapply(stixel_ids, function(stxl_id){
#     stxl_dat <- dat %>%
#       dplyr::filter(stixel_id == stxl_id)
#
#   })
# }


#' tixelwise
#'
#' @param fun function to apply to `dat` tixelwise. Must accept data as its first positional argument, and 'stixel_id' as a keyword argument.
#' @param dat data frame, which is passed as the first function to `fun`. required columns:
#' @param pass_tixel_id boolean; pass tixel number as an arg to `fun`?
#'
#' @export
tixelwise <- function(dat1, dat1_name, dat2=NULL, dat2_name=NULL, fun, tixels=NULL, pass_tixel_id_as=NULL, args=NULL){
  stopifnot('tixel' %in% names(dat1))
  if(is.null(tixels)){
    tixels <- unique(dat1$tixel)
  }

  l <- lapply(tixels, function(tix){
    message('tixel: ', tix)
    tdat1 <- dat1 %>% dplyr::filter(tixel==tix)
    argz <- c(list(`dat`=tdat1), args)
    names(argz)[[1]] <- eval(dat1_name)
    if(!is.null(pass_tixel_id_as)){
      argz <- c(argz,  list(`tixel`=tix))
      names(argz)[[length(argz)]] <- pass_tixel_id_as
    }
    if(!is.null(dat2)){
      stop('TODO: BUILD THIS OUT')
      stopifnot('tixel' %in% names(dat2))
      stopifnot( all(tixels %in% unique(dat2$tixel)) )
      tdat2 <- dat2 %>% dplyr::filter(tixel==tix)
      # do.call(fun(tdat1, tdat2, c(pass_tixel_id_as, ...))
    } else{
      do.call(fun, argz)
    }
  })
  names(l) <- sprintf('tixel_%s', tixels)

  # If all elements are data.frames, return one single data frame by row-binding:
  class_l <- sapply(l, function(x){
    "data.frame" %in% class(x)
  })

  if(all(class_l)){
    return(dplyr::bind_rows(l))
  }

  return(l)
}


#' @export
make_tixels <- function(doy, start_doy){
  if (start_doy < 1 || start_doy > 366) {
    stop("Invalid start day of year. Must be in the range [1, 366].")
  }

  groups <- list()

  for (i in 1:13) {
    group_start <- ((i - 1) * 28) %% 366 + start_doy
    group_end <- group_start + 27

    groups[[i]] <- cycle_date(group_start:group_end)
  }
  names(groups) <- c(1:13)
  df <- do.call(rbind,
                lapply(names(groups), function(name){
                  data.frame(day=groups[[name]], bin=name)
                }))
  leaveouts <- doy[!doy %in% df$day] %>% unique()
  df <- bind_rows(df,
                  do.call(rbind, lapply(leaveouts, function(lo){
                    # find day that is closest and use whatever that group is
                    return(data.frame(day=lo, bin=df$bin[which.min(abs(df$day - lo))]))
                  }))
  )
  bins <- sapply(doy, function(d){
    df$bin[match(d, df$day)]
  })
  stopifnot(!any(is.na(bins)))
  return(bins)
}

#' Map habitats
#'
#' Cluster habitats using kmeans, map them onto a colored plot and describe clusters
#' @export
map_habitats <- function(path_erd,
                         n_clusters=8,
                         lat_min=44.13,
                         lat_max=44.46,
                         lon_min=-74.22,
                         lon_max=-73.33){

  cluster_vars <- cmc.utils::names2chr(mcd12q1_lccs1_c1_pland,
                                       mcd12q1_lccs1_c11_pland,
                                       mcd12q1_lccs1_c13_pland,
                                       mcd12q1_lccs1_c14_pland,
                                       mcd12q1_lccs1_c15_pland,
                                       mcd12q1_lccs1_c21_pland,
                                       mcd12q1_lccs1_c22_pland,
                                       mcd12q1_lccs1_c31_pland,
                                       mcd12q1_lccs1_c32_pland,
                                       mcd12q1_lccs1_c42_pland,
                                       mcd12q1_lccs1_c43_pland,
                                       mcd12q1_lccs2_c25_pland,
                                       mcd12q1_lccs2_c35_pland,
                                       mcd12q1_lccs2_c36_pland,
                                       mcd12q1_lccs3_c27_pland,
                                       mcd12q1_lccs3_c50_pland,
                                       astwbd_c2_pland, astwbd_c3_pland,
                                       gsw_c2_pland, gsw_c3_pland,
                                       elevation_250m_median, elevation_250m_sd,
                                       elevation_30m_median, elevation_30m_sd)

  erd_checklists <- arrow::read_parquet(path_erd) %>%
    dplyr::filter(dplyr::between(latitude, lat_min, lat_max),
                  dplyr::between(longitude, lon_min, lon_max)) %>%
    select(all_of(c('latitude', 'longitude', cluster_vars)))



  n_clusterings = length(n_clusters)

  kclust <- kmeans(x=erd_checklists %>% select(all_of(cluster_vars)), centers=6)
  erd_checklists$hab_cluster <- kclust$cluster

  # sa_map <- OpenStreetMap::openmap(c(lat2, lon1), c(lat1, lon2), zoom = 10,
  #                   type = "esri-topo", mergeTiles = TRUE)
  #
  # # reproject onto WGS84
  # sa_map2 <- openproj(sa_map)
cl_sf <- erd_checklists %>% cmc.utils::dat2sf() %>% sf::st_transform(crs=3857)
basemap_ggplot(ext=sf::st_bbox(cl_sf),
               map_service='osm',
               map_type='streets')
 p <- erd_checklists %>%
    cmc.utils::dat2sf() %>%
    ggplot() +
    geom_sf(aes(color=factor(hab_cluster)))

descr <- as.data.frame(kclust$centers)
# find which ones have no variance
not_meaningful <- which(sapply(names(descr), function(colname){
  var(descr[[colname]]) == 0
})) %>% names()

descr <- descr %>% select(-all_of(not_meaningful))
ebirdst_varnames(names(descr))

kclust$iter

kclust$ifault


descr %>%
  mutate(cluster=1:6) %>%
  tidyr::pivot_longer(cols=-c('cluster', 'latitude', 'longitude'))
  ggplot(aes(x))

}


#' Run trends propensity score model
#'
#' @param model_type one of "classification", "grf"
#' returns predictions and a spatial plot if
#' @export
run_ps_model <- function(dat, model_type, response, features, importance, probability=NULL){
  stopifnot(importance %in% c('none', 'impurity'))
  stopifnot(!(model_type == 'grf' & !is.null(probability)))
  stopifnot(!(model_type == 'classification' & is.null(probability)))
  missings <- setdiff(c(response, features), names(dat))
  if(length(missings) > 0){
    stop("Missing from data: ", paste(missings, collapse=", "))
  }
  stopifnot(all(c(response, features) %in% names(dat)))
  n_trees_grf <- 2000
  W <- dat[[response]]

  if(model_type=='grf'){
    w_forest <- grf::regression_forest(dat[, features],
                                           W,
                                           num.trees = n_trees_grf,
                                           tune.parameters = "all")
    W_hat <- stats::predict(w_forest)$predictions
  } else if(model_type=='classification'){
    w_forest <- ranger::ranger(data=dat,
                                   formula=paste0(sprintf("%s ~", response),
                                                  paste(features, collapse=' + ')),
                                   num.trees=n_trees_grf,
                                   probability=probability,
                               importance=importance)
    W_hat <- cmc.utils::preds_from_modal_probs(w_forest$predictions, rf=w_forest)
  }else{
    stop("unknown propensity model type", model_type)
  }
  if(all(c('latitude', 'longitude') %in% names(dat))){
    p_spatial <- dat %>%
      select(latitude, longitude) %>%
      bind_cols(data.frame(W=W, What=W_hat)) %>%
      dat2sf() %>%
      ggplot2::ggplot() +
      ggplot2::geom_sf(alpha=0.5, ggplot2::aes(color=factor(W))) +
      ggplot2::facet_wrap(~What) +
        ggplot2::labs(title="Propensity",
             color=sprintf('true %s', response))
  } else{
    p_spatial <- NULL
  }

  p_calibration <- data.frame(W=W, What=W_hat) %>%
    ggplot2::ggplot(ggplot2::aes(x=W, y=What)) +
    ggplot2::geom_abline(slope=1, color='red', linetype='dashed') +
    ggplot2::geom_hex() +
    ggplot2::scale_x_continuous(breaks=2012:2022) +
    ggplot2::scale_y_continuous(breaks=2012:2022)


  return(list(What=W_hat,
              p_spat = p_spatial,
              p_calib = p_calibration,
              model=w_forest))
}

##############################################################
### Handle tweaks inside a function called prepare_features ##
##############################################################
#' @export
prepare_features <- function(config_path, aggregated, include_atlas, exclude_moon_etc,
                             exclude_landcover){
  p <- arrow::read_json_arrow(config_path)

  propensity_features <- p$propensity_features %>% unlist()
  marginal_features <- p$marginal_features %>% unlist()
  # variables that will be fixed to a set quantile:
  eff_covs <- setdiff(c((p$effort_features %>% unlist()),
                        (p$aggregated_effort_features %>% unlist())),
                      "day_of_year")

  # Propensity features refinements:
  if(include_atlas){
    propensity_features <-  c(propensity_features,'atlas')
    # variables that will be fixed to a set quantile:
    eff_covs <- c(eff_covs, 'atlas')
  }
  if(exclude_moon_etc){
    propensity_features <- setdiff(
      # replace day_of_year with is_weekend:
      c(propensity_features, 'is_weekend'),
      c(grep('cds_*', propensity_features, value = TRUE),
        'moon_fraction', 'moon_altitude', 'day_of_year')
    )
  }
  if(exclude_landcover){
    propensity_features <- setdiff(propensity_features,
                                   grep('mcd*|eastness_*|northness_*|elevation_*|astwbd_*|shoreline_*|road_density_*',
                                        propensity_features, value=TRUE))
  }

  treatment_features <- p$treatment_features %>% unlist()

  if(!aggregated){
    message("removing aggregated features")
    # remove aggregated features
    propensity_features <- setdiff(propensity_features,
                                   grep("mean_*", propensity_features, value=TRUE)
                                   )
    marginal_features <- setdiff(marginal_features,
                                    grep("mean_*", marginal_features, value=TRUE)
    )
    treatment_features <- setdiff(treatment_features,
                                    grep("mean_*", treatment_features, value=TRUE)
    )
    eff_covs <- setdiff(eff_covs,
                        grep("mean_*", eff_covs, value=TRUE))
  }
  return(list('propensity_features'=propensity_features,
              'marginal_features'=marginal_features,
              'treatment_features'=treatment_features,
              'eff_covs'=eff_covs))
}

#' @export
ebirdst_varnames <- function(vars){
  dat <- dplyr::tribble(
    ~var, ~varname,
    'mcd12q1_lccs1_c1_pland', 'barren',
    'mcd12q1_lccs1_c2_pland', 'permanent snow and ice',
    'mcd12q1_lccs1_c11_pland', 'evergreen needleleaf forests',
    'mcd12q1_lccs1_c12_pland', 'evergreen broadleaf forests',
    'mcd12q1_lccs1_c13_pland', 'deciduous needleleaf forests',
    'mcd12q1_lccs1_c14_pland', 'deciduous broadleaf forests',
    'mcd12q1_lccs1_c15_pland', 'mixed broadleaf/needleleaf forests',
    'mcd12q1_lccs1_c16_pland', 'mixed broadleaf evergreen/deciduous forests',
    'mcd12q1_lccs1_c21_pland', 'open forests',
    'mcd12q1_lccs1_c22_pland', 'sparse forests',
    'mcd12q1_lccs1_c31_pland', 'dense herbaceous',
    'mcd12q1_lccs1_c32_pland', 'sparse herbaceous',
    'mcd12q1_lccs1_c41_pland', 'dense shrublands',
    'mcd12q1_lccs1_c42_pland', 'shrubland/grassland mosaics',
    'mcd12q1_lccs1_c43_pland', 'sparse shrublands',
    'mcd12q1_lccs1_c255_pland', 'unclassified',
    'mcd12q1_lccs2_c9_pland', 'urban and built-up lands',
    'mcd12q1_lccs2_c25_pland', 'forest/cropland mosaics',
    'mcd12q1_lccs2_c35_pland', 'natural herbaceous/cropland mosaics',
    'mcd12q1_lccs2_c36_pland', 'herbaceous croplands',
    'mcd12q1_lccs3_c27_pland', 'woody wetlands',
    'mcd12q1_lccs3_c50_pland', 'herbaceous wetlands',
    'mcd12q1_lccs3_c51_pland', 'tundra',
    'mcd12q1_lccs1_c1_ed', 'ED barren',
    'mcd12q1_lccs1_c2_ed', 'ED permanent snow and ice',
    'mcd12q1_lccs1_c11_ed', 'ED evergreen needleleaf forests',
    'mcd12q1_lccs1_c12_ed', 'ED evergreen broadleaf forests',
    'mcd12q1_lccs1_c13_ed', 'ED deciduous needleleaf forests',
    'mcd12q1_lccs1_c14_ed', 'ED deciduous broadleaf forests',
    'mcd12q1_lccs1_c15_ed', 'ED mixed broadleaf/needleleaf forests',
    'mcd12q1_lccs1_c16_ed', 'ED mixed broadleaf evergreen/deciduous forests',
    'mcd12q1_lccs1_c21_ed', 'ED open forests',
    'mcd12q1_lccs1_c22_ed', 'ED sparse forests',
    'mcd12q1_lccs1_c31_ed', 'ED dense herbaceous',
    'mcd12q1_lccs1_c32_ed', 'ED sparse herbaceous',
    'mcd12q1_lccs1_c41_ed', 'ED dense shrublands',
    'mcd12q1_lccs1_c42_ed', 'ED shrubland/grassland mosaics',
    'mcd12q1_lccs1_c43_ed', 'ED sparse shrublands',
    'mcd12q1_lccs1_c255_ed', 'ED unclassified',
    'mcd12q1_lccs2_c9_ed', 'ED urban and built-up lands',
    'mcd12q1_lccs2_c25_ed', 'ED forest/cropland mosaics',
    'mcd12q1_lccs2_c35_ed', 'ED natural herbaceous/cropland mosaics',
    'mcd12q1_lccs2_c36_ed', 'ED herbaceous croplands',
    'mcd12q1_lccs3_c27_ed', 'ED woody wetlands',
    'mcd12q1_lccs3_c50_ed', 'ED herbaceous wetlands',
    'mcd12q1_lccs3_c51_ed', 'ED tundra',
    'astwbd_c1_pland', 'ocean',
    'astwbd_c2_pland', 'river',
    'astwbd_c3_pland', 'lakes',
    'astwbd_c1_ed', 'ED ocean',
    'astwbd_c2_ed', 'ED river',
    'astwbd_c3_ed', 'ED lakes',
    'gsw_c2_pland', 'seasonal water',
    'gsw_c3_pland', 'permanent water',
    'gsw_c2_ed', 'ED seasonal water',
    'gsw_c3_ed', 'ED permanent water',
    'road_density_c1', 'highways',
    'road_density_c2', 'primary roads',
    'road_density_c3', 'secondary roads',
    'road_density_c4', 'tertiary roads',
    'road_density_c5', 'local roads',
    'cds_u10', 'wind speed E/W',
    'cds_v10', 'wind speed N/S',
    'cds_i10fg', 'instantaneous wind gust',
    'cds_d2m', 'dewpoint temperature',
    'cds_t2m', 'temperature',
    'cds_cbh', 'cloud base height',
    'cds_hcc', 'high cloud cover',
    'cds_mcc', 'medium cloud cover',
    'cds_lcc', 'low cloud cover',
    'cds_msl', 'sea level pressure mean',
    'cds_slc', 'sea level pressure change',
    'cds_sf', 'snow water equivalent',
    'cds_rf', 'rainfall',
    'cds_tp', 'total precipitation',
    'bathymetry_elevation_median', 'Elevation summarized from 15-arcsecond grid (median)',
    'bathymetry_elevation_sd', 'Elevation summarized from 15-arcsecond grid (standard deviation)',
    'bathymetry_slope_median', 'Topographic slope summarized from 15-arcsecond grid (median)',
    'bathymetry_slope_sd', 'Topographic slope summarized from 15-arcsecond grid (standard deviation)'
  )
  if(any(vars %in% dat$var)){
    sapply(vars, function(var){
      if(var %in% dat$var){
        return(dat$varname[dat$var==var])
      }
      else{
        return(var)
      }
    }, USE.NAMES = FALSE)
  } else{
    sapply(vars, function(varname){
      if(varname %in% dat$varname){
        return(dat$var[dat$varname==varname])
      }else{
        return(varname)
      }
    })
  }

}


#' @param fit_predict output from `ebirdstwf::fit_model_stixel`
#' @param sp species code (maybe not necessary?)
#' @param erd2srd data frame that associates checklist_ids with srd_ids
#' @param srd_day which srd day do we want to predict to?
#' @param cci_count_varname one of `c_cci0`, `cci_occ_stx`, `prop_obs_gte0.50`, `mean_z_ppoisrrf`, or `none`
srd_pred_countcci <- function(fit_predict,
                              sp,
                              #erd2srd,
                              srd_day,
                              dir_srd_weekly = "~/data/erd_2023/srd_3km_week.parquet",
                              path_srd_year = "~/data/erd_2023/srd_3km_year.parquet",
                              cci_count_varname){
  stopifnot(dir.exists(dir_srd_weekly))
  stopifnot(file.exists(path_srd_year))
  stopifnot((cci_count_varname != "none") == (cci_count_varname %in% fit_predict$model$model_count$forest$independent.variable.names))

  m <- fit_predict
  sx <- unique(m$model_summary$stixel_id)
  # Check if srd preds exist already:
  f1 <- file.path('~/Documents/cci/cci_count/srd_preds',
            sprintf('%s_srd%s_%s_opt%s_%s.rds',
                    sp,
                    srd_day,
                    stringr::str_replace_all(cci_count_varname, "_", ""),
                    "TRUE",
                    sx))
  f2 <- file.path('~/Documents/cci/cci_count/srd_preds',
                  sprintf('%s_srd%s_%s_opt%s_%s.rds',
                          sp,
                          srd_day,
                          stringr::str_replace_all(cci_count_varname, "_", ""),
                          "FALSE",
                          sx))
  if( (file.exists(f1) & file.exists(f2)) | (cci_count_varname == "none" & file.exists(f2)) ){
    message("SRD file exists. Exiting now.")
    return()
  }

  # Generate params:
  param_objects <- c("PREDICTOR_LIST", "HABITAT_COVARS", "OBS_COVARS", "WEATHER_COVARS", "srd")
  params <- lapply(param_objects, function(pname){
    union_stixel_params(dir_interactive_stx = '~/data/erd_2023/results/interactive',
                        stxl_ids = sx,
                        pname)
  })
  names(params) <- param_objects

  # Setup for cci count values to predict to:
  cci_count_values <- list(c_cci0 = list(default_val = 1.85,
                                         valid_range = c(0,2)),
                           cci_occ_stx = list(default_val = 1.85,
                                              valid_range = c(0, 2)),
                           prop_obs_gte0.50 = list(default_val = 0.5,
                                                   valid_range = c(0, 1)),
                           prop_obs_gte0.75 = list(default_val = 0.25,
                                                   valid_range = c(0, 1)),
                           mean_z_ppoisrrf = list(default_val = 0.25,
                                                  valid_range = c(-2, 2)),
                           none = list(default_val = NA,
                                       valid_range = NA)
  )

  # Get SRD info and SRD data:

  # m$model$sampled <- erd2srd %>%
  #   mutate(checklist_id = as.character(checklist_id)) %>%
  #   inner_join(m$model$sampled)

  stopifnot(nrow(m$model$sampled) > 0)
  # rm(erd2srd); gc()

  # sp_params <- readRDS(sp_params_file)
  stixel_depth <-  (366 / 13) #if (params$IS_RESIDENT) 366 else (366 / 13)
  stixel_params <- ebirdstwf::parse_stixel_id(sx, stixel_depth = stixel_depth)
  # srd_days <- subset_days(m$model$sampled$closest_srd_day %>% unique(),
  #                         stixel_params$start_day, stixel_params$end_day)
  # stopifnot(srd_day %in% srd_days)

  ### SRD PREDS

  message("Reading in SRD files")
  # The problem with simply using sp_params$srd is that it might be missing some predictors that we actually fit the model on
  srd <- # weekly:
    arrow::read_parquet(file.path(dir_srd_weekly,
                                  sprintf("day_of_year=%s/part-0.parquet", srd_day))) %>%
  #   filter(srd_id %in% m$model$sampled$srd_id) %>%
  #   # annual:
  #   inner_join(arrow::read_parquet(path_srd_year))
    inner_join(params$srd %>%
    mutate(srd_id = checklist_id) %>%
    distinct(srd_id, .keep_all = TRUE))

  stopifnot(nrow(srd) > 0)

  # Dealing with CCI count predictive value:
  count_params <- cci_count_values[[cci_count_varname]]
  # have one version of optimized cci_count var, and one option that uses the
  # default value
  for(optimize in c(TRUE, FALSE)){
    message("species: ", sp, " optimize: ", optimize, " cci_count: ", cci_count_varname, " srd day: ", srd_day)
    if(optimize == TRUE & cci_count_varname == "none"){
      next()
    }
    outfile <- file.path('~/Documents/cci/cci_count/srd_preds',
                         sprintf('%s_srd%s_%s_opt%s_%s.rds',
                                 sp,
                                 srd_day,
                                 stringr::str_replace_all(cci_count_varname, "_", ""),
                                 optimize,
                                 sx))
    if(file.exists(outfile)){
      message("file exists. exiting now")
      return()
    }
    srd_effort <- ebirdstwf::optimize_srd_effort(model = m$model,
                                                 data = m$model$sampled,
                                                 srd_day = srd_day,
                                                 params = params,
                                                 optimize_cci_count = optimize,
                                                 cci_count_varname=cci_count_varname,
                                                 cci_count_default = count_params$default_val,
                                                 cci_count_valid_range = count_params$valid_range)

    cci_count_value <- srd_effort[[cci_count_varname]]
    # Predict to SRD using determined effort values
    message("predicting to SRD...")
    srd_preds_w <- predict_srd(model = m$model,
                               srd = srd, srd_day = srd_day,
                               srd_effort = srd_effort, params = params) %>%
      mutate(srd_day=srd_day,
             cci_count_value = cci_count_value,
             cci_count_varname = cci_count_varname,
             cci_count_optimized = as.character(optimize))

    message("Writing ", outfile)
    saveRDS(srd_preds_w,
            outfile)
  }
}

#' Combine up to multiple SRD predictions objects to a single raster that can be plotted
#' @export
srds2raster <- function(list_srd_predictions, srd_day){
  #browser()
  # one big srd pred object
  srd_preds <- do.call(dplyr::bind_rows, list_srd_predictions)

  # Template:
  r_w <- terra::rast("~/data/erd_2023/srd_3km_mask_ocean.tif")
  r_w[r_w == 1] <- NA

  r_w_occ <- r_w
  r_w_range <- r_w
  r_w_count <- r_w

  # fill in with values
  terra::values(r_w_occ)[srd_preds$checklist_id] <- srd_preds$pred_occ
  terra::values(r_w_range)[srd_preds$checklist_id] <- srd_preds$pred_range
  terra::values(r_w_count)[srd_preds$checklist_id] <- srd_preds$pred_count

  r_w_rel_abd <- r_w_occ * r_w_count
  r_w_rel_abd_thresh <- r_w_rel_abd
  r_w_range_thresh <- r_w_range
  r_w_range_thresh[r_w_range_thresh == 0] <- NA
  r_w_rel_abd_thresh <- terra::mask(r_w_rel_abd, r_w_range_thresh)
  terra::varnames(r_w_rel_abd_thresh) <- names(r_w_rel_abd_thresh) <- sprintf("rel_abd_thresh_%s", srd_day)
  terra::varnames(r_w_rel_abd) <- names(r_w_rel_abd) <-sprintf("rel_abd_%s", srd_day)
  terra::varnames(r_w_occ) <- names(r_w_occ) <- sprintf("occ_%s", srd_day)

  terra::varnames(r_w_count) <- names(r_w_count) <- sprintf("count_%s", srd_day)


  return(terra::rast(list(r_w_occ, r_w_count, r_w_rel_abd, r_w_rel_abd_thresh)))


}

#' @param rstrs raster object to plot
#' @param layer_to_plot index of raster to actually plot
#' @param plot_dir directory name where map should be written
#' @param type one of "rel_abd", "occ", or "count" (used for plot labels and file names only)
#' @param optimized whether count CCI was optimized when predicting to the SRD
#' @param cci_count_varname name of variable representing a count CCI
#' @param cci_count_val value of CCI `cci_count_varname` used in predicting to SRD
#' @param v_breaks vector of values to use to generate map color palette breaks. If NULL, will use values found in current raster layer.
map_raster_countcci <- function(rstr, layer_to_plot, species,
                        plot_dir,
                        type,
                        optimized,
                        cci_count_varname,
                        cci_count_vals,
                        v_breaks = NULL
){
  # mapping:
  r <- rstr[[layer_to_plot]]

  # Kind of irrelevant since we only plot one layer at a time:
  if(terra::nlyr(r) > 4){
    r_agg <- mean(r, na.rm=TRUE) %>%
      aggregate(fact=9)
  } else{
    r_agg <- r
  }

  ext_sinu <- ebirdstwf::raster_bounding_box(r_agg, crs=terra::crs(r_agg))
  crs <- r_agg |>
    ebirdstwf::raster_bounding_box() |>
    ebirdstwf::define_projection()
  ext_crs <- ebirdstwf::raster_bounding_box(r_agg, crs=crs)

  v <- r %>%
    terra::crop(ext_sinu) %>%
    terra::values(na.rm=TRUE) %>%
    as.vector()

  v <- v[v>0]
  if(is.null(v_breaks)){
    breaks <- quantile(v, probs = seq(0, 1, by=0.05), na.rm=TRUE)
  } else{
    breaks <- quantile(v_breaks, probs = seq(0, 1, by=0.05), na.rm=TRUE)
  }

  pal <- ebirdstwf::map_palette(length(breaks) - 1, type = "abundance")
  labels <- c(breaks[2], median(breaks), breaks[length(breaks) - 1])
  rt <- ebirdstwf::generate_raster_template(crs = crs$proj4string,
                                 ext = terra::ext(ext_crs),
                                 res = terra::res(r))
  basemap <- ebirdstwf::projected_basemap(crs)

  # map a single layer (week)
  r_proj <- ebirdstwf::project_raster(r, rt$template, method = "near")

  png(file.path(plot_dir, sprintf("%s_%s_srd_%s_%s_opt%s.png", species, type,
                                  names(r),
                                  stringr::str_replace_all(cci_count_varname, '_', ''),
                                  optimized)
                ),
      width = 2400, height = 2400)
  par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))

  ebirdstwf::create_basemap(basemap, bbox = sf::st_bbox(rt$map_bbox)*1.2,
                 # this will add the background zeros in a light grey
                 zeros = NULL)#r_proj)

  # the non-zero abundance data
  terra::plot(r_proj, col = pal, breaks =breaks,
       maxcell = terra::ncell(r_proj),
       legend = FALSE, axes = FALSE, bty = "n",
       reset = FALSE, add = TRUE)

  # state/country lines
  ebirdstwf::add_boundaries(basemap)

  # labels and legend
  ebirdstwf::add_map_text(sprintf('CCI count metric: %s', cci_count_varname),
               x = 0.020, y = (1-0.035),#0.17,
               cex = 6, pos = 4)
  ebirdstwf::add_map_text(sprintf("Optimized: %s; %s", optimized, paste(round(sort(cci_count_vals),1), collapse=", ")),
               x = 0.020, y = (1-0.08),#0.125,
               cex = 6, pos = 4)
  ebirdstwf::add_map_text(names(r),
               x = 0.020, y = 0.08, cex = 6, pos = 4)
  ebirdstwf::add_map_text(species,#"Wood Thrush",
               x = 0.020, y = 0.035, cex = 9, pos = 4, font = 2)
  ebirdstwf::add_map_legend(palette = pal,
                 title = type,#"Probability of detection",
                 labels = labels,
                 # use this if you're using min not 5th percentile as bottom label
                 # less_than = FALSE,
                 x = 0.88, y = 0.25, width = 0.03, height = 0.5,
                 scale_axis_text=2.5)

  dev.off()

}
