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
  dat <- tribble(
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
    'mcd12q1_lccs1_c13_ed', 'EDdeciduous needleleaf forests',
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
    'mcd12q1_lccs2_c25_ed', 'ED forest/croed mosaics',
    'mcd12q1_lccs2_c35_ed', 'ED natural herbaceous/croed mosaics',
    'mcd12q1_lccs2_c36_ed', 'ED herbaceous croeds',
    'mcd12q1_lccs3_c27_ed', 'ED woody wetlands',
    'mcd12q1_lccs3_c50_ed', 'ED herbaceous wetlands',
    'mcd12q1_lccs3_c51_ed', 'ED tundra',
    'astwbd_c1_pland', 'ocean',
    'astwbd_c2_pland', 'river',
    'astwbd_c3_pland', 'lakes',
    'gsw_c2_pland', 'seasonal water',
    'gsw_c3_pland', 'permanent water',
    'gsw_c2_ed', 'ED seasonal water',
    'gsw_c3_ed', 'ED permanent water',
    'road_density_c1', 'highways',
    'road_density_c2', 'primary roads',
    'road_density_c3', 'secondary roads',
    'road_density_c4', 'tertiary roads',
    'road_density_c5', 'local roads',
    'cds_u10', 'wind speed e/w',
    'cds_v10', 'wind speed n/s',
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
    'cds_tp', 'total precipitation'
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
