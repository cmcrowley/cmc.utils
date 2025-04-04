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
