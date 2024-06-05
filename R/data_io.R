#' Convert data frame to sf object with ESPG 4326 SRD
#' @export
dat2sf <- function(dat, coord_cols=c('longitude', 'latitude')){
  dat %>%
    sf::st_as_sf(coords=coord_cols) %>%
    sf::st_set_crs("+proj=longlat +datum=WGS84")
}

#' Get a data frame from a .db file
#
#' @param dir_db path to .db file
#'
#' @return data frame or list of data.frames
#' @export
read_db_table <- function(path_db){
  con <- DBI::dbConnect(RSQLite::SQLite(), path_db)
  tables <- DBI::dbListTables(con)
  if(length(tables)==1){
    dat <- DBI::dbReadTable(con, tables)
  } else{
    stop("Write me to read more than one table!")
  }
  DBI::dbDisconnect(con)
  return(dat)
}