#' @name get_missing_dates
#' @description Given a vector of dates, return any that are 'skipped', ie missing. Assumes neither the first nor the last date in the sequence is missing.
#' @export
get_missing_dates <- function(dates){
  #TODO test
  #TODO handle NAs in date vector
  #TODO handle when `dates` is not class Date
  stopifnot(class(dates) == 'Date')
  seqd <- seq.Date(from=min(dates), to=max(dates), by=1)
  return(seqd[!(seqd %in% dates)])
}

#' get quantiles on vector of dates
#'
#' @param dates Date vector
#' @param probs numeric; quantiles for which you want exceedance probabilities
#'
#' @export
get_date_quantiles <- function(dates, probs){
  as.Date(quantile(as.numeric(dates), probs), origin="1970-01-01")
}

#' Convert a vector of dates to values on [0,1]
#' @param dates; Date
#'
#' @export
date_ecdf <- function(dates){
  n <- as.numeric(dates)
  cdf <- ecdf(n)
  cdf(n)
}

#' @export
doy2date <- function(year, doy){
  as.Date(doy -1, origin=sprintf('%d-01-01', year))
}

#' @export
cycle_date <- function(d){
  sapply(d, function(d){
    if(d > 366){
      return(d-366)
    } else{
      return(d)
    }
  })
}