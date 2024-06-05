#' @export
r_squared <- function(observed, predicted){
  ssr <- sum( (observed - predicted)^2, na.rm=TRUE )
  sst <- sum( (observed - mean(observed))^2, na.rm=TRUE )
  1 - (ssr/sst)
}

#' @description A residual weighting that downweights by the magnitudes of the values.
#' Good for when residual weights should roughly half in magnitude as order of magnitude in y increases by 1 (log10 scale)
#' For example, a residual of 3 is weighted higher for 4-1 than it does for 43-40. It also needs to be symmetric, that is,
#' a residual of 4-1 should have the same magnitude as a residual of 1-4, which is why we weight by the mean of y and yhat (rather than by either one).

wresid3 <- function(yhat, y){
  r = y-yhat
  mean_y_yhat <- sapply(1:length(yhat), function(i){mean(c(yhat[i], y[i])) })
  r/(log10(mean(c(y, yhat)) + 1)) # The plus one makes it a pseudo 'base-1' logarithm.
}


wresid4 <- function(yhat, y, base){
  r = y-yhat
  # As the residual gets larger, the asymmetry will be larger.
  r/(log(yhat+1, base=base)) # The plus one makes it a pseudo 'base-1' logarithm.
}


#' Difference as a fraction of expected
#' @export
mpe <- function(yhat, y){
  yhat1 <- yhat+1
  (y-yhat)/(yhat1)
}


#' @export
mpe_inverse <- function(mpe, yhat){
  mpe*(yhat+1) + yhat
}

# Scale so that percent errors are weighted higher for larger numbers
#' @export
mpe_logscaled <- function(yhat, y){
  yhat1 <- yhat+1
  log(yhat1)*(y-yhat)/(yhat1)
}

#' @export
logratio <- function(yhat, y){
  log((y+1)/(yhat+1))
}

#' @export
logratio_inverse <- function(logratio, yhat){
  exp(logratio)*(yhat+1)-1
}

#' Calculate a weighted mean and a standard error of the mean in multiple ways
#'
#' @export
weighted_mean_se <- function(x, obs_ss, type='2', cap_ss_value=NULL){
  # x = valid cci estimates of length n_folds
  # obs_ss = observer sample size that went into that estimate
  # type 1 true weighted mean, where we calculate an effective sample size based on the weights
  # type 2 each point gets mostly the same weight, but the uncertainties on each are very different
  
  
  ####################
  # Create the weights
  ####################
  #The optimal allocation of weights is to make each point's weight inversely proportional to the square of its
  # uncertainty
  # ie, in inverse proportion to its variance
  #(lower sample size =higher uncertainty = lower weight)
  # Create weights that are capped if specified
  obs_ss <- obs_ss[!is.na(x)]
  x <- x[!is.na(x)]
  if(!is.null(cap_ss_value)){
    obs_ss[obs_ss > cap_ss_value] <- cap_ss_value
  }
  #observe the measured uncertainties
  # w = 1/variance --> variance = 1/weight
  # Assuming CCI is in and of itself a mean of IID normal variables (z-scores)
  # (which is not quite true at least because there's correlation, not independence)
  # -- to deal with this, change sqrt(n) in the denominator to something smaller
  # then we should assume the variance decreases as the sqrt(n) increases
  w <- sqrt(obs_ss)
  
  ################
  # Sample sizes
  ################
  # I think we should use type 2 but we can try both
  if(type=='1'){
    n = n_effective(w) # basically dealing across folds and weighted-summing the sample sizes within folds
  } else if(type=='2'){
    n = length(x[!is.na(x)]) # counting only the valid folds
  } else{
    stop('Must specify type "1" or "2"')
  }
  correction = n / (n-1)
  
  xbar <- stats::weighted.mean(x, w)
  variance_i <- correction*(sum( w*( (x-xbar)^2 ) / sum(w) ))
  variance_mean <- variance_i / n
  return(list(
    cci_estimate=xbar,
    cci_se = sqrt(variance_mean)
  )
  )
}