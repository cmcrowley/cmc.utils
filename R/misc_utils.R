#' get_constant_vars
#' Return variable names of columns that are constant across `dat`
#' @export
get_constant_vars <- function(dat){
  is_constant <- sapply(names(dat),
         function(varx){
           constant <- length(unique(dat[[varx]])) == 1
           }
         )
  return(names(is_constant[which(is_constant)]))
}

#' @export
inverse_logit <- function(logit){
  1 / (exp(-logit) + 1)
}

#' @export
geometric_mean <- function(x){
  n <- length(x)
  prod(x)^(1/n)
}

#' @export
round_to_nearest <- function(v, constant=NULL, logbase=NULL,
                             ties='random') {
  sapply(v, function(x) {
    # TODO handle ties in a principled way (disabled for now)
    stopifnot(ties %in% c('floor', 'ceiling', 'random'))
    if(is.null(constant) == is.null(logbase)) {
      stop("Must specify either `constant` or `base`, but not both")
    }

    if(!is.null(logbase)) {
      if(x == 0){
        return(0)
      }
      # Cant have negative numbers for a log-based thing
      # Could see argument for returning 0 instead
      if(x == -1){
        # treat this like you treat a value of 1
        x <- 1
      }
      if (x < -1){
        return(NA)
      }
      # Determine constant
      # down <- logbase^floor(log(x, base=logbase))
      # up <- logbase^ceiling(log(x, base=logbase))
      constant <- logbase^round(log(x, base=logbase))
    }

    # if(x %% constant == 0.5(constant) & ties != 'floor') {
    #   up <- ceiling(x/constant)*constant
    #   if(ties == 'ceiling'){
    #     return(up)
    #   } else{
    #     # If x is halfway between two values & ties are random,
    #     # there should be an equal chance of rounding up or rounding down:
    #     down <- floor(x/constant)*constant
    #     warning(sprintf("x= %s falls exactly between two eligible values. Picking randomly between %s and %s", x, down, up))
    #     return(sample(c(down, up), size=1))
    #   }
    #
    # } else{
    #   # This equates to "floor" option:
    round(x/constant)*constant
    # }
  })

}

#' @description Found at: https://stackoverflow.com/questions/520810/does-r-have-quote-like-operators-like-perls-qw
#' @export
names2chr <- function(...){
  sapply(match.call()[-1], deparse)
}


#' For each value in one data set, compute its rarity in another
#' @param x data.frame of interest
#' @param dist data.frame from which to compute frequencies
data_ecdf <- function(x, dist){
    # For character or factor variables, compute the proportion in `dist` that match
    # For numeric variables, compute the percentile (ecdf)
   # First make sure all the types match
   # pass
  chr_ecdf_col <- function(xi, col){
    stopifnot(length(xi)==1)
    stopifnot(class(xi)==class(col))

    t <- table(col)
    v <- as.numeric(t)
    v <- v/(length(v))

    return(v[match(xi, names(t))])
  }
  # Then do the character variables
  names_chrs <- names(which(sapply(dist, function(col){
    class(col) %in% c('character', 'factor', 'logical')
  }) ))

  # then do the numeric variables
  names_num <-  names(which(sapply(dist, function(col){
    class(col) %in% c('numeric', 'integer')
  }) ))

  if(length(names_chrs) > 0){

  }
}





##### caseweighting functions
#' Create weights inversely proportional to relative class size
#'
#' @param class_col character vector of classes
#' @export
stratify2weights <- function(class_col){
  tab <- prop.table(table(class_col))

  return((1-sapply(class_col, function(class){tab[[class]]})) / min(tab))
}

#' @export
joint_caseweights <- function(dat, stratify_col, caseweight_col){
  message_time("Computing joint case weights of ", caseweight_col, " stratified by ", stratify_col)
  stopifnot(length(stratify_col) == 1 & length(caseweight_col) == 1)
  stratify_ws <- stratify2weights(dat[[stratify_col]])
  final_weights <- dat[[caseweight_col]]*stratify_ws
  return(final_weights)
}

#' Calculate an "effective" sample size when using weights
#'
#' @export
n_effective <- function(w){
  sum(w)^2 / sum(w^2)
}



# Partial dependence function (without the summarize step)
# TODO: make it so that there's a flexible response
# logit vs prob vs class
# currently, can only have one grid sise even if there are multiple variables to jointly conduct PD
#' @export
ice <- function(dat, mod, vars, grid_size){
  stopifnot( all(vars %in% names(dat)) )
  datmod <- dat %>% dplyr::select(-all_of(vars))

  # get the artificial vectors of variables of interest to compute over:
  vals <- lapply(vars, function(var){
    seq(from=min(dat[[var]], na.rm=TRUE),
        to=max(dat[[var]], na.rm=TRUE),
        length.out = grid_size)
  })
  # What's not clear to me: Does PD(x, y) mean that we compute over the combination of X and Y
  # For the moment: do the exhuastive combos while waiting for an answer on S.O.
  combos <- expand.grid(vals)
  names(combos) <- vars

  message("Generating data over combinations of vars...")
  pd_dat <- do.call(rbind, lapply(1:nrow(combos), function(i){
    d <- datmod
    for(j in 1:length(vars)){
      vn <- vars[j]
      d[[vn]] <- combos[[vn]][j]
    }
    return(d)
  }))

  # TODO: use arg predict_fun if we wanto to not be specific to ranger?
  message('Computing predictions...')
  preds <- cbind(combos, predict(mod, pd_dat)$predictions)
  names(preds) <- c(names(combos), sprintf('prob_%s', mod$forest$class.values))

  return(list('ice'=preds))
}
