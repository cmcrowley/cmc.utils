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
