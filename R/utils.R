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

#' Convert data frame to sf object with ESPG 4326 SRD
#' @export
dat2sf <- function(dat, coord_cols=c('longitude', 'latitude')){
  dat %>%
    sf::st_as_sf(coords=coord_cols) %>%
    sf::st_set_crs("+proj=longlat +datum=WGS84")
}

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

#' @description Found at: https://stackoverflow.com/questions/520810/does-r-have-quote-like-operators-like-perls-qw
#' @export
names2chr <- function(...){
  sapply(match.call()[-1], deparse)
}

#' Parse individual variable names from formula terms
#'
#' (intended to handle GAM arguments like s(foo, bs="re", "k=10"))
#'
terms2vars <- function(terms){
  sapply(terms, function(term){
    stringr::str_replace_all(terms,
                         c("s\\("= "",
                           ",\\s*[a-z]+\\s*=\\s*.+\\s*\\)$"="")
    ) %>%
      stringr::str_split(pattern="\\s*,\\s*") %>%
      unlist()
  })

}

#' @name formula2terms
#' @param formula `formula` or `character`
#' @export
formula2terms <- function(formula){
  # return(rownames(attr(terms.formula(as.formula(formula)), "factors")))
  f.gam <- mgcv:::interpret.gam(formula)
  return(c(f.gam$pred.names, f.gam$response))
}

#' Get the number of out-of-bag trees for each data point
#' @param rf ranger random forest object (with inbag.counts)
#'
#' @export
get_ntrees_oob <- function(rf){
  # nobs x ntrees matrix
  inbag_mat <- do.call(cbind, rf$inbag.counts)
  # vector of size nobs. for each obs, the number of trees for which that obs was out-of-bag
  # in other words, the number of trees that contribute to its prediction value
  rowSums(inbag_mat == 0)
}


#' @name get_var_imp_rf
#' @description Get variable importances from a `ranger` model object as a d
#' @param rf `ranger` model object
#' @value a data frame with names `var_name`, `importance_type`, and `importance_value`
#' @export
get_var_imp_rf <- function(rf){
  imp_type <- rf$importance.mode
  imp_v    <- rf$variable.importance
  df       <- data.frame(var_name=names(imp_v), importance_value=as.numeric(imp_v), importance_type=imp_type)
  return(df[order(df$importance_value, decreasing = TRUE),])
}



#' @name plot_var_imp_rf
#' @description Produce a ggplot of variable importances. If length(`rf_list`) > 1, will produce a multi-faceted plot with a panel for each `rf`
#' @param rf_list list of `ranger` model objects
#' @param n_top (int) subset of variables for which to plot importances
#' @export
plot_var_imp_rf <- function(rf_list, n_top=NULL, facet_scales="fixed", names_function=NULL){
  # do.call(rbind, lapply(rf_list, function(rf){}))
  stopifnot(is.list(rf_list))

  df_all <- do.call(rbind, lapply(seq_along(rf_list), function(i){
    df <- get_var_imp_rf(rf_list[[i]])
    df$rf_name <- names(rf_list)[i]
    df$r2_oob <- ifelse('r.squared' %in% names(rf_list[[i]]),
                        round(rf_list[[i]]$r.squared, 3),
                        round(1-rf_list[[i]]$prediction.error, 3))
    return(df)
  }))

  if(!all(is.null(df_all$rf_name))){
    df_all$rf_name <- factor(df_all$rf_name, levels=names(rf_list), ordered=TRUE)
  }

  if(!is.null(names_function)){
    df_all$var_name <- names_function(df_all$var_name)
  }
  # Compute the mean importance across all rfs to determine variable ordering
  df_all <- df_all %>%
    dplyr::group_by(var_name) %>%
    dplyr::mutate(mean_imp = mean(importance_value, na.rm=TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(var_name=factor(var_name, levels=unique(var_name[order(mean_imp)])))

  # Narrow down to top_n
  if(is.null(n_top)){
    n_top <- length(unique(df_all$var_name))
  }
  vars_include <- df_all %>% dplyr::arrange(desc(mean_imp)) %>% dplyr::distinct(var_name) %>% dplyr::pull(var_name)
  vars_include <- vars_include[1:n_top] %>% droplevels()
  df_sub <- df_all %>% dplyr::filter(var_name %in% vars_include) %>% droplevels()

  labels <- df_sub %>% distinct(rf_name, .keep_all=TRUE) %>% select(rf_name, r2_oob)
  max_imp <- max(df_sub$importance_value)
  p <- ggplot2::ggplot() +
    ggplot2::geom_dotplot(data=df_sub, binaxis = "y",
                          ggplot2::aes(x=round(importance_value), y=var_name,
                                       fill=var_name, color=var_name)) +
    ggplot2::geom_text(data=labels, ggplot2::aes(label=sprintf("R2=\n %s", as.character(r2_oob))),
                       x=0.8*max_imp,
                       y=3) +
    ggplot2::facet_wrap(~rf_name, nrow=2, scales=facet_scales) +
    ggplot2::labs(y='Variable',
         x='Importance Value')

  return(ggtheme(p, scale_x_date=FALSE) + ggplot2::theme(legend.position='none')
)
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



# x <- sort(round(runif(50, 1, 100)))
# y <- x*log(x)
# y <- log(x)/x
# plot(x, y)
# abline(a=0, b=1)
#
# y <- sort(round(runif(50, 1, 100)))
# r <- c(0, 1, 10, 20, 50, -50, -20, -10, -1)
# df <- expand.grid(list(y, r)); names(df) <- c('y', 'r')
# df$yhat <- df$y - df$r
# df <- df %>% filter(yhat >= 0)
#
# df <- df[, c(1, 3, 2)]
# df$rfrac <- mpe(df$yhat, df$y)
# df$logrfrac <- mpe_logscaled(df$yhat, df$y)
# df$logratio <- log(df$y/df$yhat)
# View(df)
#
# df %>%
#   mutate(resid=r) %>%
#   tidyr::pivot_longer(cols=c('rfrac', 'logrfrac', 'logratio'), names_to='error_metric', values_to = 'error_value') %>%
#   ggplot(aes(x=yhat, y=error_value)) +
#   geom_point(aes(color=error_metric)) +
#   facet_wrap(~r, scales='free')

#' @export
ggtheme <- function(p, scale_x_date=TRUE, text_multiplier = 1){
  # title_text_col <- "#222222"
  text_col <- text_color()
  if(scale_x_date){
    p <- p + ggplot2::scale_x_date(date_breaks='1 year',
                                   date_labels="'%y")
  }
  p + ggplot2::theme(axis.title=ggplot2::element_text(size=15*text_multiplier,
                                                      color=text_col,
                                                      family="Helvetica"#,
                                                      #face='bold'
  ),
  axis.text.y=ggplot2::element_text(size=14*text_multiplier,
                                    color=text_col,
                                    family="Helvetica"),
  axis.text.x=ggplot2::element_text(size=14*text_multiplier,
                                    color=text_col,
                                    family="Helvetica"),
  plot.title=ggplot2::element_text(family="Helvetica",
                                   face="bold",
                                   color=text_col,
                                   size=20*text_multiplier,
                                   hjust=0,
                                   # top,right,bottom,left
                                   margin=ggplot2::margin(0,0,5,0)),
  plot.subtitle=ggplot2::element_text(family="Helvetica",
                                      color=text_col,
                                      hjust=0,
                                      size=16*text_multiplier,
                                      margin=ggplot2::margin(t=4, b=7)),
  strip.text.x = ggplot2::element_text(size=14*text_multiplier, hjust=0,
                                       color=text_col,
                                       family="Helvetica"),
  strip.text.y = ggplot2::element_text(size=14*text_multiplier,
                                       color=text_col,
                                       family="Helvetica",
                                       angle=-90,
                                       hjust=0.5),
  strip.background = ggplot2::element_rect(fill=NA,#'grey95',
                                           color=NA),
  panel.background = ggplot2::element_rect(fill='white', color=NA),
  panel.border = ggplot2::element_rect(fill=NA, color='#545454'),
  panel.grid=ggplot2::element_line(color='grey92'),
  legend.key = ggplot2::element_rect(fill='white', color=NA),
  legend.text=ggplot2::element_text(color=text_col,
                                    size=14*text_multiplier),
  legend.title=ggplot2::element_text(color=text_col,
                                     size=15*text_multiplier),
  complete=TRUE)
}

#' @export
text_color <- function(){
  '#434545'
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


#' @export
doy2date <- function(year, doy){
  as.Date(doy -1, origin=sprintf('%d-01-01', year))
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
