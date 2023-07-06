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

#' @description Found at: https://stackoverflow.com/questions/520810/does-r-have-quote-like-operators-like-perls-qw
#' @export
names2chr <- function(...){
  sapply(match.call()[-1], deparse)
}


#' @name formula2terms
#' @param formula `formula` or `character`
#' @export
formula2terms <- function(formula){
  return(rownames(attr(terms.formula(as.formula(formula)), "factors")))
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
plot_var_imp_rf <- function(rf_list, n_top=NULL){
  # do.call(rbind, lapply(rf_list, function(rf){}))
  stopifnot(is.list(rf_list))

  df_all <- do.call(rbind, lapply(seq_along(rf_list), function(i){
    df <- get_var_imp_rf(rf_list[[i]])
    df$rf_name <- names(rf_list)[i]
    return(df)
  }))

  if(!all(is.null(df_all$rf_name))){
    df_all$rf_name <- factor(df_all$rf_name, levels=names(rf_list), ordered=TRUE)
  }


  # Compute the mean importance across all rfs to determine facet ordering
  df_all <- df_all %>%
    dplyr::group_by(var_name) %>%
    dplyr::mutate(mean_imp = mean(importance_value)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(var_name=factor(var_name, levels=unique(var_name[order(mean_imp)])))

  # Narrow down to top_n
  if(is.null(n_top)){
    n_top <- length(unique(df_all$var_name))
  }
  vars_include <- df_all %>% dplyr::arrange(desc(mean_imp)) %>% dplyr::distinct(var_name) %>% dplyr::pull(var_name)
  vars_include <- vars_include[1:n_top] %>% droplevels()
  df_sub <- df_all %>% dplyr::filter(var_name %in% vars_include) %>% droplevels()

  p <- ggplot2::ggplot(df_sub, ggplot2::aes(x=round(importance_value), y=var_name)) +
    ggplot2::geom_dotplot(binaxis = "y", ggplot2::aes(fill=var_name, color=var_name)) +
    ggplot2::facet_wrap(~rf_name, nrow=2) +
    ggplot2::theme(legend.position='bottom') +
    ggplot2::labs(y='Variable',
         x='Importance Value',
         title= "RF importances by tixel",
         subtitle="for top 10 mean importance across tixels")

  return(p)
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
  p + ggplot2::theme(axis.title=ggplot2::element_text(size=14*text_multiplier,
                                                      color=text_col,
                                                      family="Helvetica"#,
                                                      #face='bold'
  ),
  axis.text.y=ggplot2::element_text(size=13*text_multiplier,
                                    color=text_col,
                                    family="Helvetica"),
  axis.text.x=ggplot2::element_text(size=9*text_multiplier,
                                    color=text_col,
                                    family="Helvetica"),
  plot.title=ggplot2::element_text(family="Helvetica",
                                   face="bold",
                                   color=text_col,
                                   size=20*text_multiplier,
                                   hjust=0),
  plot.subtitle=ggplot2::element_text(family="Helvetica",
                                      color=text_col,
                                      hjust=0,
                                      size=15*text_multiplier,
                                      margin=ggplot2::margin(t=4, b=7)),
  strip.text.x = ggplot2::element_text(size=13.5*text_multiplier, hjust=0,
                                       color=text_col,
                                       family="Helvetica"),
  strip.background = ggplot2::element_rect(fill=NA,#'grey95',
                                           color=NA),
  panel.background = ggplot2::element_rect(fill='white', color=NA),
  panel.border = ggplot2::element_rect(fill=NA, color='#545454'),
  panel.grid=ggplot2::element_line(color='grey92'),
  legend.key = ggplot2::element_rect(fill='white', color=NA),
  legend.text=ggplot2::element_text(color=text_col,
                                    size=13*text_multiplier),
  legend.title=ggplot2::element_text(color=text_col,
                                     size=14*text_multiplier),
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
get_date_from_doy <- function(year, doy){
  as.Date(doy -1, origin=sprintf('%d-01-01', year))
}
