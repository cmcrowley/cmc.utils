#' @name get_missing_dates
#' @description Given a vector of dates, return any that are 'skipped', ie missing. Assumes neither the first nor the last date in the sequence is missing.
#' @export
get_missing_dates <- function(dates){
  #TODO test
  #TODO handle NAs in date vector
  #TODO handle when `dates` is not class Date
  stopifnot(is.Date(dates))
  seqd <- seq.Date(from=min(dates), to=max(dates), by=1)
  return(seqd[!(seqd %in% dates)])
}


#' @name formula2terms
#' @export
formula2terms <- function(formula){
  return(rownames(attr(terms.formula(formula), "factors")))
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
#' @param top_n (int) subset of variables for which to plot importances
#' @export
plot_var_imp_rf <- function(rf_list, top_n=NULL){
  # do.call(rbind, lapply(rf_list, function(rf){}))
  stopifnot(is.list(rf_list))

  df_all <- do.call(rbind, lapply(seq_along(rf_list), function(i){
    df <- get_var_imp_rf(rf_list[[i]])
    df$rf_name <- names(rf_list)[i]
    return(df)
  }))

  df_all$rf_name <- factor(df_all$rf_name, levels=names(rf_list), ordered=TRUE)

  # Compute the mean importance across all rfs to determine facet ordering
  df_all <- df_all %>%
    dplyr::group_by(var_name) %>%
    dplyr::mutate(mean_imp = mean(importance_value)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(var_name=factor(var_name, levels=unique(var_name[order(mean_imp)])))

  # Narrow down to top_n
  vars_include <- df_all %>% dplyr::arrange(desc(mean_imp)) %>% dplyr::distinct(var_name) %>% dplyr::pull(var_name)
  vars_include <- vars_include[1:top_n] %>% droplevels()
  df_sub <- df_all %>% dplyr::filter(var_name %in% vars_include) %>% droplevels()

  p <- ggplot2::ggplot(df_sub, ggplot2::aes(x=round(importance_value), y=var_name)) +
    ggplot2::geom_dotplot(binaxis = "y", ggplot2::aes(fill=var_name, color=var_name)) +
    ggplot2::facet_wrap(~rf_name, nrow=2) +
    ggplot2::theme(legend.position='bottom') +
    ggplot2::labs(y='Variable',
         x='Importance Value',
         title= "RF importances by tixel",
         subtitle="for top 10 mean importance across tixels")
}
