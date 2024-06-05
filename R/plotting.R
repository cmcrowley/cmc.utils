#' @export
ggtheme <- function(p, scale_x_date=TRUE, text_multiplier = 1){
  # title_text_col <- "#222222"
  text_col <- text_color()
  if(scale_x_date){
    p <- p + ggplot2::scale_x_date(date_breaks='1 year',
                                   date_labels="'%y")
  }
  p <- p + ggplot2::theme(axis.title=ggplot2::element_text(size=15*text_multiplier,
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
  
  return(p)
}

#' @export
text_color <- function(){
  '#434545'
}

#' @name plot_var_imp_rf
#' @description Produce a ggplot of variable importances. If length(`rf_list`) > 1, will produce a multi-faceted plot with a panel for each `rf`
#' @param rf_list list of `ranger` model objects
#' @param n_top (int) subset of variables for which to plot importances
#' @param ... other arguments to pass to ggtheme()
#' @export
plot_var_imp_rf <- function(rf_list, n_top=NULL, facet_scales="fixed",
                            names_function=NULL, ...) {
  # do.call(rbind, lapply(rf_list, function(rf){}))
  stopifnot(is.list(rf_list))
  
  df_all <- do.call(rbind, lapply(seq_along(rf_list), function(i){
    df <- get_var_imp_rf(rf_list[[i]]) %>%
      slice_max(order_by = importance_value, n=n_top)
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
  
  ordering <- df_all %>%
    group_by(rf_name) %>%
    slice_max(order_by = importance_value, n=n_top) %>%
    droplevels() %>%
    do(data.frame(al=levels(reorder(interaction(.$rf_name, .$var_name, drop=TRUE),
                                    .$importance_value)))) %>%
    pull(al)
  
  varname <- gsub("^.*\\.", "", ordering)
  
  # R^2 labels and x coordinate of placement:
  label_importance <- df_all %>%
    group_by(rf_name) %>%
    summarize(median_importance = median(importance_value))
  r2_labels <- df_all %>%
    distinct(rf_name, .keep_all=TRUE) %>%
    select(rf_name, r2_oob) %>%
    full_join(label_importance, by='rf_name')
  
  p <- df_all %>%
    mutate(var = factor(interaction(rf_name, var_name), levels=ordering)) %>%
    # drops unused levels in a way that actually works here:
    filter(!is.na(var)) %>%
    ggplot(aes(x=round(importance_value), y = var)) +
    facet_wrap(~rf_name, scales=facet_scales) +
    # geom_line(aes(color=var_name, group=var_name)) +
    geom_dotplot(binaxis = 'y',
                 aes(fill=var_name, color=var_name)) +
    scale_y_discrete(breaks=ordering, labels = varname) +
    ggplot2::geom_text(data=r2_labels,
                       # TODO: this positioning doesn't work. Fix!
                       ggplot2::aes(label=sprintf("R2=\n %s",
                                                  as.character(r2_oob)),
                                    x=median_importance),
                       
                       y=3) +
    ggplot2::labs(y='Variable',
                  x='Importance Value',
                  title = sprintf("Top %s variables by model", n_top))
  
  return(ggtheme(p, scale_x_date=FALSE, ...) +
           ggplot2::theme(legend.position='none') +
           ggplot2::theme(axis.text.y = element_text(hjust=1))
  )
}
