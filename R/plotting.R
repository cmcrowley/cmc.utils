#' plot_pd
#' @param pd data frame with columns "predictor", "x", "summarized_response", "q25_response"
#' Returns a partial dependence plot given output from `calculate_pd()`
#' @export
plot_pd <- function(pd, include_iqr = TRUE){

  p <- ggplot(pd) +
    aes(x = value, y = response_mean) +
    geom_point() +
    geom_line() +
    geom_line(aes(y = q_0.05),
              linetype = "dashed") +
    geom_line(aes(y = q_0.25)) +
    geom_line(aes(y = q_0.75)) +
    geom_line(aes(y = q_0.95),
              linetype = "dashed") +
    scale_color_viridis_c() +
    facet_wrap(~ predictor, ncol = 1, scales = "free") +
    theme(legend.position = "bottom",
          legend.key.width = unit(0.2, "npc"))
return(p)
}

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
#' Produce a ggplot of variable importances. If length(`rf_list`) > 1,
#' will produce a multi-faceted plot with a panel for each `rf`. In this case to
#' ease comparison among objects, points are colored by variable name.
#' @param rf_list list of `ranger` model objects
#' @param n_top (int) subset of variables for which to plot importances
#' @param ... other arguments to pass to ggtheme()
#' @export
plot_var_imp_rf <- function(rf_list, n_top=NULL, facet_scales="free",
                            names_function=NULL, text_multiplier = 1,...) {
  stopifnot(is.list(rf_list))

  # Determine depth of list (& therefore structure of plot)
  if(class(rf_list[[1]]) == "list"){
    multilevel <- TRUE
    # We have two levels of nesting, so we'll use facet_grid,
    # taking the first level as the columns and second level as the rows
    cols <- names(rf_list)
    rows <- names(rf_list[[1]])
    # Check that each element of cols() has the same set of rows
    all_second_level <- unique(unlist(lapply(rf_list, names), use.names = FALSE))
    if(!setequal(rows, all_second_level)){
      error("mismatch in units")
    }
  } else{
    multilevel <- FALSE
  }

  parse_rfl <- function(rfl){
    df_all <- do.call(rbind, lapply(seq_along(rfl), function(i){
      df <- get_var_imp_rf(rfl[[i]]) %>%
        slice_max(order_by = importance_value, n=n_top)
      df$row_var <- names(rfl)[i]
      df$r2_oob <- ifelse('r.squared' %in% names(rfl[[i]]),
                          round(rfl[[i]]$r.squared, 3),
                          round(1-rfl[[i]]$prediction.error, 3))
      return(df)
    }))
    return(df_all)
  }
  ll <- lapply(rf_list, parse_rfl)
  if(multilevel){
     df_plot <- do.call(rbind, lapply(1:length(ll), function(i){
       df_imp <- ll[[i]]
       df_imp[["col_var"]] <- cols[i]
      return(df_imp)
  }) )
  } else{
    df_plot <- ll[[1]]
    df_plot[["col_var"]] <- NA
  }

  df_plot <- df_plot %>%
    dplyr::mutate(rf_name = sprintf("%s_%s", col_var, row_var)) %>%
    group_by(col_var, row_var) %>%
    arrange(importance_value) %>%
    mutate(order = seq_along(importance_value)) %>%
    ungroup()


  # if(!all(is.null(df_all$rf_name))){
  #   df_all$rf_name <- factor(df_all$rf_name, levels=names(rf_list), ordered=TRUE)
  # }

  if(!is.null(names_function)){
    df_plot$var_name <- names_function(df_plot$var_name)
  }

  # ordering <- df_plot %>%
  #   group_by(rf_name) %>%
  #   slice_max(order_by = importance_value, n=n_top) %>%
  #   droplevels() %>%
  #   do(data.frame(al=levels(reorder(interaction(.$rf_name, .$var_name, drop=TRUE),
  #                                   .$importance_value)))) %>%
  #   pull(al)

  # varname <- gsub("^.*\\.", "", ordering)


  # R^2 labels and x coordinate of placement:
  label_importance <- df_plot %>%
    group_by(row_var, col_var, rf_name) %>%
    summarize(x_placement = median(c(min(importance_value),
                                     max(importance_value))))

  r2_labels <- df_plot %>%
    distinct(row_var, col_var, rf_name, r2_oob) %>%
    full_join(label_importance)

  # order_by <- cols[1]
  p <- df_plot %>%
    group_by(col_var, row_var) %>%
    arrange(col_var, row_var, importance_value) %>%
    ungroup() %>%
    mutate(var_ordered = paste(var_name, row_var, col_var, sep = "_")) %>%
    mutate(var_ordered = factor(var_ordered, levels = unique(.data[["var_ordered"]]))) %>%
    # drops unused levels in a way that actually works here:
    # filter(!is.na(var)) %>%
    ggplot(aes(x=round(importance_value),
               y = var_ordered))+# factor(var_ordered, labels = gsub("_.*", "", var_ordered)))) +
    ggh4x::facet_grid2(row_var ~
                         factor(col_var, levels = c(cols[1], setdiff(cols, cols[1]))),
                       scales = "free", independent = "y") +
    # geom_line(aes(color=var_name, group=var_name)) +
    geom_dotplot(binaxis = 'y',
                 aes(fill=var_name, color=var_name)) +
    scale_y_discrete(labels = function(x) gsub("_0.*", "", x)) + # Remove suffix from labels

    ggplot2::geom_text(data=r2_labels,
                       ggplot2::aes(label=sprintf("R2 =\n %s",
                                                  as.character(r2_oob)),
                                    x=x_placement),
                       # size = 14*text_multiplier,
                       y=2) +

    # facet_grid(row_var~col_var, scales="free") +
    ggplot2::labs(y='Variable',
                  x='Importance Value',
                  title = sprintf("Top %s variables by model", n_top))

  return(ggtheme(p, scale_x_date=FALSE, text_multiplier = text_multiplier) +
           ggplot2::theme(legend.position='none') +
           ggplot2::theme(axis.text.y = element_text(hjust=1))
  )
}
