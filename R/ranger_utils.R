#' @name calculate_pd
#' Calculate partial dependencies
#' (based on function with same name in ebirdstwf, but doesn't require ebirdst
#' fit_predict output)
#' @export
calculate_pd <- function(predictor, model, data, x_res, n,
                         quantile_grid = FALSE,
                         x_range = c(-Inf, Inf),
                         trim_right = FALSE){
  # values of focal predictor
  v <- data[[predictor]]

  keeps <- v >= x_range[1] & v <= x_range[2]
  v <- v[keeps]

  if (trim_right && !(predictor == "solar_noon_diff")) {
    # only use values below 95% quantile of non-zero values
    v_nozero <- v[v > 0]
    if (length(v_nozero) >= 30) {
      v <- v[v <= stats::quantile(v_nozero, probs = 0.95, na.rm = TRUE)]
    } else {
      # if there are < 30 non-zero values, don't produce a pd
      return(NULL)
    }
  }
  if (quantile_grid) {
    # define quantile prediction grid
    x_grid <- stats::quantile(v,
                              probs = seq(from = 0, to = 1, length = x_res),
                              na.rm = TRUE)
  } else {
    # define evenly spaced prediction grid
    rng <- range(v, na.rm = TRUE)
    x_grid <- seq(rng[1], rng[2], length.out = x_res)
  }
  # remove duplicates
  dupes <- duplicated(signif(x_grid, 8))
  x_grid <- x_grid[!dupes]
  x_grid <- unname(unique(x_grid))

  # if x-grid only has a single unique value, don't produce a pd
  if (length(x_grid) == 1) {
    return(NULL)
  }

  grid <- data.frame(predictor = predictor, x = x_grid)
  names(grid) <- c("predictor", predictor)
  # subsample training data
  n <- min(n, nrow(data))
  s <- sample(seq_len(nrow(data)), size = n, replace = FALSE)
  data <- data[s, ]
  # drop focal predictor from data
  data <- data[names(data) != predictor]

  # Handle derived variables:
  if(predictor == "solar_noon_diff"){
    # recalculate abs_solar_noon_diff
    data$abs_solar_noon_diff <- NULL
    grid$abs_solar_noon_diff <- abs(grid$solar_noon_diff)
  }
  grid <- merge(grid, data, all = TRUE)

  # summarize
  pd <- grid[, c("predictor", predictor)]
  names(pd) <- c("predictor", "value")
  # predict
  # TODO: calling .$predictions might be specific to ranger.
  # Maybe use the broom package to access model predictions?
  pd$response <- stats::predict(model, grid)$predictions
  probs <- seq(0.05, 0.95, by = 0.05)

  pds <- dplyr::group_by(pd, .data$predictor, .data$value) %>%
    dplyr::reframe(prob = seq(0.05, 0.95, by = 0.05),
                   quantile = quantile(.data$response, prob),
                   response_mean = mean(.data$response),
                   response_median = median(.data$response)) %>%
    tidyr::pivot_wider(names_from = prob,
                       names_prefix = "q_",
                       values_from = quantile)

  # assign weights based on number of checklists within each interval
  pds <- dplyr::arrange(pds, .data$value)

  v_binned <- cut(v, breaks = pds$value, include.lowest = TRUE, right = FALSE)
  pds$n <- c(as.integer(table(v_binned)), NA_integer_)
  return(pds)
}

#' @name get_var_imp_rf
#' @description Get variable importances from a `ranger` model object as a d
#' @param rf `ranger` model object
#' @return  a data frame with names `var_name`, `importance_type`, and `importance_value`
#' @export
get_var_imp_rf <- function(rf){
  imp_type <- rf$importance.mode
  imp_v    <- rf$variable.importance
  df       <- data.frame(var_name=names(imp_v), importance_value=as.numeric(imp_v), importance_type=imp_type)
  return(df[order(df$importance_value, decreasing = TRUE),])
}


#' Get tree depth for each tree
#' a vector of length num.trees that represents the number of nodes in the longest continuous branch, inclusive of the first and last
#' #' @export
tree_depth <- function(rf){
  sapply(1:rf$num.trees, function(t){
    ti <- ranger::treeInfo(rf, tree=t)
    # from tree info, get the largest node ID and trace back to count how many
    # nodes are on its branch
    max_node <- max(ti$nodeID)
    current_node <- max_node
    n_nodes <- 1
    while(current_node > 0){
      # it can only be the right or the left, never both:
      l <- with(ti, nodeID[which(leftChild==current_node)])
      r <- with(ti, nodeID[which(rightChild==current_node)])
      current_node <- ifelse(length(l) > 0, l, r)
      n_nodes <- n_nodes + 1
    }
    return(n_nodes)
  })
}

#' Get the number of out-of-bag trees for each data point
#' @param rf ranger random forest object (with inbag.counts)
#'
#' @export
get_ntrees_oob <- function(rf){
  # vector of size nobs. for each obs, the number of trees for which that obs was out-of-bag
  # in other words, the number of trees that contribute to its prediction value
  if(is.null(rf$inbag_counts)){
    return(NULL)
  } else{
    # nobs x ntrees matrix
    inbag_mat <- do.call(cbind,rf$inbag_counts)
    rowSums(inbag_mat == 0)
  }
}

#' Get predicted classifications from weighted mean of class assignment probability
#' (instead of maximum probability across the classes)
#' @param pred_mat matrix of classification probabilities (n_obs x n_classes)
#' @export
preds_from_weighted_mean_probs <- function(pred_mat){
  classes <- as.numeric(colnames(pred_mat))
  stopifnot(!is.null(classes))
  # for each row of the matrix, calculate a weighted mean of year using the row values as weights
  sapply(1:nrow(pred_mat), function(i){
    w <- pred_mat[i,]
    weighted.mean(classes, w)
  })

}

#' @export
preds_from_modal_probs <- function(pred_mat, rf=NULL){
  if(is.null(colnames(pred_mat))){
    class_vals <- rf$forest$class.values
    if(is.null(class_vals) | is.null(dim(pred_mat))){
      stop("Error: Not a probability forest?")
    }
    # get them from the model object
    classes <- rf$forest$class.values
    colnames(pred_mat) <- classes
  }
  if(is.null(colnames(pred_mat)) & is.null(rf)){
    stop("Need values for the col classes")
  }
  colnames(pred_mat)[apply(pred_mat, 1, which.max)] %>% as.numeric()
}

