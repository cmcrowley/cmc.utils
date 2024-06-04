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

