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