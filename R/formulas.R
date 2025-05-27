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

#' Remove terms from list that are constant in the data
#' @name pare_formula
#' @param dat data frame to be modeled
#' @param predictors character vector of predictor names
#' @export
pare_formula <- function(dat, predictors){
  get_constants <- function(dat){
    constant <- sapply(names(dat), function(var){
      length(unique(dat[[var]])) == 1
    })
    names(dat)[constant]
  }
  exclusions <- c(
    # vars that have already been excluded for whatever reason
    setdiff(predictors, names(dat)),
    # plus any features that are constant in the data
    get_constants(dat)
  )
  setdiff(predictors, exclusions)
}
