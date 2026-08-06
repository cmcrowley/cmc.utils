compare_pois_nbinom <- function(dispersion, lambda){
  df <- data.frame(nbinoms=rnbinom(n=1000, size=dispersion, mu=lambda),
                   npoiss = rpois(n=1000, lambda=lambda)) %>%
    tidyr::pivot_longer(cols=tidyselect::everything(),
                        names_to='dist')

  p <- ggplot2::ggplot(data=df,
                       ggplot2::aes(x=value)) +
    ggplot2::geom_histogram(ggplot2::aes(fill=dist)) +
    ggplot2::facet_wrap(~dist)

  return(p)
}
