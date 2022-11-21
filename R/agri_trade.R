#' Agricultural trade multiplex network
#'
#' An undirected multiplex network containing trade volumes for 13
#' highly traded agricultural products for the year 2010, collected by
#' the Food and Agriculture Organization of the United Nations (FAO).
#' The original data set can be downloaded from Manlio DeDomenico's website.
#' Array entries are in units of tonnes (metric tons) of bilateral trade
#' of a given agricultural product. For further documentation and
#' product definitions see
#' \url{https://www.fao.org/faostat/en/#definitions/}.
#'
#' @usage data(agri_trade)
#'
#' @format An array of dimension \eqn{145 \times 145 \times 13}.
#'
#' @source \url{https://manliodedomenico.com/}; \url{https://www.fao.org/faostat/en/#data/}
#'
#' @references
#' DeDomenico et al. (2015) \href{https://www.nature.com/articles/ncomms7864/}{Nature Communications}
"agri_trade"
