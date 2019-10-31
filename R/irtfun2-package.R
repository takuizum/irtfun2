#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom stats dnorm
#' @importFrom stats na.omit
#' @importFrom stats optimise
#' @importFrom stats optim
#' @importFrom stats runif
#' @importFrom stats cor
#' @importFrom stats qnorm
#' @importFrom stats dlnorm
#' @importFrom stats dbeta
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 stat_function
#' @importFrom ggplot2 labs
#' @importFrom latex2exp TeX
#' @useDynLib irtfun2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

.onUnload = function(libpath) {
  library.dynam.unload("irtfun2", libpath)
}

