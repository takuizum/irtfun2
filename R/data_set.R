# data set

#'Simulated item response data for 1PLM
#'
#'This data set contains 30 item responses of 3,000 subjects
#'@format a dataset with 31 cols, first col is ID line.
"sim_data_1"

#'Simulated item response data for 2PLM
#'
#'This data set contains 30 item responses of 3,000 subjects
#'@format a dataset with 31 cols, first col is ID line.
"sim_data_2"

#'Simulated item response data for 3PLM
#'
#'This data set contains 30 item responses of 3,000 subjects
#'@format a dataset with 31 cols, first col is ID line.
"sim_data_3"

#'Simulated multi group item response data for 2PLM
#'
#'This data set contains 30 item responses of 9,000 subjects of 3 diffelent level groups.
#'@format a dataset with 31 cols, first col is ID line and second col is group ID line.
"sim_data_4"

#'Estimated 2PLM item parameter.
#'
#'A list of estimated result from \code{\link{estip}}.
#'@format See \code{\link{estip}}.
#'@examples
#'sim_param$para
#'sim_param$SE
"sim_param"

#'true 2PLM item parameter.
#'
#'A list of true parameter .
#'@format See \code{\link{sim_param}}.
#'@examples
#' set.seed(0204)
#' theta <- rnorm(3000)
#' a <- rlnorm(30, sdlog = 0.25)
#' b <- rnorm(30)
#' list(theta,a,b)
"true_param"

#'Estimated theta parameter under 2PLM.
#'
#'A list of estimated result from \code{\link{estheta}}.
#'@format See \code{\link{estheta}}.
#'@examples
#'sim_eap$res
"sim_eap"

