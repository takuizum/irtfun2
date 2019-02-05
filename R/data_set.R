# data set

#'Simulated item response data for 1PLM
#'
#'This data set contains 30 item responses of 3,000 subjects
#'@format a dataset with 31 cols, first col is ID line.
#'@examples
#' set.seed(0204)
#' # single test data: subjects=3,000 and item=30, 1PLM
#' sim_data_1 <- irtfun2::sim_gen(theta=rnorm(3000), b=rnorm(30), a=rep(1,30))
"sim_data_1"

#'Simulated item response data for 2PLM
#'
#'This data set contains 30 item responses of 3,000 subjects
#'@format a dataset with 31 cols, first col is ID line.
#'@examples
#'# single test data: subjects=3,000 and item=30, 2PLM
#' set.seed(0204)
#' theta <- rnorm(3000)
#' a <- rlnorm(30, sdlog = 0.25)
#' b <- rnorm(30)
#' sim_data_2 <- irtfun2::sim_gen(theta=theta, b=b, a=a)
#'
"sim_data_2"

#'Simulated item response data for 3PLM
#'
#'This data set contains 30 item responses of 3,000 subjects
#'@format a dataset with 31 cols, first col is ID line.
#'@examples
#'set.seed(0204)
#'# single test data: subjects=3,000 and item=30, 3PLM
#'# data.frame(x=c(0:1)) %>%
#'#   ggplot(aes(x=x))+stat_function(fun=dbeta, args = list(shape1=1, shape2=6))
#' sim_data_3 <- irtfun2::sim_gen(theta=rnorm(3000),
#'   b=rnorm(30), a=rlnorm(30, sdlog=0.25), c=rbeta(30, shape1=1, shape2=6))
#'
"sim_data_3"

#'Simulated multi group item response data for 2PLM
#'
#'This data set contains 30 item responses of 9,000 subjects of 3 diffelent level groups.
#'@format a dataset with 32 cols, first col is ID line and second col is group ID line.
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
#' #list(theta,a,b)
"true_param"

#'Estimated theta parameter under 2PLM.
#'
#'A list of estimated result from \code{\link{estheta}}.
#'@format See \code{\link{estheta}}.
#'@examples
#'head(sim_eap$res)
"sim_eap"

#'Simulated item response data for GIRT model
#'
#'This data set contains 30 item responses of 3,000 subjects for GIRT model
#'@format a dataset with 31 cols, first col is ID line.See \code{\link{estGip}}.
#'@examples
#'# GIRT estimated parameter
#'set.seed(0204)
#'theta <- rnorm(3000)
#'phi <- rinvchi(3000, max = 2)
#'a <- rlnorm(30, sdlog = 0.25)
#'b <- rnorm(30)
#'sim_dat_girt <- sim_gen(theta=theta, phi=phi, a=a, b=b)
"sim_dat_girt"

#'true 2PLM item parameter of 'sim_dat_st'.
#'
#'A list of true parameter .
#'@format See \code{\link{sim_dat_st}}.
"true_param_st"

#'Simulated multi group item response data for scaling test design.
#'
#'This data set contains 70 item responses of 2,500 subjects for scaling test design.
#'@format a dataset with 72 cols, first col is ID line and second col is group ID line.
"sim_dat_st"
