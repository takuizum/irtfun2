# observed score equating

# packages to need to install (without library)
# installed.packages(c("tibble", "purrr"))

# Define pgrm function
# graded response model
pgrm <- function(theta, a, b, k, D){
  K <- length(b) + 1
  if(k == K){
    p1 <- 0
    p0 <- 1/(1+exp(-a*(theta - b[k-1])))
  } else if(k == 1){
    p1 <- 1/(1+exp(-a*(theta - b[k])))
    p0 <- 1
  } else if(k > 1 && k < K){
    p1 <- 1/(1+exp(-a*(theta - b[k])))
    p0 <- 1/(1+exp(-a*(theta - b[k-1])))
  } else if(is.na(k)){
    p1 <- 0
    p0 <- 1
  } else {
    p1 <- 0
    p0 <- 0
  }
  p0 - p1 # probability observed category k response
}

# generalized partial credit model
pgpcm <- function(theta, a, b, k, D){
  K <- length(b)
  G <- rep(1,K+1)
  for(v in 1:K) G[v+1] <- exp(sum(D*a*(theta-b[1:v])))
  p <- G[k]/sum(G) # k is start from 1 not 0
  p
}

# Obserbed score distribution function (sub routine assigned for apply function)
obs_sub <- function(para, theta, min_category, model, D){
  J <- length(para) # the number of items
  K <- unlist(purrr::map(para, function(x){length(x)}))
  # model selection
  if(model == "GRM"){
    pf <- pgrm
  } else if(model == "GPCM"){
    pf <- pgpcm
  } else if(model == "NRM"){
    stop("Error in model definition. 'GRM' or 'GPCM' is only available now.")
    # pf <- NULL
  } else {
    stop("Error in model definition. 'GRM' or 'GPCM' is only available now.")
  }
  # item 1, special routine.
  prb <- rep(0, K[1]) # too short vector still at tis point
  for(k in 1:K[1]) prb[k] <- pf(theta, para[[1]][1], para[[1]][-1], k, D)
  # after item 2
  for(j in 2:J){
    maxrr <- j:sum(K[1:j]) - j + 1
    f0 <- prb
    prb <- numeric(length(maxrr)) # refresh
    for(k in 1:K[j]){
      p <- pf(theta, para[[j]][1], para[[j]][-1], k, D)
      # cat(f0, "\n", p, "\n") # for debugging
      prb <- prb + c(rep(0,k-1), c(f0 * p), rep(0,K[j]-k))
      # cat(prb, "!!!\n", c(rep(0,k-1), c(f0 * p), rep(0,K[j]-k)), "|||\n") # for debugging
    }
  }
  prb
}

obs_equating <- function(para, theta, weight = NULL, model = "GRM", min_category = 1, D = 1.0){
  # theta's weight check
  if(is.null(weight)) weight <- rep(1, length(theta))/length(theta) # uniform dist
  # parameterfile check
  # if(!is.list(para)){}
  prb_matrix <- apply(matrix(theta), 1, obs_sub, para = para, min_category = min_category, D = D, model = model)
  prb <- prb_matrix %*% weight
  J <- length(para) # the number of items
  K <- unlist(purrr::map(para, function(x){length(x)}))
  maxrr <- J:sum(K) - J + 1
  tibble::tibble(SCORE = maxrr+J*min_category-1, prb = as.vector(prb))
}


# validation----

# # GRM parameters
# para <- list(j1 = c(1.197, -1.906, 0.103, 1.713),
#              j2 = c(1.029, -2.094, -0.208, 2.020),
#              j3 = c(1.627, -2.335, -1.481, -0.803, -0.197, 0.551, 1.670)
# )
#
# para <- list(j1 = c(1.197, -1.906, 0.103, 1.713),
#              j2 = c(1.029, -2.094, -0.208, 2.020),
#              j3 = c(1.627, -2.335, -1.481)
# )
# dichotomous parameter
para <- list(j1 = c(1.197, -1.906),
             j2 = c(1.029, -2.094),
             j3 = c(1.627, -2.335),
             j4 = c(1.223,  2.456),
             j5 = c(1.223,  2.456)
)

# simulation data
theta <- seq(-4, 4, length.out = 31)
# weight <- dnorm(theta)/sum(dnorm(theta))

# polytomous observed score distribution
obs_equating(para, theta = theta, min_category = 0, model = "GRM")

# dichotomous observed score distribution
obs1 <- obscore_dist(theta, c(1.197, 1.029, 1.627, 1.223, 1.223), c(-1.906, -2.094, -2.335, 2.456, 2.456), c(0,0,0,0,0), 1.0)
table(obs1)/length(obs1)
