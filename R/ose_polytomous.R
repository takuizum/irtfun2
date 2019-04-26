# observed score equating

para <- list(j1 = c(1.197, -1.906, 0.103, 1.713),
              j2 = c(1.029, -2.094, -0.208, 2.020),
              j3 = c(1.627, -2.335, -1.481, -0.803, -0.197, 0.551, 1.670)
              )

para <- list(j1 = c(1.197, -1.906, 0.103, 1.713),
             j2 = c(1.029, -2.094, -0.208, 2.020),
             j3 = c(1.627, -2.335, -1.481)
)
para <- list(j1 = c(1.197, -1.906),
             j2 = c(1.029, -2.094),
             j3 = c(1.627, -2.335),
             j4 = c(1.223,  2.456),
             j5 = c(1.223,  2.456)
)

# Define pgrm function
pgrm <- function(theta, a, b, k){
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

# Obserbed score distribution function
obs_sub <- function(para, theta, min_category){
  J <- length(para) # the number of items
  N <- length(unlist(para)) # maximum total score()
  K <- unlist(purrr::map(para, function(x){length(x)}))
  # item 1, special routine.
  prb <- rep(0, K[1]) # too short vector still at tis point
  for(k in 1:K[1]) prb[k] <- pgrm(theta, para[[1]][1], para[[1]][-1], k)
  # after item 2
  for(j in 2:J){
    maxrr <- j:sum(K[1:j]) - j + 1
    f0 <- prb
    prb <- numeric(length(maxrr)) # refresh
    for(k in 1:K[j]){
      p <- pgrm(theta, para[[j]][1], para[[j]][-1], k)
      # cat(tmp+k, "\n", f0, "\n", p, "\n")
      prb <- prb + c(rep(0,k-1), c(f0 * p), rep(0,K[j]-k))
      # cat(prb, "!!!\n", c(rep(0,k-1), c(f0 * p), rep(0,K[j]-k)), "|||\n")
    }
  }
  prb
}

obs_equating <- function(para, theta, weight = NULL, min_category = 1){
  if(is.null(weight)) weight <- rep(1, length(theta)) # uniform dist
  prb_matrix <- apply(matrix(theta), 1, obs_sub, para = para, min_category = min_category)
  prb <- prb_matrix %*% weight
  J <- length(para) # the number of items
  K <- unlist(purrr::map(para, function(x){length(x)}))
  maxrr <- J:sum(K) - J + 1
  tibble::tibble(SCORE = maxrr+J*min_category-1, prb = as.vector(prb))
}


# validation
theta <- seq(-4, 4, length.out = 31)
weight <- dnorm(theta)/sum(dnorm(theta))
test <- obs_equating(para, theta = theta, weight = weight, min_category = 1)

