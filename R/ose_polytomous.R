# observed score equating

para <- list(j1 = c(1.197, -1.906, 0.103, 1.713),
              j2 = c(1.029, -2.094, -0.208, 2.020),
              j3 = c(1.627, -2.335, -1.481, -0.803, -0.197, 0.551, 1.670)
              )

para <- list(j1 = c(1.197, -1.906, 0.103, 1.713, 2),
             j2 = c(1.029, -2.094, -0.208, 2.020, 3),
             j3 = c(1.627, -2.335, -1.481, -0.803)
)


# Define pgrm function
pgrm <- function(theta, a, b, k){
  K <- length(b)
  if(k == K){
    p1 <- 0
    p0 <- 1/(1+exp(-a*(theta - b[k-1])))
  } else if(k == 1){
    p1 <- 1/(1+exp(-a*(theta - b[k])))
    p0 <- 1
  } else if(k > 1 && k < K){
    p1 <- 1/(1+exp(-a*(theta - b[k])))
    p0 <- 1/(1+exp(-a*(theta - b[k-1])))
  } else {
    p1 <- 0
    p0 <- 1
  }
  p0 - p1 # probability observed category k response
}

# Obserbed score distribution function
obs_sub <- function(para, theta){
  J <- length(para) # the number of items
  N <- length(unlist(para)) - J # maximum total score
  prb <- matrix(rep(0, N))
  K <- unlist(purrr::map(para, function(x){length(x) - 1}))
  # item 1, special routine.
  for(k in 1:K[1]){
    prb[k] <- pgrm(theta, para[[1]][1], para[[1]][-1], k)
  }
  # after item 2
  tmp <- 0 # initialize
  for(j in 2:J){
    tmp <- tmp + K[j-1] # max raw score that can be scored until item j
    f0 <- prb[(j-1):tmp]
    for(k in 1:K[j]){
      p <- pgrm(theta, para[[j]][1], para[[j]][-1], k)
      cat((j+k-1):(j+k+tmp-1), "\n", f0, "\n", p, "\n")
      prb[(j+k-1):(k+tmp)] <- prb[(j+k-1):(k+tmp)] + f0 * p
    }
  }
  prb
}

obs_equating <- function(para, theta){
  apply(matrix(theta), 1, obs_sub, para = para)
}

theta <- rnorm(1)

obs_equating(para, theta = 1)
