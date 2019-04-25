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
      # cat((j+k-1):(j+k+tmp-1), "\n", f0, "\n", p, "\n")
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



# # -----
# # Define GRM & GPCM function----
# # thetaはひとつずつ流すので，applyは考えなくて良い
# pgrm <- function(theta, a, b, k, D){
#   K <- length(b)
#   if(k == K){
#     p1 <- 0
#     p0 <- 1/(1+exp(-D*a*(theta - b[k])))
#   } else if(k == 0){
#     p1 <- 1/(1+exp(-D*a*(theta - b[1])))
#     p0 <- 1
#   } else if(k > 0 && k < K){
#     p1 <- 1/(1+exp(-D*a*(theta - b[k+1])))
#     p0 <- 1/(1+exp(-D*a*(theta - b[k])))
#   } else {
#     p1 <- 0
#     p0 <- 1
#   }
#   p0 - p1 # probability observed category k response
# }
#
# # Define GPCM function
# pgpcm <- function(theta, a, b, k, D){
#   K <- length(b)
#   G <- rep(1,K+1)
#   for(v in 1:K) G[v+1] <- exp(sum(D*a*(theta-b[1:v])))
#   p <- G[k+1]/sum(G)
#   p
# }
#
#
# # Obserbed score distribution function
# obs_sub_grm <- function(para, theta, model, D){
#   J <- length(para) # the number of items
#   N <- length(unlist(para)) - J # maximum total score
#   prb <- matrix(rep(0, N))
#   K <- unlist(purrr::map(para, function(x){length(x)-1})) # K is n of b parameter(max category score, not catogory)
#   # item 1, special routine.
#   for(k in 0:K[1]){
#     prb[k+1] <- pgrm(theta, para[[1]][1], para[[1]][-1], k = k, D = D)
#   }
#   # after item 2
#   tmp <- 0 # initialize
#   for(j in 2:J){
#     tmp <- tmp + K[j-1] # max raw score that can be scored until item j
#
#     f0 <- prb[1:(tmp+1)]
#     for(k in 0:K[j]){
#       p <- pgrm(theta, para[[j]][1], para[[j]][-1], k = k, D = D)
#       minr <- 0
#       maxr <- j+k+tmp-1
#       cat(minr:maxr, "\n", f0, "\n", p, "\n")
#       prb[minr:maxr] <- prb[minr:maxr] + f0 * p
#     }
#   }
#   prb
# }
#
# # GPCM version
# obs_sub_gpcm <- function(para, theta, model, D){
#   J <- length(para) # the number of items
#   N <- length(unlist(para)) - J # maximum total score
#   prb <- matrix(rep(0, N))
#   K <- unlist(purrr::map(para, function(x){length(x)})) # K is n of b parametet + 1
#   # item 1, special routine.
#   for(k in 0:K[1]){
#     prb[k] <- pgpcm(theta, para[[1]][1], para[[1]][-1], k = k, D = D)
#   }
#   # after item 2
#   tmp <- 0 # initialize
#   for(j in 2:J){
#     tmp <- tmp + K[j-1] # max raw score that can be scored until item j
#     f0 <- prb[(j-1):tmp]
#     for(k in 0:K[j]){
#       p <- pgpcm(theta, para[[j]][1], para[[j]][-1], k = k, D = D)
#       minr <- j+k-1
#       maxr <- j+k+tmp-1
#       cat(minr:maxr, "\n", f0, "\n", p, "\n")
#       prb[minr:maxr] <- prb[minr:maxr] + f0 * p
#     }
#   }
#   prb
# }
#
# #---- wrapper function----
# obs_equating <- function(para, theta, D = 1.0,  model = "GRM"){
#   if(model == "GRM"){
#     apply(matrix(theta), 1, obs_sub_grm, para = para, model = model, D = D)
#   } else if(model == "GPCM"){
#     apply(matrix(theta), 1, obs_sub_gpcm, para = para, model = model, D = D)
#   }
# }
#
# # validation
# para <- list(j1 = c(1.197, -1.906),
#              j2 = c(1.029, -2.094),
#              j3 = c(1.627, -2.335)
# )
#
# theta <- rnorm(2)
# obs_equating(para, theta = 0, model = "GRM")
#
