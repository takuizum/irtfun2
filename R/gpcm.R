# GPCM

gpcm_sub <- function(theta, a, b, D, k){
  K <- length(b)
  G <- rep(1,K+1)
  for(v in 1:K) G[v+1] <- exp(sum(D*a*(theta-b[1:v])))
  p <- G[k+1]/sum(G)
  p
}

#'Graded Item Response Model, Generalized Percial Credit Model
#'
#'@inheritParams ptheta
#'@param b a vector of transition parameter.
#'@param k a number of category.
#'@examples
#'tp <- c(-1.5, 0, 1)
#'ggplot2::ggplot(data = data.frame(x=c(-4:4)),
#'                ggplot2::aes(x=x))+
#'  ggplot2::stat_function(fun=gpcm,
#'    args = list(a=1.5, b=tp, D=1.702, k=0),
#'    ggplot2::aes(colour="category0"))+
#'  ggplot2::stat_function(fun=gpcm,
#'    args = list(a=1.5, b=tp, D=1.702, k=1),
#'    ggplot2::aes(colour="category1"))+
#'  ggplot2::stat_function(fun=gpcm,
#'    args = list(a=1.5, b=tp, D=1.702, k=2),
#'    ggplot2::aes(colour="category2"))+
#'  ggplot2::stat_function(fun=gpcm,
#'    args = list(a=1.5, b=tp, D=1.702, k=3),
#'    ggplot2::aes(colour="category3"))+
#'  ggplot2::labs(x=latex2exp::TeX("$\\theta$"),
#'                y=latex2exp::TeX("$P(\\theta)$"),
#'                colour="Category")
#'@export
gpcm <- function(theta, a, b, D, k){
  apply(as.matrix(theta), 1, gpcm_sub, a=a,b=b,D=D,k=k)
}


# gpcm function (incomplete)
#-------------------------------------------------------------------------------------------------------------


# データセットは因子型のデータフレーム，もしくは整数型のdata.frameを想定(tibbleに対応できればなお良い)
library(irtoys)
library(tidyverse)
library(magrittr)

Science[c(1,3,4,7)] %>% purrr::map_df(as.integer)

X <- x %>% switch("f" = x %>% purrr::map_df(as.integer), # if factor type
                  "i" = x)

X <- Science[c(1,3,4,7)] %>% purrr::map_df(as.integer) %>% purrr::map_df(function(x){x-1}) %>% purrr::map_df(as.integer)
X %<>% dplyr::bind_rows(X) # %>% dplyr::mutate(id = 1:nrow(.)) # bind row (add more data lines)

N <- nrow(X)
J <- ncol(X)
Item <- colnames(X)

# list up all cat and min and max cat
cat_item <- X %>% purrr::map(unique) %>% purrr::map(sort) # %>% purrr::map(function(x){x[-1]})
cat_item_n <- X %>% purrr::map(dplyr::n_distinct) # %>% purrr::map(~ .-1)
max_cat_item <- cat_item %>% purrr::map(max, na.rm = T)
max_cat_all <- max_cat_item %>% unlist() %>% max()
min_cat_item <- cat_item %>% purrr::map(min, na.rm = T)
min_cat_all <- min_cat_item %>% unlist() %>% min()

# add group id column
set.seed(0204)
X %<>% dplyr::mutate(group = as.integer(sample(c(1,2), nrow(X), replace = T)))
# X %<>% dplyr::mutate(group = as.integer(rep(1), nrow(X)))

# add count column
X %<>% dplyr::mutate(count = 1) %>% dplyr::arrange(group) # and arrange

# unique response patterns and them count
Xl <- X %>% group_by_if(is.integer) %>% dplyr::summarise(count = sum(count)) %>% dplyr::arrange(group)

# n of sunjects in each groups
group_N <- X %>% dplyr::select(group, count) %>% group_by(group) %>% dplyr::summarise(count = sum(count))
ng <- nrow(group_N)

# design matrix for each groups
design_each_g <- X %>% dplyr::select(-count) %>% group_by(group) %>% dplyr::summarise_all(sum, na.rm = T) %>% tibble::column_to_rownames("group")
design_each_g[design_each_g != 0] <- 1

# design matrix for all response pattern
design_all <- Xl %>% dplyr::select(-count, -group)
design_all[is.na(design_all)] <- NA


# intial value
beta0 <- matrix(0, nrow = J, ncol = length(1:max_cat_all))
k0 <- matrix(0, nrow = J, ncol = length(1:max_cat_all))
b0 <- numeric(J)
a0 <- rep(1, J)
cat_count <- X %>% purrr::map(~ as.vector(table(.)))
cat_count <- X %>% purrr::map(table)

for(j in 1:J){
  cat <- dimnames(cat_count[[j]])[[1]][-1] %>% as.integer()
  cat_j <- cat_count[[j]][-1] %>% as.integer()
  total_x <- rowSums(X %>% dplyr::select(-count, -group))
  a0[j] <- cor(X[,j], total_x)
  prob <- cat_j / (group_N[design_each_g[,j] == 1,] %>% dplyr::select(count) %>% sum())
  beta0[j, cat] <- -log(prob / (1 - prob))
  b0[j] <- mean(beta0[j, cat])
  k0[j, cat] <- beta0[j, cat] - b0[j]
}

rm(prob)

rl <- Xl %>% dplyr::arrange(group)
rl <- rl$count
group_id <- Xl %>% dplyr::arrange(group)
group_id <- group_id$group
Xl_tbl <- Xl # save tibble data as
Xl %<>% dplyr::arrange(group) %>% dplyr::select(-count, -group) %>% as.matrix()
L <- nrow(Xl) # n of response
sf <- max_cat_item %>% unlist() # scoring function "T_j"


# max category in each item vector and repeat it  n of subjects
n_cat_vec <- max_cat_item %>% unlist(use.names = F) %>% rep.int(times = L)

create_ujk <- function(response, max_category){ # both of them is vector
  res <- rep(0, max_category)
  res[response] <- 1
  res
}

# tidyr for long type data.frame
Xl_long <- Xl_tbl %>%
  dplyr::arrange(group) %>%
  tibble::rowid_to_column("id") %>%
  dplyr::select(-count) %>%
  tidyr::gather(key = item, value = response, -group, -id) %>%
  dplyr::arrange(id)
# Xl_long %<>% tibble(max_category = n_cat_vec)

# group vector
group_vec <- Xl_tbl %>% dplyr::arrange(group) %$% group

# create ujk
# ujk <- mapply(create_ujk, Xl_long$response, n_cat_vec) %>% as.vector()
# ujk <- tibble::tibble(u = as.integer(ujk))

# prior of theta
M <- 31 # of node
min_th <- -4
max_th <- 4
xm <- seq(min_th, max_th, length.out = M)
wm <- dnorm(xm)/sum(dnorm(xm)) %>% rep.int(times = ng) %>% matrix(nrow = M, ncol = ng, byrow = T)

# EM algorithm starts

# E step
#----
# X_long <- data_frame(group = )

# logit and probability
Pgjkm <- array(0, dim = c(ng, J, max_cat_all, M))
Zgjkm <- array(0, dim = c(ng, J, max_cat_all, M))
Zgjkm_cum <- array(0, dim = c(ng, J, max_cat_all, M))
Zgjm_bar <- array(0, dim = c(ng, J, M))
Tgjm <- array(0, dim = c(ng, J, M))
Zgjm <- array(0, dim = c(ng, J, M))
for(g in 1:ng){
  item_list <- which(design_each_g[g, ] == 1) # pattern location of group g in data.frame
  for(j in item_list){
    cat_list <- sort(cat_item[[j]])[-1]    #[-cat_item_n[[j]]]
    for(m in 1:M){
      for(k in cat_list){ # calculate Z_jk
        Zgjkm[g, j, k, m] <- 1.702 * a0[j] * (xm[m] - b0[j] + k0[j, k])
        # Zgjkm_cum[g, j, k, m] <- 1.702 * a0[j] * (k * (xm[m] - b0[j]) + sum(k0[j, 1:k]))
      }
      temp <- cumsum(Zgjkm[g,j,,m])
      Zgjkm_cum[g, j,, m] <- temp
      # log_sum_exp
      log_p_const <- (temp - max(temp)) %>% exp() %>% sum() %>% log() %>% magrittr::add(max(temp))
      # const <- exp((temp))
      for(k in cat_list){ # calculate Zjk # ここ，forいらないかも
        Pgjkm[g, j, k, m] <- exp(Zgjkm_cum[g, j, k, m]) / (1 + exp(log_p_const))
      }
      Zgjm[g, j, m] <- Pgjkm[g, j, , m] %*% temp
      Tgjm[g, j, m] <- Pgjkm[g, j, , m] %*% cat_list
    }
  }
}


# likelihood
Lw <- array(0, dim = c(ng, L, M))
Ngjm <- array(0, dim = c(ng, J, M))
rLw <- array(0, dim = c(ng, L, M))
const <- array(0, dim = c(ng, L, M))
for(g in 1:ng){
  pattern_location <- which(group_vec == g) # pattern location of group g in data.frame
  for(l in pattern_location){
    item_list <- which(design_all[l,] != 0) # item location that not NA
    for(m in 1:M){
      p <- wm[m, g]
      for(j in item_list){
        p <- p * Pgjkm[g,j,Xl[l,j],m]
      }
      Lw[g, l, m] <- p
      rLw[g, l, m] <- p * rl[l]
      const[g, l,] <- const[g, l,] + p
    }
  }
}

# expected frequency of subjects that responsed category k in item j
rgjkm <- array(0, dim = c(ng, J, max_cat_all, M))
for(g in 1:ng){
  pattern_location <- which(group_vec == g) # pattern location of group g in data.frame
  for(l in pattern_location){
    item_list <- which(design_all[l,] != 0) # item location that not NA
    for(m in 1:M){
      for(j in item_list){
        rgjkm[g, j, Xl[l, j], m] <- rgjkm[g, j, Xl[l, j], m] + rLw[g, l, m] / const[g, l, m]
      }
    }
  }
}

rjkm <- apply(rgjkm, c(2,3,4), sum)

# # provisional sumple size
Ngm <- matrix(0, ncol = M, nrow = ng)
for(g in 1:ng){
  Ngm[g, ] <- (colSums(rLw[g,group_id == g,]/const[g, group_id == g,]))
}

Nm <- apply(Ngm, 2, sum)


#----
# M step algorithm
#----
M1_gpcm <- function(a0, b0, k0, cate, xm, rr, NN){
  # First M step
  K <- length(cate)
  M <- length(xm)
  Zcum <- Prob <- Zbar <- Tbar <- numeric(K * M)
  Tk <- rep(cate, M)
  NN <- NN %>% rep(rep(K,M))
  rr <- rr %>% as.vector()
  Z <- numeric(K)
  t1 <- 1
  for(m in 1:M){
    Z <- 1.702 * a0 * (xm[m] - b0 + k0)
    temp <- cumsum(Z)
    Zcum[t1:(t1 + K -1)] <- temp
    # log_sum_exp
    log_p_const <- (temp - max(temp)) %>% exp() %>% sum() %>% log() %>% magrittr::add(max(temp))
    Prob[t1:(t1 + K -1)] <- exp(temp) / (1 + exp(log_p_const))
    Zbar[t1:(t1 + K -1)] <- Prob[t1:(t1 + K -1)] %*% temp %>% rep.int(times = K) #
    Tbar[t1:(t1 + K -1)] <- Prob[t1:(t1 + K -1)] %*% cate %>% rep.int(times = K) # cumsum
    t1 <- t1 + K
  }
  a <- a0^(-1) * sum(rr*(Zcum - Zbar))
  b <- a0 * sum(rr*(-Tk + Tbar))
  aa <- a0^(-2) * sum(NN * Prob * (Zcum - Zbar)^2)
  bb <- a0^2 * sum(NN * Prob * (-Tk + Tbar)^2)
  ab <- sum(NN * Prob * (-Tk + Tbar) * (Zcum - Zbar))
  gr <- c(a, b)
  V <-  matrix(c(aa,ab,ab,bb), 2,2)
  # new parameter
  par <- c(a0, b0) + solve(V) %*% gr
  list(par = par, gr = gr, V = V,
       tab = tibble::tibble(k = Tk, P = Prob, r = rr, N = NN, Zcum = Zcum, Tbar = Tbar, Zbar = Zbar) %>% dplyr::arrange(k))
}

M2_gpcm <- function(a0, b0, k0, cat_item, xm, rr, NN){
  # Second M step
  k1 <- k0
  J <- length(b0)
  K <- cat_item %>% purrr::map(.f = function(x){x[-1]}) %>% unlist() %>% length()
  M <- length(xm)
  Prob <- numeric(K * M)
  t1 <- 1
  for(m in 1:M){
    for(j in 1:J){
      KK <- max_cat_item[[j]]#length(cat_item[[j]])
      Z <- 1.702 * a0[j] * (xm[m] - b0[j] + k0[j,])
      temp <- cumsum(Z)
      # log_sum_exp
      log_p_const <- (temp - max(temp)) %>% exp() %>% sum() %>% log() %>% magrittr::add(max(temp))
      Prob[t1:(t1 + KK -1)] <- exp(temp) / (1 + exp(log_p_const)) # categpry 0 to K
      t1 <- t1 + KK
    }
  }
  # contain the probability of categpry 0
  XX <- tibble::tibble(node = c(1:M) %>% rep(rep(sum(unlist(max_cat_item)), M)),
                      item = c(1:J) %>% rep.int(times = (unlist(max_cat_item, use.names = F))) %>% rep.int(times = M),
                      category = cat_item %>% purrr::map(.f = function(x){x[-1]}) %>% unlist(use.names = F) %>% rep.int(times = M),
                      P = Prob)

  r_bar <- apply(rr, c(1,3), sum)
  for(j in 1:J){
    gr_k <- numeric(max_cat_item[[j]])
    V_k <- matrix(0, nrow = max_cat_item[[j]], ncol = max_cat_item[[j]])
    rr_j <- rr[j,,]
    for(k in 1:max_cat_item[[j]]){
      rr_jk <- rr_j[k:max_cat_item[[j]],]
      X1 <- XX %>% dplyr::filter(item == j, category >= k) %>% dplyr::mutate(rr = rr_jk %>% as.vector)
      # gr_k <- sum(a0[j] * (X1$rr - X1$P * r_bar[j, X1$node])) # equation 6.26
      gr_k[k] <- sum(a0[j] * (X1$rr - X1$P * r_bar[j, X1$node])) # equation 6.26
      X1 <- X1 %>% group_by(node) %>% dplyr::summarise(cumP = sum(P))
      # X0 <- XX %>% dplyr::filter(item == j, category >= k) %>% group_by(node) %>% dplyr::summarise(cumP = sum(P))
      # V_k <- sum(NN * a0[j]^2 * X1$cumP * (1 - X0$cumP)) # equation 6.27
      # k1[j,k] <- k0[j,k] + V_k^(-1)*gr_k
      col_list <- 1:max_cat_item[[j]]
      col_list <- col_list[col_list <= k]
      for(l in col_list){
        X0 <- XX %>% dplyr::filter(item == j, category >= l) %>% group_by(node) %>% dplyr::summarise(cumP = sum(P))
        V_k[k,l] <- sum(NN * a0[j]^2 * X1$cumP * (1 - X0$cumP)) # equation 6.27
      }
    }
    # V_k[upper.tri(V_k)] <- V_k[lower.tri(V_k)]
    k1[j,] <- k0[j,] + solve(V_k) %*% gr_k
  }
  list(res = k1, V = V_k, g = gr_k)
}

# intial value
{beta0 <- matrix(0, nrow = J, ncol = length(1:max_cat_all))
k0 <- matrix(0, nrow = J, ncol = length(1:max_cat_all))
b0 <- numeric(J)
# a0 <- rep(1, J)
# a0 <- c(0.8, 0.8, 2.2, 0.7)/1.702
cat_count <- X %>% purrr::map(~ as.vector(table(.)))
cat_count <- X %>% purrr::map(table)

for(j in 1:J){
  cat <- dimnames(cat_count[[j]])[[1]][-1] %>% as.integer()
  cat_j <- cat_count[[j]][-1] %>% as.integer()
  total_x <- rowSums(X %>% dplyr::select(-count, -group))
  R <- cor(X[,j], total_x)
  a0[j] <- R/(1-R^2)
  prob <- cat_j / (group_N[design_each_g[,j] == 1,] %>% dplyr::select(count) %>% sum())
  beta0[j, cat] <- -log(prob / (1 - prob))
  b0[j] <- mean(beta0[j, cat])
  k0[j, cat] <- beta0[j, cat] - b0[j]
}

rm(prob)}

a1 <- numeric(J)
b1 <- numeric(J)
mstep_iter <- TRUE
while(mstep_iter){
  # first block M step
  for(j in 1:J){
    t1 <- M1_gpcm(a0[j], b0[j], k0[j,], cate = cat_item[[j]][-1], xm = xm, rr = rjkm[j,,], NN = Nm)
    a1[j] <- t1$par[1]
    b1[j] <- t1$par[2]
  }
# second block M step
  k1 <- M2_gpcm(a1, b1, k0, cat_item, xm, rjkm, Nm)$res

  # convergence check
  if(all(abs(a0 - a1) < 0.001) && all(abs(b0 - b1) < 0.001) &&  all(abs(k0 - k1) < 0.001)){
    mstep_iter <- FALSE
  } else {
    a0 <- a1
    b0 <- b1
    t0 <- t1
  }
}


# log likelihood function
# N! = gamma(n + 1)


