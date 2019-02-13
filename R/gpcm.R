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


# データセットは因子型のデータフレーム，もしくは整数型のdata.frameを想定(tibbleに対応できればなお良い)
library(irtoys)
library(tidyverse)
library(magrittr)

Science %>% purrr::map_df(as.integer)

X <- x %>% switch("f" = x %>% purrr::map_df(as.integer), # if factor type
                  "i" = x)

X <- Science %>% purrr::map_df(as.integer)#%>% as.matrix()
X %<>% dplyr::bind_rows(X) # %>% dplyr::mutate(id = 1:nrow(.)) # bind row (add more data lines)

N <- nrow(X)
J <- ncol(X)
Item <- colnames(X)

# list up all cat and min and max cat
cat_item <- X %>% purrr::map(unique) %>% purrr::map(sort)
cat_item_n <- X %>% purrr::map(dplyr::n_distinct)
max_cat_item <- cat_item %>% purrr::map(max, na.rm = T)
max_cat_all <- max_cat_item %>% unlist() %>% max()
min_cat_item <- cat_item %>% purrr::map(min, na.rm = T)
min_cat_all <- min_cat_item %>% unlist() %>% min()

# mapply(create_ujk, )

# add group id column
set.seed(0204)
X %<>% dplyr::mutate(group = as.integer(sample(c(1,2), nrow(X), replace = T)))

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
design_all[is.na(design_all)] <- 0

# intial value
beta0 <- matrix(0, nrow = J, ncol = length(min_cat_all:max_cat_all))
k0 <- matrix(0, nrow = J, ncol = length(min_cat_all:max_cat_all))
a0 <- rep(1, J)
cat_count <- X %>% purrr::map(~ as.vector(table(.)))
cat_count <- X %>% purrr::map(table)

for(j in 1:J){
  cat <- dimnames(cat_count[[j]])[[1]] %>% as.integer()
  cat_j <- cat_count[[j]] %>% as.integer()
  prob <- cat_j / (group_N[design_each_g[,j] == 1,] %>% dplyr::select(count) %>% sum())
  beta0[j, cat] <- -log(prob / (1 - prob))
}
b0 <- rowMeans(beta0)
# b0 <- rep(0, J)
k0 <- beta0 - matrix(b0, nrow = 7, ncol = max_cat_all)
k0 %<>% scale()
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
group_vec <- Xl_tbl %>% dplyr::arrange(group) %>% .$group

# create ujk
ujk <- mapply(create_ujk, Xl_long$response, n_cat_vec) %>% as.vector()
ujk <- tibble::tibble(u = as.integer(ujk))
# prior of theta
M <- 31 # of node
min_th <- -4
max_th <- 4
xm <- seq(min_th, max_th, length.out = M)
wm <- dnorm(xm)/sum(dnorm(xm)) %>% rep.int(times = ng) %>% matrix(nrow = M, ncol = ng, byrow = T)

# EM algorithm starts

# E step
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
    cat_list <- sort(cat_item[[j]])
    for(m in 1:M){
      for(k in cat_list){ # calculate Z_jk
        Zgjkm[g, j, k, m] <- a0[j] * (xm[m] - b0[j] + k0[j, k])
      }
      temp <- cumsum(Zgjkm[g, j, , m])
      Zgjkm_cum[g, j, , m] <- temp
      # log_sum_exp
      log_p_const <- (temp - max(temp)) %>% exp %>% sum() %>% log() %>% magrittr::add(max(temp)) %>% exp() %>% magrittr::add(1) %>% log()
      # log_p_const <- log(1 + exp(log_p_const))
      # p_const <- 1 + sum(exp(temp))
      for(k in cat_list){ # calculate Zjk
        # Pgjkm[g, j, k, m] <- (1 + exp(Zgjkm_cum[g, j, k, m])) / p_const
        Pgjkm[g, j, k, m] <- exp(log(1 + exp(Zgjkm_cum[g, j, k, m])) - log_p_const)
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

# # provisional sumple size
# Nf <- numeric(M)
# # expected frequency of subjects that responsed category k in item j
# rjkm <- array(0, dim = c(J, max_cat_all, M))
# Pjkm <- array(0, dim = c(J, max_cat_all, M))
# Zjkm_cum <- array(0, dim = c(J, max_cat_all, M))
# Zjm <- matrix(0, nrow = J, ncol = M)
# Tjm <- matrix(0, nrow = J, ncol = M)
# for(g in 1:ng){
#   Nf <- Nf + (colSums(rLw[g,group_id == g,]/const[g, group_id == g,]))
#   rjkm <- rjkm + rgjkm[g,,,]
#   Pjkm <- Pjkm + Pgjkm[g,,,]
#   Zjkm_cum <- Zjkm_cum + Zgjkm_cum[g,,,]
#   Zjm <- Zjm + Zgjm[g,,]
#   Tjm <- Tjm + Tgjm[g,,]
# }
# Pjkm %<>% divide_by(ng)
# Zjkm_cum %<>% divide_by(ng)
# Zjm %<>% divide_by(ng)
# Tjm %<>% divide_by(ng)
Ngm <- matrix(0, ncol = M, nrow = ng)
for(g in 1:ng){
  Ngm[g, ] <- (colSums(rLw[g,group_id == g,]/const[g, group_id == g,]))
}


# M step algorithm

# jごとにベクトル化して，そのベクトルを受け取った関数でgradientや情報行列を計算すればよい
# group_N$group %>% rep.int(times = rep(sum(unlist(max_cat_item)) * M, ng)) # of group
# c(1:J) %>% rep.int(times = unlist(max_cat_item, use.names = F) * M) # %>% rep.int(times = ng) # of item
# cat_item %>% unlist(use.names = F) %>% rep.int(times = rep(M, length(.))) #%>% rep.int(ng) # of category

# sort
t <- 0
rgjkm_vec <- numeric(length = length(as.vector(rgjkm)))
Ngm_vec <- numeric(length = length(rjkm_vec))
Pgjkm_vec <- numeric(length = length(rjkm_vec))
Zgjkm_cum_vec <- numeric(length = length(rjkm_vec))
Zgjm_vec <- numeric(length = length(rjkm_vec))
Tgjm_vec <- numeric(length = length(rjkm_vec))
for(g in 1:ng){
  for(j in 1:J){
    for(k in 1:max_cat_all){
      for(m in 1:M){
        # cat(t,"\r")
        t <- t + 1
        rgjkm_vec[t] <- rgjkm[g,j,k,m]
        Ngm_vec[t] <- Ngm[g,m]
        Pgjkm_vec[t] <- Pgjkm[g,j,k,m]
        Zgjkm_cum_vec[t] <- Zgjkm_cum[g,j,k,m]
        Zgjm_vec[t] <- Zgjm[g,j,m]
        Tgjm_vec[t] <- Tgjm[g,j,m]
      }
    }
  }
}


X_long <- dplyr::data_frame(group = c(1:ng) %>% rep.int(times = rep(sum(unlist(max_cat_item, use.names = F) * M), ng)),
                            item = c(1:J) %>% rep.int(times = unlist(max_cat_item, use.names = F) * M) %>% rep.int(times = ng),
                            category = cat_item %>% unlist(use.names = F) %>% rep.int(times = rep(M, length(.))) %>% rep.int(times = ng),
                            node = c(1:M) %>% rep.int(times = length(unlist(cat_item))) %>% rep.int(times = ng),
                            r = rgjkm_vec,
                            N = Ngm_vec,
                            P = Pgjkm_vec,
                            Z_cum = Zgjkm_cum_vec,
                            Z_bar = Zgjm_vec,
                            T_bar = Tgjm_vec
                            )



a1 <- numeric(J)
b1 <- numeric(J)
k1 <- matrix(0, nrow = J, ncol = max_cat_all)
for(j in 1:J){
  # 1st Mstep
  cat(j,"\n")
  response_group <- which(design_each_g[,j] == 1)
  X_temp <- numeric(0)
  for(g in response_group){
    X_temp <- X_temp %>% dplyr::bind_rows(X_long %>% dplyr::filter(item == j) %>% dplyr::filter(group == g))
  }
  X_temp %<>% dplyr::slice(-1)
  convergence <- TRUE
  t <- 0
  while(convergence){
    t <- t + 1
    cat(t,"  ")
    gr <- numeric(2)
    gr[1] <- a0[j]^(-1) * sum(X_temp$r * (X_temp$Z_cum - X_temp$Z_bar))
    gr[2] <- a0[j] * sum(X_temp$r * (-X_temp$category + X_temp$T_bar))
    V <- matrix(0,2,2)
    V[1,1] <- a0[j]^(-2) * sum(X_temp$N * X_temp$P * (X_temp$Z_cum - X_temp$Z_bar)^2)
    V[2,2] <- a0[j]^2 * sum(X_temp$N * X_temp$P * (-X_temp$category + X_temp$T_bar)^2)
    V[1,2] <- V[2,1] <- sum(X_temp$N * X_temp$P * (X_temp$Z_cum - X_temp$Z_bar) * (-X_temp$category + X_temp$T_bar))

    t1 <- c(a0[j], b0[j]) + solve(V) %*% gr
    a1[j] <- t1[1]
    b1[j] <- t1[2]
    if(abs(a0[j] - a1[j]) < 0.001 && abs(b0[j] - b1[j]) < 0.001){
      convergence <- FALSE
    } else {
      cat("a ", abs(a0[j] - a1[j]), "b", abs(b0[j] - b1[j]), "\r")
      a0[j] <- a1[j]
      b0[j] <- b1[j]
    }
  }


  # 2nd Mstep
  temp_N <- Ngm[response_group,] %>% t() %>% as.vector()
  for(k in cat_item[[j]]){
    rjkm_sum <- X_temp %>% dplyr::select(group, node, r) %>% dplyr::group_by(group, node) %>% dplyr::summarise(r = sum(r)) %>% .[,3] %>% apply(1, rep.int, times = 4) %>% as.vector()
    gr2 <- a0[j] * sum(X_temp$r - X_temp$P * rjkm_sum)
    Pjkm_sum <- X_temp %>% dplyr::select(group, node, category, P) %>% dplyr::filter(category >= k) %>% dplyr::group_by(group, node) %>% dplyr::summarise(P = sum(P)) %>% .[,3]
    Pjk1m_sum <- X_temp %>% dplyr::select(group, node, category, P) %>% dplyr::filter(category >= k-1) %>% dplyr::group_by(group, node) %>% dplyr::summarise(P = sum(P)) %>% .[,3]
    V2 <- a0[j]^2 * sum(temp_N *Pjkm_sum * (1 - Pjk1m_sum))
    k1[j, k] <- k0[j,k] + gr2/V2
  }
}

