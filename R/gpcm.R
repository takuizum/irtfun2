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

if(data_type == "f")
Science %>% purrr::map_df(as.integer)

X <- x %>% switch("f" = x %>% purrr::map_df(as.integer), # if factor type
                  "i" = x)

X <- Science %>% purrr::map_df(as.integer)#%>% as.matrix()
X %<>% dplyr::bind_rows(X) # %>% dplyr::mutate(id = 1:nrow(.)) # bind row (add more data lines)

N <- nrow(X)
J <- ncol(X)
Item <- colnames(X)

# list up all cat and min and max cat
cat_item <- X %>% purrr::map(unique)
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
b0 <- matrix(0, nrow = J, ncol = length(min_cat_all:max_cat_all))
a0 <- rep(1, J)
cat_count <- X %>% purrr::map(~ as.vector(table(.)))
cat_count <- X %>% purrr::map(table)

for(j in 1:J){
  cat <- dimnames(cat_count[[j]])[[1]] %>% as.integer()
  cat_j <- cat_count[[j]] %>% as.integer()
  prob <- cat_j / (group_N[design_each_g[,j] == 1,] %>% dplyr::select(count) %>% sum())
  b0[j, cat] <- -log(prob / (1 - prob))
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
group_vec <- Xl_tbl %>% dplyr::arrange(group) %>% .$group

# create ujk
ujk <- mapply(create_ujk, Xl_long$response, n_cat_vec) %>% as.vector()
ujk <- tibble(u = as.integer(ujk))


# prior of theta
M <- 31 # of node
min_th <- -4
max_th <- 4
xm <- seq(min_th, max_th, length.out = M)
wm <- dnorm(xm)/sum(dnorm(xm)) %>% rep.int(times = ng) %>% matrix(nrow = M, ncol = ng, byrow = T)

# EM algorithm starts

# E step

# logit and probability
Pjk <- array(0, dim = c(ng, J, max_cat_all, M))
Zjk <- array(0, dim = c(ng, J, max_cat_all, M))
Zjk_cum <- array(0, dim = c(ng, J, max_cat_all, M))
Zj_bar <- array(0, dim = c(ng, J, M))
Tj <- array(0, dim = c(ng, J, M))
Zj <- array(0, dim = c(ng, J, M))
for(g in 1:ng){
  item_list <- which(design_each_g[g, ] == 1) # pattern location of group g in data.frame
  for(j in item_list){
    cat_list <- sort(cat_item[[j]])
    for(m in 1:M){
      for(k in cat_list){ # calculate Z_jk
        Zjk[g, j, k, m] <- 1.702 * a0[j] * (xm[m] - b0[j, k])
      }
      temp <- cumsum(Zjk[g, j, , m])
      Zjk_cum[g, j, , m] <- temp
      p_const <- 1 + sum(exp(temp))
      for(k in cat_list){ # calculate Zjk
        Pjk[g, j, k, m] <- (1 + exp(Zjk_cum[g, j, k, m])) / p_const
      }
      Zj[g, j, m] <- Pjk[g, j, , m] %*% temp
      Tj[g, j, m] <- Pjk[g, j, , m] %*% cat_list
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
        p <- p * Pjk[g,j,Xl[l,j],m]
      }
      Lw[g, l, m] <- p
      rLw[g, l, m] <- p * rl[l]
      const[g, l,] <- const[g, l,] + p
    }
  }
}

# system.time(
#   for(l in 1:L){
#     for(g in 1:ng){
#       if(group_id[l] != g) next
#       # cat("pattern ", l, "\r")
#       for(m in 1:M){
#         p <- wm[m, g]
#         for(j in 1:J){
#           if(design_all[l, j] == 0) next
#           # p <- p * gpcm(xm[m], a = a0[j], b = b0[j,], k = Xl[l, j], D = 1.702) # irp of category k
#           p <- p * Pjk[g,j,Xl[l,j],m]
#         }
#         Lw[g, l, m] <- p
#         rLw[g, l, m] <- p * rl[l]
#         const[g, l,] <- const[g, l,] + p
#       }
#     }
#   }
# )


# provisional sumple size
Nf <- numeric(M)
for(g in 1:ng){
  Nf <- Nf + (colSums(rLw[g,group_id == g,]/const[g, group_id == g,]))
}

# expected frequency of subjects that responsed category k in item j
rgjkm <- array(0, dim = c(ng, J, max_cat_all, M))

# for(l in 1:L){
#   for(g in 1:ng){
#     if(group_id[l] != g) next
#     cat("pattern ", l, "\r")
#     for(m in 1:M){
#       p <- wm[m, g]
#       for(j in 1:J){
#         if(design_all[l, j] == 0) next
#         rgjkm[g, j, Xl[l, j], m] <- rgjkm[g, j, Xl[l, j], m] + rLw[g, l, m] / const[g, l, m]
#       }
#     }
#   }
# }
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

rjkm <- array(0, dim = c(J, max_cat_all, M))
for(g in 1:ng){
  rjkm <- rjkm + rgjkm[g,,,]
}


# M step algorithm

for(j in 1:J){
}
