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

if(data_type == "f")
Science %>% purrr::map_df(as.integer)

X <- x %>% switch("f" = x %>% purrr::map_df(as.integer), # if factor type
              "i" = x)

X <- Science %>% purrr::map_df(as.integer)#%>% as.matrix()

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

# add group id column
X %<>% dplyr::mutate(group = as.integer(1))

# add count column
X %<>% dplyr::mutate(count = 1)

# unique response patterns and them count
Xl <- X %>% group_by_if(is.integer) %>% dplyr::summarise(count = sum(count))

# intial value
b_init <- matrix(0, nrow = J, ncol = length(min_cat_all:max_cat_all))
a_init <- rep(1, J)
cat_count <- X %>% purrr::map(~ as.vector(table(.)))
cat_count <- X %>% purrr::map(table)

for(j in 1:J){
  cat <- dimnames(cat_count[[j]])[[1]] %>% as.integer()
  cat_j <- cat_count[[j]] %>% as.integer()
  prob <- cat_j / N
  b_init[j, cat] <- -log(prob) / (1 - prob)
}

rl <- Xl$count
group <- Xl$group
Xl %<>% dplyr::select(-count, -group) %>% as.matrix()

