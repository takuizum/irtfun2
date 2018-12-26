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

