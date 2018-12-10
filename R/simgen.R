#' Calculate subject log likelihood.
#'
#' @param theta theta vector.
#' @param x item response matrix or df.
#' @param a slope parameter.
#' @param b location parameter.
#' @param c asymptote parameter.


Flol <- function(theta,x,a,b,c){
  #cat("識別力が0の項目を削除します。¥n")
  x <- x[,a!= 0]
  #cat(sum((a == 0)*1),"個の項目が削除されました。¥n")
  c <- c[a!=0]
  b <- b[a!=0]
  a <- a[a!=0]
  dat <- cbind(theta,x)
  lolF <- function(dat,a,b){
    theta <- dat[1]
    u <- dat[-1]
    p <- c+(1-c)/(1+exp(-1.702*a*(theta-b)))
    LL <- sum(u*log(p)+(1-u)*log(1-p),na.rm = T)
  }
  apply(dat, 1, lolF, a=a, b=b)
}


#' Generate binary data from response probability vector
#'
#' @param prob response probabirity vector
#' @param power a power of probability of NA. for example power = 1/5
#' @export
subfunc <- function(prob,power){
  if(prob < runif(1)) res <- 0
  else res <- 1
  sample(x=c(res,NA),size=1,prob=c(prob^power,1-prob^power))
}

#' Generate binary data from response probability matrix.
#'
#' @param prob response probabirity matrix
#' @param power a power of probability of NA. for example power = 1/5
#' @export


resfunc <- function(prob,power){
  # 反応確率と一様乱数から01データを発生させる関数。
  # 受検者一人分の正答確率を与える。（apply関数などで）
  # powerの値をいじることで，NAの発生率を変更できる。powerの値が小さいほど，発生率も小さい。
  prob <- matrix(prob, ncol = 1)
  res <- apply(prob, 1, subfunc,power=power)
  return(res)
}


#' Generate simulation binary data in IRT.
#'
#' This function contain `resfunc` and `subfunc`.
#' @param theta theta vector
#' @param phi vector of hyperparameter of phi in GIRT model. Default is `NULL`
#' @param a slope parameter.
#' @param b location parameter.
#' @param c asymptote parameter. Default is `NULL`
#' @param item Character. item code.
#' @param power a power of probability of NA. for example power = 1/5
#' @param D a factor constant.
#' @export

sim_gen <- function(theta, phi=NULL, a, b, c=NULL, item = 'A', power=0, D=1.702){
  if(is.null(phi)){
    # IRT
    if(is.null(c)) c <- rep(0,length(a))
    prob <- t(apply(matrix(theta, ncol = 1), 1, ptheta, a=a, b=b, c=c, D=D))
  }else{
    # GIRT
    prob <- t(mapply(gptheta, matrix(theta, ncol = 1), matrix(phi, ncol=1), MoreArgs = list(a=a, b=b, D=D)))
  }

  # generate item respinse pattern
  test <- t(apply(prob, 1, resfunc, power=power))
  # change item name
  itemn <- formatC(c(1:length(a)), width = 3, flag = 0)
  itemn <- apply(matrix(item, ncol = 1), 1, paste0, itemn)

  test <- cbind(c(1:length(theta)),test)
  colnames(test) <- c("ID",itemn)
  test <- as.data.frame(test)
  return(test)
}

