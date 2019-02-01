
#' a function generates conditional probability for fixed ability on 2-parameter logistic model
#'
#' @param trait a vector of theta in IRT.
#' @param a a slope parameter
#' @param b a location parameter
#' @param c a lower asymptote parameter
#' @param D a factor constant
#' @author The original Fortran77 program was developed by Inoue,S., December 1990., extended by Shibayama,T., January 1991., translated into R by Shibayama,T. September 2008., functionalized by Itamiya, C., & Shibuya. T., June 2018.
#' @references Kolen, M. J., & Brennan, R. L. (2014). Test Equating, Scaling, and Linking. Springer.
#' @export

probability <- function(trait,a,b,c,D=1.702){
  m <- a %>% as.vector() %>% length()
  m1 <- m+1
  ptheta <- matrix(0,m,1)
  qtheta <- matrix(0,m,1)
  prb <- matrix(0,m1,1)
  #
  ptheta<- c+(1-c)/(1+exp(-D*a*(trait-b)))
  qtheta<- 1-ptheta
  #
  prb[1] <- qtheta[1]
  prb[2] <- ptheta[1]
  # recursive formula
  for(j in 2:m){
    l <- j -1
    j1 <- j+1
    l1 <- l+1
    prb[j1] <- prb[l1]*ptheta[j]

    for(i in l:1){
      k <- i -1
      i1 <- i+1
      k1 <- k+1
      prb[i1]<-prb[k1]*ptheta[j]+prb[i1]*qtheta[j]
    }

    prb[1] <- prb[1]*qtheta[j]
  }
  #
  probability <- prb
  #
}



#' A function calculates IRT observed score using recursion formula.
#'
#' @param theta a vector of theta estimator EAP, MAP or MLE...
#' @param a a slope parameter.
#' @param b a location parameter
#' @param c a lower asymptote parameter
#' @param D a factor constant
#' @param output int. if 1 score vector, if 2 cumulative distribution plot.
#' @param name a plot title
#' @param color a plot color.
#' @export


obscore_dist <- function(theta,a,b,c,D=1.702,name="test",color="cyan", output=1){

  a <- as.matrix(a)
  b <- as.matrix(b)
  c <- as.matrix(c)
  theta <- as.matrix(theta)


  #--------------------------------------------------------------#
  # Number of Items and Subjects
  #--------------------------------------------------------------#
  m <- nrow(a)
  m1 <- m+1
  n <- nrow(theta)

  #--------------------------------------------------------------#
  # Distribution of Ability
  #--------------------------------------------------------------#
  i <- sort.list(theta)
  theta <- theta[i]
  i <- matrix(1:n,1,n)


  #--------------------------------------------------------------#
  # Conditional Probabilities of Test Scores
  #--------------------------------------------------------------#

  prbtestscore<-matrix(0,n,m1)
  for(i in 1:n){
    prbtestscore[i,]<- t(probability(theta[i],a,b,c,D=D))
  }

  #--------------------------------------------------------------#
  # Marginal Probabilities of Test Score
  #--------------------------------------------------------------#
  freq <- t(prbtestscore) %*% matrix(1,n,1)
  freq <- cbind(matrix(0:m,m1,1),freq)
  temp <- round(freq[,2])
  score <- rep(freq[,1],temp)

  if(output == 1){
    return(score)
  }else if(output==2){
    mx <- max(score)
    graphics::hist(score,freq=FALSE,ylim=c(0,0.15),breaks=seq(-0.5,(mx+0.5),1),col=color,main=name,cex.main=1.5)
  }
}


# Comditional Probabilities of True Scores
truescore<-function(trait,a,b,c,D){
  ptheta<- c+(1-c)/(1+exp(-D*a*(trait-b)))
  truescore <- sum(ptheta, na.rm = T)
}


#' A function calculates IRT true score.
#'
#' @param theta a vector of theta estimator EAP, MAP or MLE...
#' @param a a slope parameter.
#' @param b a location parameter
#' @param c a lower asymptote parameter
#' @param D a factor constant
#' @export

tscore_dist <- function(theta,a,b,c,D=1.702){
  #--------------------------------------------------------------#
  a <- as.matrix(a)
  b <- as.matrix(b)
  c <- as.matrix(c)
  theta <- as.matrix(theta)

  # Number of Items and Subjects
  m <- length(a) #n of items
  m1 <- m+1
  n <- length(theta) #n of subjects

  # Distribution of Ability
  i <- sort.list(theta)
  theta <- theta[i]
  i <- matrix(1:n,1,n)

  tscore <- matrix(0,n,1)

  for(i in 1:n){
    tscore[i,] <- t(truescore(theta[i],a,b,c,D=D))
  }

  return(as.vector(round(tscore)))
}

dist_f <- function(x, mxc = NULL, mnc = NULL){
  x <- as.vector(x) # vectorization
  if(is.null(mnc)) mnc <- 0
  if(is.null(mxc)) mxc <- max(x)
  score <- mnc:mxc
  m <- length(score)
  if(min(x) < mnc || max(x) < mxc) stop("max or min of score is incorrect!")

  res <- data.frame(score = score, freq = rep(0,m))
  f <- numeric(m)
  for(i in score){
    f[i+1] <- sum(x == i)
  }
  res$freq <- f

  res <- res %>%
    dplyr::mutate(cum_freq = cumsum(res$freq)) %>%
    dplyr::mutate(percent = res$freq/sum(res$freq)*100) %>%
    dplyr::mutate(cum_percent = cumsum(percent)) %>%
    dplyr::mutate(ddf = res$freq/length(x)) %>%
    dplyr::mutate(cddf = cumsum(ddf))
  res
}


prf <- function(q,table){

  xast <- round(q)
  x1 <- q+1 #
  x1ast <- round(x1)
  cddf <- table$cddf

  # パーセンタイルランクの計算
  if(x1ast == 1){    # 0以下の度数は存在しないため0と置いた。
    prf <- (0+(q-xast+0.5)*(cddf[x1ast]-0))*100

  } else if(x1ast >= 2){
    prf <- (cddf[x1ast-1]+(q-xast+0.5)*(cddf[x1ast]-cddf[x1ast-1]))*100

  } else if(x1ast == 0){ #　-0.5以下の場合は0
    prf <- 0
  }
  # 最大値＋0以上の場合は100。
  # ただし，最大値＋0.5以上の値を入れるとprfがNAを返す性質を利用した。
  if (is.na(prf)) prf <- 100
  prf
}

pfU <- function(p,tabley2){
  if(0 <= p && p < 100){
    a <- tabley2[tabley2$cddf > (p/100),]
    yU <- a$score[1]
    FyU <- a$cddf[1]
    b <- tabley2[tabley2$cddf <= (p/100),]
    FyU1 <- b$cddf[yU]
    if(nrow(b) == 0) FyU1 <- 0
    pfU <- (p/100 - FyU1)/(FyU-FyU1) + yU - 0.5
  }else{
    pfU <- tabley2$score %>% max() %>% magrittr::add(0.5)
  }
  pfU
}


# The inverse of percentile rank function which uses the largest integer score with cumulative persent that is less then `p`.

pfL <- function(p,tabley2){
  if(0 < p && p <= 100){
    a <- tabley2[tabley2$cddf >= (p/100),]
    yL1 <- a[1,1]
    FyL1 <- a$cddf[1]
    b <- tabley2[tabley2$cddf < (p/100),]
    FyL <- b$cddf[yL1]
    if(nrow(b) == 0) FyL <- 0
    pfL <-(p/100 - FyL)/(FyL1-FyL) + (yL1 - 1) + 0.5
  }else if(p == 0){
    pfL <- -0.5
  }
  pfL
}


#' Equipercentile equating function of row score.
#' This function generates the Form Y equipercentle equivalent of score x on FormX, eY(x).
#' @param x integer vector. test score of Form X
#' @param y integer vector. test score of Form Y
#' @param type character. if "U", result percentile is calculated by the smmallest integer score with a cum_percent that is  greater then p,
#' "L", the largest score that is less than p or "both" output both of them.
#' @author Takumi Shibuya.
#' @examples
#' set.seed(0204)
#' X <- round(rnorm(1000) * 10 + 50)
#' Y <- round(rnorm(900) * 9 + 40)
#' res <- epe(x = X, y = Y)
#'
#' set.seed(0507)
#' X <- round(rnorm(1000) * 10 + 40)
#' X[X < 0] <- 0
#' Y <- round(rnorm(900) * 10 + 40)
#' res2 <- epe(x = X, y = Y)
#' @export

epe <- function(x, y, type = "both"){

  # frequency distribution table
  tablex <- dist_f(x)
  tabley <- dist_f(y)

  #n of item on Form X and Y
  nx <- max(tablex$score)
  ny <- max(tabley$score)

  resultx <- matrix(0, nrow(tablex), 1)
  resulty <- matrix(0, nrow(tabley), 1)

  for (i in 0:nx) resultx[i+1, 1] <- prf(i,tablex)
  for (i in 0:ny) resulty[i+1, 1] <- prf(i,tabley)

  tabley2 <- data.frame(score = tabley$score, cddf = tabley$cddf, prf = resulty)

  eYx <- matrix(0, nrow(tablex), 2)
  for (i in 1:nrow(tablex)){
    eYx[i,1] <- pfU(resultx[i],tabley2)
    eYx[i,2] <- pfL(resultx[i],tabley2)
  }

  result <- cbind(tablex$score, eYx)
  if(type == "both"){
    result <- data.frame(x = tablex$score, ey_U = eYx[,1], ey_L = eYx[,2])
  } else if(type == "L"){
    result <- data.frame(x = tablex$score, ey = eYx[,2])
  } else if(type == "U"){
    result <- data.frame(x = tablex$score, ey = eYx[,1])
  }

  list(type = type, table = result, freq_x = tablex, freq_y = tabley, pr = tabley2)
}
