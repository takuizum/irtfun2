
#' a function generates conditional probability for fixed ability on 2-parameter logistic model
#' The original Fortran77 program was developed by Inoue,S., December 1990., extended by Shibayama,T., January 1991., translated into R by Shibayama,T. September 2008., functionalized by Itamiya, C., & Shibuya. T., June 2018.
#' @param trait a vector of theta on IRT.
#' @param a a slope parameter
#' @param b a location parameter
#' @param D a factor constant
#' @export

probability <- function(trait,a,b,D=1.702){
  m <- nrow(a)
  m1 <- m+1
  ptheta <- matrix(0,m,1)
  qtheta <- matrix(0,m,1)
  prb <- matrix(0,m1,1)
  #
  ptheta<- 1/(1+exp(-D*a*(trait-b)))
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
#' @param D a factor constant
#' @param output int. if 1 score vector, if 2 cumulative distribution plot.
#' @param name a plot title
#' @param color a plot color.
#' @export


obscore_dist <- function(theta,a,b,D=1.702,name="test",color="cyan", output=1){

  a <- as.matrix(a)
  b <- as.matrix(b)　　　
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
    prbtestscore[i,]<- t(probability(theta[i],a,b,D=D))
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
truescore<-function(trait,a,b,D){

  m <- nrow(a)
  m1 <- m+1
  ptheta<- 1/(1+exp(-D*a*(trait-b)))
  truescore <- sum(ptheta)
}


#' A function calculates IRT true score.
#'
#' @param theta a vector of theta estimator EAP, MAP or MLE...
#' @param a a slope parameter.
#' @param b a location parameter
#' @param D a factor constant
#' @export

tscore_dist <- function(theta,a,b,D=1.702){
  #--------------------------------------------------------------#
  a <- as.matrix(a)
  b <- as.matrix(b)　　　
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
    tscore[i,] <- t(truescore(theta[i],a,b,D=D))
  }

  return(round(tscore))
}



#' Equipercentile equating function of row score.
#' This function generates the Form Y equipercentle equivalent of score x on FormX, eY(x).
#' @param x test score of Form X
#' @param y test score of Form Y
#' @export

epe <- function(x, y){

  #----------------------------------------------------------------------------------
  #
  #		This function generates the Form Y equipercentle equivalent of score x on FormX, eY(x).
  #		eY(x) is called equating function.
  #
  #		x : 			test score of Form X
  #		nx:			number of items on Form X
  #		y : 			test score of Form Y
  #		ny			number of items on Form Y
  #		dist.f:			function for generating the frequency distribution
  #		prfx:			function for generating the percentile rank on score x in test form X
  #		prfy:			function for generating the percentile rank on score y in test form Y
  #		pfL:			the inverse of the percentile rank function.
  #		pfU:			the inverse of the percentile rank function.
  #		eYx:			the function of the equipercentile equivalent of score x on the scale of Form Y.
  #		eXy:			the function of the equipercentile equivalent of score y on the scale of Form X.
  #		tablex:
  #		tabley:
  #		resultx:
  #		resulty:
  #----------------------------------------------------------------------------------

  tablex <- x;tabley <- y

  prfx <- function(q){
    x <- q
    x1 <- q+1
    rx <- round(x1)
    ddf <- tablex[,3]
    cddf <- tablex[,5]
    prf <- (cddf[rx]+(x-rx+0.5)*(cddf[rx]-cddf[rx-1]))*100
  }

  prfy <- function(r){
    x <- r
    x1 <- r+1
    rx <- round(x1)
    ddf <- tabley[,3]
    cddf <- tabley[,5]
    prf <- (cddf[rx]+(x-rx+0.5)*(cddf[rx]-cddf[rx-1]))*100
  }

  #n of item on Form X and Y

  nx <- max(x$score)
  ny <- max(y$score)

  resultx <- matrix(0, nrow(x), 1)
  resulty <- matrix(0, nrow(y), 1)
  resultx[1,1] <- x[1,3]
  resulty[1,1] <- y[1,3]

  z<- 0
  for (i in 1:nx){
    z <- z +1
    resultx[(z+1), 1] <- prfx(z)
  }

  z<- 0
  for (i in 1:ny){
    z <- z +1
    resulty[(z+1), 1] <- prfy(z)
  }

  #----------------------------------------------------------------------------------
  # equating function
  # test Form Y equipercentile equivalent of score x on Form X
  #----------------------------------------------------------------------------------

  tabley2 <- matrix(0, nrow(y), 3)
  tabley2[,1] <- y$score
  tabley2[,2] <- y[,5]      #relative.p
  tabley2[,3] <- resulty    #percentile rank
  tabley2 <- as.data.frame(tabley2)

  #	paste(tabley2)

  #----------------------------------------------------------------------------------
  #^ "tabley2 " is used for "pf "
  #Dont forget to make "table2" as data frame.
  #
  #define the percentile function that contains yL
  #(This value is smallest value corresponds with above cumulative percent "p"%
  #----------------------------------------------------------------------------------

  pfU <- function(p){
    if(p/100 > tabley2$V2[1]){
      a <- tabley2[tabley2$V2 > (p/100),]
      yU <- a[1,1]
      FyU <- a$V2[1]
      b <- tabley2[tabley2$V2 < (p/100),]
      FyU1 <- b$V2[yU]				#yU-1 is not proper value because table b couont "0".
      pfU <- (p/100 - FyU1)/(FyU-FyU1) + yU - 0.5
    }else{
      pfU <- tabley2$V2[1]
    }
  }

  #----------------------------------------------------------------------------------
  #Next, define the percentile function that contains yL
  #(This value is largest value corresponds with below cumulative percent "p" percent)
  #----------------------------------------------------------------------------------

  pfL <- function(p){
    if(p/100 > tabley2$V2[1]){
      a <- tabley2[tabley2$V2 > (p/100),]
      yL1 <- a[1,1]
      FyL1 <- a$V2[1]
      b <- tabley2[tabley2$V2 < (p/100),]
      FyL <- b$V2[yL1]
      pfL <-(p/100 - FyL)/(FyL1-FyL) + (yL1 - 1) + 0.5
    }else{
      pfL <- tabley2$V2[1]
    }

  }

  eYx <- matrix(0, nrow(x), 2)
  for (i in 1:nrow(x)){
    eYx[i,1] <- pfU(resultx[i])
    eYx[i,2] <- pfL(resultx[i])
  }

  result <- cbind(x$score, eYx)
  colnames(result) <- c("X raw score", "pfU", "pfL")
  result
}
