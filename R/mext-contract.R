
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

#' a table for Equipercentile equating function of row score.
#' The return table contains roe score, frequency, cumrative frequency, percent, cumlative percent, discrete distribution and cumlative discrete distribution.
#' @param x test score vector
#' @author Takumi Shibuya.
#' @export
dist.f <- function(x){
  mxc <- max(x)
  m <- length(0:mxc)

  freq <- data.frame(rep(0,m), row.names = as.character(0:mxc))
  colnames(freq) <- "freq"
  for(i in 0:mxc){
    f <- sum(1 * (x == i))
    freq[i+1,] <- f
  }
  score <- 0:mxc
  cum.freq <- cumsum(freq)
  percent <- freq/sum(freq)*100
  cum.pcnt <- cumsum(percent)
  ddf <- freq/length(x)
  cddf <- cumsum(ddf)

  result <- cbind(score,freq, cum.freq, percent, cum.pcnt, ddf, cddf)
  colnames(result) <- c("score","freq", "cum.freq", "percent", "cum.pcnt", "ddf", "cddf")
  #result <- data.frame(score=m, freq=freq, cum.freq=cum.freq, percent=percent, cum.pcnt=cum.pcnt, ddf=ddf, cddf=cddf)
  as.data.frame(result)
}

#' percentile rank function.
#'
#' @param q test raw score vector.
#' @param table This must be distribution table generated by `dist.f` function.
#' @author Takumi Shibuya.
#' @export
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

#' The inverse of percentile rank function which uses the smallest integer score with cumulative persent that is greater then `p`.
#' This function generates `p` percentile of the score in `tabley2`.
#' @param p percentile.
#' @param tabley2 distribution table(must be data.frame). `V2` is cumulative discrete distribution.
#' @author Takumi Shibuya.
#' @export
pfU <- function(p,tabley2){
  if(p/100 > tabley2$V2[1]){
    a <- tabley2[tabley2$V2 >= (p/100),]
    yU <- a[1,1]
    FyU <- a$V2[1]
    b <- tabley2[tabley2$V2 < (p/100),]
    FyU1 <- b$V2[yU]
    pfU <- (p/100 - FyU1)/(FyU-FyU1) + yU - 0.5
  }else if(p/100 == 0){
    pfU <- 0
  }else{
    pfU <- p/100 / tabley2$V2[1] + tabley2$V1[1]-0.5
  }
  pfU
}


# The inverse of percentile rank function which uses the largest integer score with cumulative persent that is less then `p`.

pfL <- function(p,tabley2){
  if(p/100 > tabley2$V2[1]){
    a <- tabley2[tabley2$V2 >= (p/100),]
    yL1 <- a[1,1]
    FyL1 <- a$V2[1]
    b <- tabley2[tabley2$V2 < (p/100),]
    FyL <- b$V2[yL1]
    pfL <-(p/100 - FyL)/(FyL1-FyL) + (yL1 - 1) + 0.5
  }else if(p/100 == 0){
    pfL <- 0
  } else {
    pfL <- p/100 / tabley2$V2[1] + tabley2$V2[1] - 1 + 0.5
  }
  pfL
}


#' Equipercentile equating function of row score.
#' This function generates the Form Y equipercentle equivalent of score x on FormX, eY(x).
#' @param x test score of Form X
#' @param y test score of Form Y
#' @author Takumi Shibuya.
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
  #----------------------------------------------------------------------------------

  # frequency distribution table
  tablex <- dist.f(x)
  tabley <- dist.f(y)

  #n of item on Form X and Y
  nx <- max(tablex$score)
  ny <- max(tabley$score)

  resultx <- matrix(0, nrow(tablex), 1)
  resulty <- matrix(0, nrow(tabley), 1)
  #resultx[1,1] <- tablex$percent[1]
  #resulty[1,1] <- tabley$percent[1]

  for (i in 0:nx){
    resultx[i+1, 1] <- prf(i,tablex)
  }

  for (i in 0:ny){
    resulty[i+1, 1] <- prf(i,tabley)
  }

  #----------------------------------------------------------------------------------
  # equating function
  # test Form Y equipercentile equivalent of score x on Form X
  #----------------------------------------------------------------------------------

  tabley2 <- matrix(0, nrow(tabley), 3)
  tabley2[,1] <- tabley$score
  tabley2[,2] <- tabley$cddf      #relative.p
  tabley2[,3] <- resulty    #percentile rank
  colnames(tabley2) <- c("V1","V2","V3")
  tabley2 <- as.data.frame(tabley2)

  #	paste(tabley2)

  #----------------------------------------------------------------------------------
  #^ "tabley2 " is used for "pf "
  #Dont forget to make "table2" as data frame.
  #
  #define the percentile function that contains yL
  #(This value is smallest value corresponds with above cumulative percent "p"%
  #----------------------------------------------------------------------------------

  eYx <- matrix(0, nrow(tablex), 2)
  for (i in 1:nrow(tablex)){
    eYx[i,1] <- pfU(resultx[i],tabley2)
    eYx[i,2] <- pfL(resultx[i],tabley2)
  }

  result <- cbind(tablex$score, eYx)
  colnames(result) <- c("X", "eYx_U", "eYx_L")
  as.data.frame(result)

}
