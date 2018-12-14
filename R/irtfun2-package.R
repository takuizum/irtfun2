#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom stats dnorm
#' @importFrom stats na.omit
#' @importFrom stats optimise
#' @importFrom stats optim
#' @importFrom stats runif
#' @importFrom stats cor
#' @importFrom stats qnorm
#' @importFrom stats dlnorm
#' @importFrom tidyr gather
#' @useDynLib irtfun2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

.onUnload = function(libpath) {
  library.dynam.unload("irtfun2", libpath)
}

#' The ICC of IRT 1~3PLM
#'
#' @param theta the person ability parameter
#' @param a the slope parameter
#' @param b the location parameter
#' @param c the guessing parameter
#' @param D a scale constant
#' @export
#'
# P(theta) in two-parameter logisticmodel
ptheta <- function(theta,a,b,c,D=1.702){
  c+(1-c)/(1+exp(-D*a*(theta-b)))
}

# 対数尤度
#' The log likelihood function of IRT 1~3PLM
#'
#' @param u the item response pattern
#' @param theta the person ability parameter
#' @param a the slope parameter
#' @param b the location parameter
#' @param c the guessing parameter
#' @param D a scale constant
#' @export
#'
LL <- function(u,theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  sum(u*log(p)+(1-u)*log(1-p),na.rm = T)
}

# 対数尤度関数
#' The log likelihood function of Bayesian IRT 1~3PLM
#' a prior distribution is Gaussian(normal) distribution
#' @param u the item response pattern
#' @param theta the person ability parameter
#' @param a the slope parameter
#' @param b the location parameter
#' @param c the guessing parameter
#' @param mu a hyperparameter of prior distribution.
#' @param sigma same as above.
#' @param D a scale constant
#' @export
#'
LL_b <- function(u,theta,a,b,c,mu,sigma,D){
  p <- ptheta(theta,a,b,c,D)
  sum(log(p)*u+log(1-p)*(1-u)-0.5*((theta-mu)/sigma)^2,na.rm = T)
}

# 一階偏微分
fpdLPD <- function(xi,theta,a,b,c,mu,sigma,D){
  p <- ptheta(theta,a,b,c,D)
  D*sum(a*(xi - p)*(p-c)/(p*(1-c)), na.rm = T)-1/sigma^2*(theta-mu)
}

# 二階偏微分
spdLPD <- function(xi,theta,a,b,c,sigma,D){
  p <- ptheta(theta,a,b,c,D)
  D^2*sum(a^2*(p-c)*(xi*c-p^2)*(1-p)/(p^2*(1-c)^2),na.rm = T) - 1/sigma
}

# 一階偏微分
fpd <- function(xi,theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  D*sum(a*(xi - p)*(p-c)/(p*(1-c)), na.rm = T)
}

# 二階偏微分
spd <- function(xi,theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  D^2*sum(a^2*(p-c)*(xi*c-p^2)*(1-p)/(p^2*(1-c)^2),na.rm = T)
}


#' Test Information Function for IRT 1~3PLM
#'
#' @param theta the person ability parameter
#' @param a the slope parameter
#' @param b the location parameter
#' @param c the guessing parameter
#' @param D a scale constant
#' @export
#'
tif <- function(theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  D^2*sum(a^2*(1-p)*(p-c)^2/((1-c)^2*p), na.rm = T) # I
}

FI <- function(theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  I <- D^2*sum(a^2*(1-p)*(p-c)^2/((1-c)^2*p), na.rm = T) # I
  1/sqrt(I)
}


#
#' MAP posterior information for IRT 1~3PLM
#'
#' @param dat the matrix or data.frame which has theta and item response pattern vector
#' @param a the slope parameter
#' @param b the location parameter
#' @param c the guessing parameter
#' @param sigma a hyper parameter
#' @param D a scale constant
#' @export
#'
pitheta <- function(dat,a,b,c,sigma,D){
  # 2PLM irf
  theta <- dat[1]
  xi <- dat[-1]
  p <- ptheta(theta,a,b,c,D)
  I <- -D^2*sum(a^2*(p-c)*(xi*c-p^2)*(1-p)/(p^2*(1-c)^2),na.rm = T) + 1/sigma
  1/sqrt(I)
}


# likelihood function with weight
WL <- function(xi, theta, a, b,c,D){
  p <- ptheta(theta,a,b,c,D)
  pi <- tif(theta,a,b,c,D)
  W <- sqrt(pi)
  #L1 <- sum(xi*log(p)+(xi-1)*log(1-p),na.rm = T)
  #L2 <- sum(xi*log(p),na.rm=T)
  LL <- LL(xi,theta,a,b,c,D)
  LL + log(W)
}

# E_3(theta) : expectation
E3 <- function(theta, a, b,c,D){
  p <- ptheta(theta,a,b,c,D)
  sum(p * log(p) + (1-p) * log(1-p),na.rm = T)
}

# sigma3
S3 <- function(theta, a, b,c,D){
  p <- ptheta(theta,a,b,c,D)
  sqrt(sum(p * (1 - p) * (log(p/(1-p)))^2,na.rm = T))
}

# person fit index: z_3 statistics
pfit <- function(dat,a,b,c,D){
  #decompose dat
  theta <- dat[1]
  xi <- dat[-1]
  (LL(xi, theta, a, b,c,D) - E3(theta, a, b,c,D)) / S3(theta, a, b,c,D)
}

# log of likelifood
lolF <- function(dat,a,b,c,D){
  theta <- dat[1]
  u <- dat[-1]
  p <- c+(1-c)/(1+exp(-D*a*(theta-b)))
  LL <- sum(u*log(p)+(1-u)*log(1-p),na.rm = T)
}

# Max. Density of Posterior Distributions
Fmaxpdc <- function(xi,theta,a,b,c,D){
  theta <- as.matrix(theta)
  xi <- as.numeric(xi)
  # 項目反応データをapplyで与えて，対数尤度を計算する。
  LLm <- apply(theta,1,LL,u=xi,a=a,b=b,c=c,D=D)
  exp(LLm)
}

#the function for "fg"
Ffg <- function(xi,theta,a,b,c,mu,sigma,D){
  p <- ptheta(theta,a,b,c,D)
  exp(sum(xi*log(p) + (1-xi)*log(1-p), na.rm=T))*dnorm(theta,mean=mu,sd=sigma)
}

FthetaEAP <- function(xi,theta,w,a,b,c,D){
  theta <- as.matrix(theta,byrow=T)
  xi <- as.numeric(xi)

  # 対数尤度の計算
  LLm <- apply(theta, 1, LL, u=xi, a=a, b=b, c=c, D=D)
  Lm <- exp(LLm)

  # 事後分布の分母にあたる，定数項
  const <- sum(Lm*w,na.rm = T)

  # 事後分布の重み
  Gm <- Lm*w/const

  # EAP
  eap <- sum(theta*Gm,na.rm = T)

  # 事後標準誤差
  SE <- sqrt(sum((theta-eap)^2*Gm))

  # 重みに対応する分点の値をかけて，和を取る＝EAP推定値
  return(c(eap, SE ,const))
}

FthetaMAP <- function(xi,a,b,c,mu,sigma,D,groupitem, maxtheta=6,mintheta=-6,method="NR"){

  gid <- xi[1]
  xi <- xi[-1]

  mm <- groupitem[gid]

  #初期値を設定。ログオッズ比を用いた。
  if(sum(xi, na.rm = TRUE) == 0){
    t0 <- log(0.5)
  }else if(sum((xi==1)*1,na.rm=T) == length(na.omit(xi))){
    t0 <- log(mm-0.5)
  }else{
    t0 <- log(sum(xi, na.rm = TRUE)/(mm-sum(xi, na.rm = TRUE)))
  }

  d <- 0
  conv1 <- 0
  if(method=="NR"){
    while(conv1==0){
      t1 <- t0 - fpdLPD(xi,t0,a,b,c,mu,sigma,D)/spdLPD(xi,t0,a,b,c,sigma,D)
      if(abs(t1-t0)<0.001 || abs(t0)>10 ||is.nan(t1)) {
        conv1 <- 1
      }else{
        t0 <- t1
        d <- d +1
        if(d > 100) { # フィッシャースコアリングの繰り返しが100回を超えたら，二分法に切り替える。
          conv2 <- 0
          p <- maxtheta ; q <- mintheta
          while(conv2 == 0){
            pf <- fpdLPD(xi,p,a,b,c,mu,sigma,D)
            qf <- fpdLPD(xi,q,a,b,c,mu,sigma,D)
            t1 <- (p+q)/2
            mf <- fpdLPD(xi,t1,a,b,c,mu,sigma,D)
            if (abs(pf-qf)<0.001 || c == 1000){
              cat("change to the bisection method. \r")
              conv2 <- 1
              conv1 <- 1
            }else{
              if(pf*mf>0) {
                p <- t1
                d <- d+1
              }else if(qf*mf>0){
                q <- t1
                d <- d+1
              } else {
                conv2 <- 1
                conv1 <- 1
              }
            }
          }
        }
      }
    }
  } else if(method=="Brent") {
    t1 <- optimise(LL_b,c(mintheta,maxtheta), maximum = T,u=xi,a=a,b=b,c=c,mu=mu,sigma=sigma,D=D)$maximum
  } else if(method=="SANN"){
    t1 <- optim(par=t0,fn=LL_b,method="SANN",control=list(fnscale=-1),u=xi,a=a,b=b,c=c,mu=mu,sigma=sigma,D=D)$par
  } else {
    t1 <-optim(par=t0,fn=LL_b,gr= function(u, theta,a,b,c,mu,sigma,D){
      p <- ptheta(theta,a,b,c,D)
      D*sum(a*(u - p)*(p-c)/(p*(1-c)), na.rm = T)-1/sigma^2*(theta-mu)
    },method=method,control=list(fnscale=-1),u=xi,a=a,b=b,c=c,mu=mu,sigma=sigma,D=D)$par
  }
  t1
}

FthetaMLE <- function(xi,a,b,c,D,groupitem, maxtheta=6,mintheta=-6,method="NR"){

  gid <- xi[1]
  xi <- xi[-1]
  d <- 0
  conv1 <- 0
  mm <- groupitem[gid]
  temp <- sum(xi==1) == sum(xi,na.rm=T)
  #初期値を設定。ログオッズ比を用いた。
  if(sum(xi, na.rm = T) == 0){
    t1 <- NA
    conv1 <- 1
  }else if(sum((xi==1)*1,na.rm=T) == length(na.omit(xi))){
    t1 <- NA
    conv1 <- 1
  }else{
    t0 <- log(sum(xi, na.rm=T)/(mm-sum(xi, na.rm=T)))
  }


  if(method=="NR"){
    while(conv1==0){
      t1 <- t0 - fpd(xi,t0,a,b,c,D)/spd(xi,t0,a,b,c,D)
      if(abs(t1-t0)<0.001 && !is.nan(t1)) {
        conv1 <- 1
      } else {
        t0 <- t1
        d <- d +1
        if(d > 100 || abs(t0)>10 ||is.nan(t1)){ # ニュートン法の繰り返しが100回を超えたら，二分法に切り替える。
          cat("change to the bisection method.\r")
          conv2 <- 0
          p <- maxtheta ; q <- mintheta
          while(conv2 == 0){
            pf <- fpd(xi,p,a,b,c,D)
            qf <- fpd(xi,q,a,b,c,D)
            t1 <- (p+q)/2
            mf <- fpd(xi,t1,a,b,c,D)
            if (abs(pf-qf)<0.001 || c == 1000){
              conv2 <- 1
              conv1 <- 1
            }else if(pf*mf>0) {
              p <- t1
              d <- d+1
            }else if(qf*mf>0){
              q <- t1
              d <- d+1
            } else {
              conv2 <- 1
              conv1 <- 1
            }
          }
        }
      }
    }
  } else if(method=="Brent" && conv1 != 1) {
    t1 <- optimise(LL,c(mintheta,maxtheta), maximum = T,u=xi,a=a,b=b,c=c,D=D)$maximum
  } else if(method=="SANN" && conv1 != 1){
    t1 <- optim(par=t0,fn=LL,method="SANN",u=xi,control=list(fnscale=-1),a=a,b=b,c=c,D=D)$par
  } else if(conv1 != 1){
    t1 <-optim(par=t0,fn=LL,gr=function(u, theta,a,b,c,D){
      p <- ptheta(theta,a,b,c,D)
      D*sum(a*(u - p)*(p-c)/(p*(1-c)), na.rm = T)
    },method=method,u=xi,control=list(fnscale=-1),a=a,b=b,c=c,D=D)$par
  }
  t1
}

FthetaWLE <- function(xi,a,b,c,D,groupitem, maxtheta=6,mintheta=-6){
  gid <- xi[1]
  xi <- xi[-1]
  d <- 0
  conv1 <- 0
  mm <- groupitem[gid]
  #初期値を設定。ログオッズ比を用いた。
  if(sum(xi, na.rm = TRUE) == 0){
    t0 <- log(0.5)
  }else if(sum((xi==1)*1,na.rm=T) == length(na.omit(xi))){
    t0 <- log(mm-0.5)
  }else{
    t0 <- log(sum(xi, na.rm = TRUE)/(mm-sum(xi, na.rm = TRUE)))
  }
  opt <- optimise(WL,interval = c(mintheta,maxtheta),maximum = T,xi=xi,a=a,b=b,c=c,D=D)
  t1 <- opt$maximum
  t1
}


#'A estimation theta function.
#'
#'@param xall a item response data
#'@param param a item parameter file.If class is df, parameter column starts from second column and the order is a, b, c.
#'@param est estimation method option. "EAP","MAP","MLE","PVs."
#'@param nofrands the number of PVs.
#'@param method iteration method option in MLE and MAP. "NR" is Newton Rapthon, "SANN" is Simulated Aannealing, "Brent" is optimization using optimise function
#'@param file a name of output csv file.
#'@param IDc colomn of ID, default is 1.Even If df has no ID columns this function works.
#'@param gc column of group or population numbers. If df has no grouo ID, set 0.
#'@param fc column of first item response, default is 3.
#'@param gh logical. If TRUE, is default, gaussHermite quadrature runs for EAP estimation.
#'@param N the number of nodes. default is 31.
#'@param D a factor constant. deafault is 1.702.
#'@param output logical. If TRUE write csv file. defauli is FALSE.
#'@param maxtheta the maximul value of theta in integration.
#'@param mintheta the minimum value of theta in integration.
#'@param mu a hyperparameter of prior distribution.
#'@param sigma same as above.
#'@param sampling_engine an option of sampling engine. if select "rejection_R", the rejection samplimg method is conducted and this engine loop witten in R lang this is so somewhat slow, but "rejection_Cpp" is for loop in C++ lang so very fast. If select "slice_R" , the slice sampling method will be conducted.
#'@param warm_up the number of iteration times for warm up in slice sampling
#'@return a list has ID, rawscore, theta, se, person fit index Z3 and log likelihood.
#'@author Takumi, Shibuya., Daisuke, Ejiri., Tadashi, Shibayama. in Tohoku University.
#'
#'@importFrom stats sd
#'@importFrom stats var
#'@importFrom utils write.csv
#'@export

estheta <- function(xall, param, est="EAP", nofrands=10, method="NR", file="default", output=FALSE, IDc=1, gc=2, fc=3,
                    gh = TRUE, N = 31, D=1.702, maxtheta = 6, mintheta = -6, mu=0, sigma=1, sampling_engine="rejection_Cpp", warm_up = 100){

  #chack the arguments
  arg_est <- c("EAP","MAP","MLE","WLE","PVs")
  arg_method <- c("NR","Brent","SANN")
  arg_engine <- c("rejection_R","rejection_Cpp","slice_R")

  if(!(est %in% arg_est)) stop("'est' argument is incorrect!")
  if(!(method %in% arg_method)) stop("'method' argument is incorrect!")
  if(!(sampling_engine %in% arg_engine)) stop("'sampling_engune' argument is incorrect!")

  #check the data
  ID <- xall[,IDc]
  if(gc == 0){
    group <- rep(1, nrow(xall))
    G <- 1
    x.all <- as.matrix(xall[,fc:ncol(xall)])
  }else{
    group <- xall[,gc]
    G <- max(as.numeric(group))
    if(is.vector(xall)){
      x.all <- xall[fc:length(xall)]
      x.all <- matrix(x.all,nrow=1)
    }
    x.all <- xall[,fc:ncol(xall)]
  }
  #remove a=0 item parameter and response column

  if(is.data.frame(param)){
    param <- param[,-1]
  }

  a <- param[,1]
  x.all <- x.all[, a != 0]

  # set item parameter data
  param <- param[param[,1] != 0, ]
  a <- param[,1]
  b <- param[,2]
  c <- param[,3]

  # Number of Subjects"
  n <- nrow(x.all)
  ng <- numeric(G)
  # number of items
  m <-length(a)
  xscore <- rowSums(x.all,na.rm=T)
  cat("n of subject is ",n,".\n")
  cat("n of item is ",m,".\n")
  cat(est,"is starting!\n")

  groupitem <- numeric(G)
  for(g in 1:G){
    key <- group==g
    temp <- x.all[key,]
    temp <- temp %>% colSums(na.rm = T)
    groupitem[g] <- temp[temp!=0] %>% length()
  }
  #cat(groupitem)

  if((est == "PVs") || (est == "EAP")){
    if(gh == TRUE){
      npoint <- N
      gh <- pracma::gaussHermite(npoint)
      xnodes <- gh$x*sqrt(2)
      weight <- gh$w/sqrt(pi)  #caribration
    }else{
      npoint <- N
      xnodes <- seq(mintheta, maxtheta, length.out = N)
      weight <- dnorm(xnodes, mu, sigma)
      weight <- weight/sum(weight)
    }
  }

  #----------------------------------------------------------
  #Eatimate EAP
  #----------------------------------------------------------

  if((est == "PVs") || (est == "EAP")){
    # 全受検者のデータをapply関数で与え，EAP推定を実行
    cat("Calculating the marginal probability and EAP.\r")
    eap_apply <- apply(x.all,1,FthetaEAP,a=a,b=b,c=c,theta=xnodes,w=weight,D=D)
    #message("周辺確率およびEAPの計算が終了しました。")
  }


  #----------------------------------------------------------
  #Eatimate MAP
  #----------------------------------------------------------

  if((est == "PVs") || (est == "MAP")){
    # 全受検者のデータをapply関数で与え，MAP推定を実行
    cat("Calculating MAP.\r")
    map_apply <- apply(cbind(group, x.all),1,FthetaMAP,a=a,b=b,c=c,mu=mu,sigma=sigma,D=D,
                       groupitem=groupitem, maxtheta=maxtheta,mintheta=mintheta,method=method)
    #message("MAP推定値の計算が終了しました。")
  }

  #----------------------------------------------------------
  #Eatimate MLE
  #----------------------------------------------------------

  if(est == "MLE"){
    # 全受検者のデータをapply関数で与え，MAP推定を実行
    cat("Calculating MLE.\r")
    mle_apply <- apply(cbind(group, x.all),1,FthetaMLE,a=a,b=b,c=c,D=D,
                       groupitem=groupitem, maxtheta=maxtheta,mintheta=mintheta,method=method)
    #message("MLE推定値の計算が終了しました。")
  }


  #----------------------------------------------------------
  #Eatimate WLE
  #----------------------------------------------------------
  if(est == "WLE"){
    # 全受検者のデータをapply関数で与え，MAP推定を実行
    cat("Calculating WLE.\r")
    wle_apply <- apply(cbind(group, x.all),1,FthetaWLE,a=a,b=b,c=c,D=D,
                       groupitem=groupitem, maxtheta=maxtheta,mintheta=mintheta)
    #message("WLEの計算が終了しました。")
  }

  if(est == "EAP"){

    # EAP
    eap_m <- eap_apply[1,]
    # estimate standard error for EAP estimator
    SE    <- eap_apply[2,] %>% round(digits = 5)
    z3 <- apply(cbind(eap_apply[1,], x.all), 1, pfit, a=a,b=b,c=c,D=D) %>% round(digits = 5)
    lol <-  apply(cbind(eap_apply[1,], x.all), 1, lolF, a=a,b=b,c=c,D=D) %>% round(digits = 5)

    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,EAP=eap_m,const=eap_apply[3,],SE=SE,z3=z3,lol=lol)
    list("res" = result, "EAPmean&sd" = c(mean(eap_m),sd(eap_m)) %>% round(digits = 5))

  }else if(est == "MAP"){
    #estimate standard error for MAP estimator
    SE <- apply(cbind(map_apply, x.all), 1, pitheta, a=a,b=b,c=c,sigma=sigma,D=D) %>% round(digits = 5)
    z3 <- apply(cbind(map_apply, x.all), 1, pfit, a=a,b=b,c=c,D=D) %>% round(digits = 5)
    lol <- apply(cbind(map_apply, x.all), 1, lolF, a=a,b=b,c=c,D=D) %>% round(digits = 5)

    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,MAP=map_apply %>% round(digits = 5), SE=SE, z3=z3, lol=lol)
    list("res" = result, "MAPmean&sd" = c(mean(map_apply), sd(map_apply)) %>% round(digits = 5))

  }else if(est == "WLE"){
    SE <- apply(matrix(wle_apply, ncol = 1), 1, FI, a=a, b=b,c=c,D=D) %>% round(digits = 5)
    z3 <- apply(cbind(wle_apply, x.all), 1, pfit, a=a,b=b,c=c,D=D) %>% round(digits = 5)
    lol <- apply(cbind(wle_apply, x.all), 1, lolF, a=a,b=b,c=c,D=D) %>% round(digits = 5)
    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,WLE=wle_apply %>% round(digits = 5),SE=SE,z3=z3,lol=lol)
    list("res" = result, "WLEmean&sd" = c(mean(wle_apply), sd(wle_apply)))

  }else if(est == "MLE"){
    SE <- apply(matrix(mle_apply, ncol = 1), 1, FI, a=a, b=b,c=c,D=D) %>% round(digits = 5)
    # person fit index
    z3 <- apply(cbind(mle_apply, x.all), 1, pfit, a=a,b=b,c=c,D=D) %>% round(digits = 5)
    lol <- apply(cbind(mle_apply, x.all), 1, lolF, a=a,b=b,c=c,D=D) %>% round(digits = 5)

    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,MLE=mle_apply %>% round(digits = 5), SE=SE, z3=z3,lol=lol)
    list("res" = result, "MLEmean&sd" = c(mean(mle_apply, na.rm = T), sd(mle_apply, na.rm = T)) %>% round(digits = 5))

  }else if( est == "PVs"){
    #--------------------------------------------------
    # Rejection Samling starts here.
    #--------------------------------------------------

    # 使用するベクトルと行列を予め作成しておく。
    pv <- matrix(0,n,nofrands)
    if(sampling_engine=="rejection_Cpp"){
      cat("Sampling Plausible Values based on von Neumann Rejection sampling.\n")
      cat("C++ version.")
      d1 <- eap_apply[1,]
      d2 <- eap_apply[3,]
      pv <- theta_pv(x=x.all, nofrands=nofrands,eap_apply=d1,const_apply=d2,map_apply=map_apply,
                     n=n,maxtheta=maxtheta,mintheta=mintheta,a=a,b=b,c=c,D=D,mu=mu,sigma=sigma)
    } else if( sampling_engine=="rejection_R"){
      cat("Sampling Plausible Values based on von Neumann Rejection sampling.\n")
      cat("R version \n")
      x.all <- as.matrix(x.all)
      for(k in 1:n){
        xi <- x.all[k,]

        #すでに計算してあるEAP推定値と，事後分布の分子を取り出して行列として展開。
        eap   <- eap_apply[1,k]
        const <- eap_apply[3,k]

        #乱数発生時の，P(θ)軸の最大値を設定。
        yheight <- Ffg(xi=xi,theta=map_apply[k],a=a,b=b,c=c,mu=mu,sigma=sigma,D=D)
        yheight <- yheight/const*1.001
        #cat("yheight is ",yheight,"\r")

        # 乱数発生時の，θ軸の最大値を設定。
        zmax <- maxtheta + eap
        zmin <- mintheta + eap

        nofpv <- 0
        while( nofpv < nofrands ){

          y <- runif(n = 1, min = 0, max = yheight)
          z <- runif(n = 1, min = zmin, max = zmax)
          fg <- Ffg(xi=xi,theta=z,a=a,b=b,c=c,mu=mu,sigma=sigma,D=D)
          fgvalue <- fg/const
          #cat(" y is",y," fgvalue is", fgvalue,"\r")
          if( y <= fgvalue){
            nofpv <- nofpv + 1
            pv[k,nofpv] <- z
          }
        }
        cat(k ," / ", n," \r")
      }
    }else if(sampling_engine=="slice_R"){
      cat("Sampling Plausible Values based via slice sampling.\n")
      cat("slice sampling in R\n")
      x.all <- as.matrix(x.all)
      sum_w <- 0 # a parrameter rerated range I(L,R)
      w <- sigma # estimate of the typical size of a slice
      for(k in 1:n){
        xi <- x.all[k,]
        times <- 0
        nofpv <- 0
        # theta軸上の初期値を決定する
        z0 <- runif(1,min=mintheta,max=maxtheta)
        while(nofpv<nofrands){
          # Draw a real value , y, uniformly from (0,fg(z0))
          fg <- Ffg(xi=xi,theta=z0,a=a,b=b,c=c,mu=mu,sigma=sigma,D=D)
          y <- runif(1,0,fg)
          # Find an interval, I=(L,R)
          # stepping out procedure
          # レンジを決めるための初期値
          m_lim <- 10#integer limitating the size of a slice to m_lim*w
          U <- runif(1)
          L <- z0 - w*U
          R <- L + w
          V <- runif(1)
          J <- floor(m_lim*V)
          K <- (m_lim-1)-J
          fL <- Ffg(xi=xi,theta=L,a=a,b=b,c=c,mu=mu,sigma=sigma,D=D)
          while(J>0 && y<fL){
            L <- L-w
            J <- J-1
            fL <- Ffg(xi=xi,theta=L,a=a,b=b,c=c,mu=mu,sigma=sigma,D=D)
          }
          fR <- Ffg(xi=xi,theta=R,a=a,b=b,c=c,mu=mu,sigma=sigma,D=D)
          while(K>0 && y<fR){
            R <- R+w
            K <- K-1
            fR <- Ffg(xi=xi,theta=R,a=a,b=b,c=c,mu=mu,sigma=sigma,D=D)
          }
          # shrinkage procedure
          Lb <- L
          Rb <- R
          repeat {
            U <- runif(1)
            z1 <- Lb+U*(Rb-Lb)
            fz <- Ffg(xi=xi,theta=z1,a=a,b=b,c=c,mu=mu,sigma=sigma,D=D)
            if(y<fz){
              times <- times + 1
              if(times > warm_up){
                nofpv <- nofpv + 1
                pv[k,nofpv] <- z1 # accept
              }
              # update parameter w and z0
              sum_w <- sum_w+abs(z0-z1)
              w <- sum_w/k
              z0 <- z1
              break
            }
            if(z1<z0) Lb <- z1
            else Rb <- z1
          }
        }
        cat(k ," / ", n," \r")
      }
    }


    #推算値の組ごとにもとめた統計量の平均を計算。なお，母集団が複数存在する場合には母集団ごとに求める。
    pvG <- data.frame(group,pv)
    group_mean <- matrix(0,G,nofrands)
    group_var <- matrix(0,G,nofrands)
    group_sd <- matrix(0,G,nofrands)

    for (i in 1:G){
      # count the number of subjects for each groups
      ng[i] <- nrow(pvG)

      # estiamte group statistics
      group_pv <- pvG %>% dplyr::filter(group==i) %>% dplyr::select(-group)
      group_mean[i,] <- apply(group_pv,2,mean)
      group_var[i,] <- apply(group_pv,2,var) # unbias variance
      group_sd[i,] <- apply(group_pv,2,sd)
    }

    PS <- list(group=c(1:G),mean=group_mean,sd=group_sd,var=group_var)


    #message("集団統計量の標準誤差の推定が終了しました。")
    #推算値の平均(Right & Wrong)
    pvmeans <- pvmeans_w <- apply(pv,1,mean)
    pvmeans_r <- apply(pv,2,mean)

    #result

    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,EAP=eap_apply[1,],MAP=map_apply,PVmeans_W=pvmeans_w,PV=pv,AREAS=eap_apply[3,])
    plausible_values <-data.frame(ID=ID,Group=group,PV=pv)
    if(output==T){
      #message("結果のCSVファイルを書き出しています。")
      cat("\n output result files.")
      write.csv(result,paste0(file,"_result.csv"),quote=F,row.names = F)
      write.csv(plausible_values,paste0(file,"_PVs.csv"),quote=F,row.names = F)
      write.csv(PS,paste0(file,"_PVS population estimator.csv"),quote=F,row.names = F)
    }
    pv <- data.frame(ID=ID,group=group,SCORE=xscore,EAP=eap_apply[1,],MAP=map_apply,PV=pv)
    list("PVs_df" = pv, "PVs_only"=pvG, "EAPmean&sd" = c(mean(eap_apply[1,]),sd(eap_apply[1,])),
         "MAPmean&sd" = c(mean(map_apply), sd(map_apply)), "PS" = PS)
  }
}

e_v_mean <- function(X){
  # 標本平均の分散＝平均の誤差分散
  n <- length(X)
  s <- var(X)*(n-1)/n
  s/n
}

e_v_var <- function(X){
  # 標本分散の分散＝分散の誤差分散
  # 不偏分散を使う
  # reference
  # http://www.crl.nitech.ac.jp/~ida/research/memo/SampleVariance/sample_variance.pdf
  n <- length(X)
  M <- mean(X)
  b4 <- (X-M)^4
  COEF <- 1/((n-1)*(n^2-3*n+3))
  term1 <- n*sum(b4)
  term2 <- (n^2-3)*var(X)^2
  COEF*(term1-term2)
}

#'Calculate the imputation variance
#'
#'@param PVs_only a df given by the result of `estheta(est="PVs")`.
#'@export
stat_pv <- function(PVs_only){
  pvG <- PVs_only
  G <- max(pvG$group)
  MEAN <- numeric(G)
  VAR <- numeric(G)
  res_M <- numeric(G)
  res_V <- numeric(G)

  for(i in 1:max(pvG$group)){
    # extract the PVs data frame for each group
    group_pv <- pvG[pvG$group==i,-1]#
    # the average of the mean
    M <- apply(pvG,2,mean) %>% as.vector()#as.vector(apply(pvG,2,mean))
    V_M <- apply(pvG,2,e_v_mean) %>% as.vector()
    V <- apply(pvG,2,var) %>% as.vector()#as.vector(apply(pvG,2,var))
    V_V <- apply(pvG,2,e_v_var) %>% as.vector()
    # imputation variance of mean
    # variance of expectation(the first term of equation 2.2.2)
    K <- length(M)
    V_imp_mean <- (1+1/K)*var(M)+mean(V_M)
    # imputation variance of variance
    # variance of expectation(the first term of equation 2.2.2)
    V_imp_var <- (1+1/K)*var(V)+mean(V_V)

    res_M[i] <- V_imp_mean
    res_V[i] <- V_imp_var
    MEAN[i] <- mean(M)
    VAR[i] <- mean(V)
  }
  data.frame(group=c(1:G),mean=MEAN,var=VAR,V_imp_mean=V_imp_mean,V_imp_var=V_imp_var)
}

cf <- function(x,af,bf,at,bt,D,q,N,w){
  A <- x[1]
  K <- x[2]
  Q1 <- 0
  Q2 <- 0
  for(m in 1:N){
    HC1 <- (ptheta(q[m],at,bt,c=0,D)-ptheta(q[m],af,bf,c=0,D))^2
    HC2 <- (ptheta(q[m],af,bf,c=0,D)-ptheta(q[m],at*A,(bt-K)/A,c=0,D))^2

    Q1 <- Q1 + sum(HC1*w[m])
    Q2 <- Q2 + sum(HC2*w[m])
  }
  Q <- Q1 +Q2
  return(Q)
}

cf2 <- function(x,af,bf,at,bt,D,q,N,w){
  A <- x[1]
  K <- x[2]
  Q1 <- 0
  Q2 <- 0
  for(m in 1:N){
    SLC1 <- (sum(ptheta(q[m],at,bt,c=0,D))-sum(ptheta(q[m],af,bf,c=0,D)))^2
    SLC2 <- (sum(ptheta(q[m],af,bf,c=0,D))-sum(ptheta(q[m],at*A,(bt-K)/A,c=0,D)))^2

    Q1 <- Q1 + SLC1*w[m]
    Q2 <- Q2 + SLC2*w[m]

  }
  Q <- Q1 +Q2
  return(Q)
}



#' Equating item parameter only sampled by common item design. FOR ONLY 2PLM
#'
#' @param T_para a parameter file equating TO.
#' @param F_para a parameter file equating FROM.
#' @param method equating method. U can select "SL" or "HB".
#' @param D a factor canstant.
#' @param N the number of nodes.
#' @param mintheta a minimum value of theta in integration.
#' @param maxtheta a maximum value of theta in integration.
#' @param output logical. if TRUE output CSV file.
#' @param Fname character. a file name of csv file.
#' @param Change int. if 0 equated parameter of common item is no transformed parameter from T_para, if 2 transformed parameter, if1 mean and geometric mean.
#' @param Easy logical. if TRUE output file is adaptive for Easy Estimation.
#' @export

CEquating <- function(T_para, F_para, method="SL", D =1.702, N = 31, mintheta=-6, maxtheta = 6,
                      output = F, Fname ="default", Change = 0, Easy = F){

  if(Easy==TRUE){
    Tpara <- T_para[T_para$V3 != 0,]   #remove item that a-parameter is 0
    Fpara <- F_para[F_para$V3 != 0,]
    CII <- Fpara[Fpara$V1 %in% Tpara$V1,1]  #Comon Item Index
    at <- bt <- af <- bf <- numeric(length(CII))

    for (i in 1:length(CII)){
      at[i] <- Tpara[Tpara$V1 == CII[i], "V3"]
      bt[i] <- Tpara[Tpara$V1 == CII[i], "V4"]
      af[i] <- Fpara[Fpara$V1 == CII[i], "V3"]
      bf[i] <- Fpara[Fpara$V1 == CII[i], "V4"]
    }
  } else {
    Tpara <- T_para[T_para$a != 0,]   #remove item that a-parameter is 0
    Fpara <- F_para[F_para$a != 0,]
    CII <- Fpara[Fpara$Item %in% Tpara$Item,1]  #Comon Item Index
    at <- bt <- af <- bf <- numeric(length(CII))

    for (i in 1:length(CII)){
      at[i] <- Tpara[Tpara$Item == CII[i], "a"]
      bt[i] <- Tpara[Tpara$Item == CII[i], "b"]
      af[i] <- Fpara[Fpara$Item == CII[i], "a"]
      bf[i] <- Fpara[Fpara$Item == CII[i], "b"]
    }
  }

  #mean & sigma method
  tm <- mean(bt)  #mean of diff para
  ts <- sqrt(mean((bt-tm)^2))    #standard deviation(It's not used unbiased estimator)
  fm <- mean(bf)
  fs <- sqrt(mean((bf-fm)^2))
  msA <- ts/fs
  msK <- tm - fm * msA

  q <- seq(mintheta,maxtheta,length.out = N)   #nodes of division quadrature
  w <- numeric(N)      #weights of division quadrature
  for (m in 1:N){    w[m] <- dnorm(q[m],0,1)  }
  w <- w/sum(w) # calibration

  if(method == "HB"){
    #criterion function
    res <- stats::nlm(f=cf, p=c(msA,msK),af=af,bf=bf,at=at,bt=bt,D=D,w=w,N=N,q=q)#initial value is equating coefficient estimated by mean & sigma

    A <- res$estimate[1]
    K <- res$estimate[2]

  } else if (method == "SL"){
    res <- stats::nlm(f=cf2, p=c(msA,msK),af=af,bf=bf,at=at,bt=bt,D=D,w=w,N=N,q=q)

    A <- res$estimate[1]
    K <- res$estimate[2]

  }

  if(Easy==TRUE){
    #equating Fpara on the scale of Tpara
    nCII <- numeric(0)
    if(nrow(Fpara) != length(CII)){
      nCII <- Fpara[!(Fpara$V1 %in% Tpara$V1),1]  #extract non comon item
      #common items
      for (n in nCII){
        Fpara[Fpara$V1==n,"V3"] <- round(Fpara[Fpara$V1==n,"V3"]/A, digits = 5)
        Fpara[Fpara$V1==n,"V4"] <- round(Fpara[Fpara$V1==n,"V4"]*A+K, digits = 5)
      }
    }
    #non comon items
    for (i in CII){
      if(Change == 0){
        Fpara[Fpara$V1==i,"V3"] <- Tpara[Tpara$V1==i,"V3"]
        Fpara[Fpara$V1==i,"V4"] <- Tpara[Tpara$V1==i,"V4"]
      } else if(Change == 1){
        Fpara[Fpara$V1 == i,"V3"] <- sqrt(Fpara[Fpara$V1 == i, "V3"]/A * Tpara[Tpara$V1 == i, "V3"]) # geometric mean
        Fpara[Fpara$V1 == i,"V4"] <- (Fpara[Fpara$V1 == i, "V4"]*A+K + Tpara[Tpara$V1 == i, "V4"])/2 # arithmetic mean
      } else if(Change == 2){
        Fpara[Fpara$V1 == i,"V3"] <- Fpara[Fpara$V1 == i, "V3"]/A
        Fpara[Fpara$V1 == i,"V4"] <- Fpara[Fpara$V1 == i, "V4"]*A+K
      }
    }
  }else{
    #equating Fpara on the scale of Tpara
    nCII <- numeric(0)
    if(nrow(Fpara) != length(CII)){
      nCII <- Fpara[!(Fpara$Item %in% Tpara$Item),1]  #extract non comon item
      #common items
      for (n in nCII){
        Fpara[Fpara$Item==n,"a"] <- round(Fpara[Fpara$Item==n,"a"]/A, digits = 5)
        Fpara[Fpara$Item==n,"b"] <- round(Fpara[Fpara$Item==n,"b"]*A+K, digits = 5)
      }
    }
    #non comon items
    for (i in CII){
      if(Change == 0){
        Fpara[Fpara$Item==i,"a"] <- Tpara[Tpara$Item==i,"a"]
        Fpara[Fpara$Item==i,"b"] <- Tpara[Tpara$Item==i,"b"]
      } else if(Change == 1){
        Fpara[Fpara$Item == i,"a"] <- sqrt(Fpara[Fpara$Item == i, "a"]/A * Tpara[Tpara$Item == i, "a"]) # geometric mean
        Fpara[Fpara$Item == i,"b"] <- (Fpara[Fpara$Item == i, "b"]*A+K + Tpara[Tpara$Item == i, "b"])/2 # arithmetic mean
      } else if(Change == 2){
        Fpara[Fpara$Item == i,"a"] <- Fpara[Fpara$Item == i, "a"]/A
        Fpara[Fpara$Item == i,"b"] <- Fpara[Fpara$Item == i, "b"]*A+K
      }
    }
  }


  result <- list(D=paste0("D = ",D), para=Fpara, METHOD=method,
                 EquatingCoefficient_A_K=c(A,K), MeanSigma_A_K=c(msA,msK),ComonItem=CII,NonComonItem=nCII)
  if(output == T){
    utils::write.table(result[1],file=paste0(Fname,"ParaEquated.csv"),sep = "",quote = F, col.names = F, row.names = F)
    utils::write.table(result[2],file=paste0(Fname,"ParaEquated.csv"),sep = ",",quote = F, append = T, col.names = F, row.names = F)
  }
  return(result)

}
