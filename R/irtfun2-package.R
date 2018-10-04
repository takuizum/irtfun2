#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom stats dnorm
#' @importFrom stats na.omit
#' @importFrom stats optimise
#' @importFrom stats runif
"_PACKAGE"

# P(theta) in two-parameter logisticmodel
ptheta <- function(theta,a,b,c,D=1.702){
  c+(1-c)/(1+exp(-D*a*(theta-b)))
}

# 対数尤度
LL <- function(u,theta,a,b,c,D){
  sum(u*log(ptheta(theta,a,b,c,D))+(1-u)*log(1-ptheta(theta,a,b,c,D)),na.rm = T)
}

# 尤度関数
LL_b <- function(u,theta,a,b,c,mu,sigma,D){
  sum(log(ptheta(theta,a,b,c,D))*u+log(1-ptheta(theta,a,b,c,D))*(1-u)-0.5*((theta-mu)/sigma)^2,na.rm = T)
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

# テスト情報量（尤度関数の二階偏微分の負の期待値）
pitheta <- function(theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  I <- D^2*sum(a^2*(1-p)*(p-c)^2/((1-c)^2*p), na.rm = T)
  1/sqrt(I)
}
pitheta_r <- function(dat,a,b,c,sigma,D){
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
  pi <- pitheta(theta,a,b,c,D)
  L1 <- sum(xi*log(p)+(xi-1)*log(1-p),na.rm = T)
  L2 <- sum(xi*log(p),na.rm=T)
  B <- sqrt(pi)
  L1 + L2 + B
}

# E_3(theta) : expectation
E3 <- function(theta, a, b,c,D){
  sum(
    ptheta(theta,a,b,c,D) * log(ptheta(theta,a,b,c,D)) +
      (1 - ptheta(theta,a,b,c,D)) * log(1- ptheta(theta,a,b,c,D))
    ,na.rm = T
  )
}

# sigma3
S3 <- function(theta, a, b,c,D){
  sqrt(
    sum(
      ptheta(theta,a,b,c,D) * (1 - ptheta(theta,a,b,c,D)) *
        (log(ptheta(theta,a,b,c,D)/(1 - ptheta(theta,a,b,c,D))))^2
      ,na.rm = T
    )
  )
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
Ffg <- function(xi,theta,a,b,c,D){
  exp(sum(xi*log( c+(1-c)/(1+exp(-D*a*(theta-b)))) + (1-xi)*log(1- (c+(1-c)/(1+exp(-D*a*(theta-b))))), na.rm=T))*dnorm(theta,mean=mu,sd=sigma)
}

MI_SE <- function(M){#var() is the function for unbias variance
  K <- length(M)
  (K+1)/K*var(M) + (K-1)/K*var(M)
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

FthetaMAP <- function(xi,a,b,c,mu,sigma,D){

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
      t1 <- t0 - fpdLPD(xi,t0,a,b,c,mu,sigma,D)/spdLPD(t0,a,b,c,sigma,D)
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
              cat("change bisection method. \r")
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

FthetaMLE <- function(xi,a,b,c,D){

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
          cat("change bisection method.\n")
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

FthetaWLE <- function(xi,a,b,c,D){
  #初期値を設定。ログオッズ比を用いた。
  if(sum(xi, na.rm = TRUE) == 0){
    t0 <- log(0.5)
  }else if (sum(xi, na.rm = TRUE) == m){
    t0 <- log(m-0.5)
  }else{
    t0 <- log(sum(xi, na.rm = TRUE)/(m-sum(xi, na.rm = TRUE)))
  }
  opt <- optimise(WL,interval = c(mintheta,maxtheta),maximum = T,xi=xi,a=a,b=b,c=c,D=D)
  t1 <- opt$maximum
  t1
}


#'A estimation theta function.
#'
#'@param xall a item response data
#'@param param a item parameter file.If class is df, parameter column starts from second column and the order is a, b, c.
#'@param est estimation method option.EAP,MAP,MLE,PVs.
#'@param nofrands the number of PVs.
#'@param method iteration method option in MLE and MAP. NR is Newton Rapthon, SANN is Simulated Aannealing, Brent is optimization using optimise function
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
#'@param counter Display progress every counter value.
#'@return a list has ID, rawscore, theta, se, person fit index Z3 and log likelihood.
#'
#'@importFrom stats sd
#'@importFrom stats var
#'@importFrom utils write.csv
#'@export

estheta <- function(xall, param, est="EAP", nofrands=10, method="NR", file="default", output=FALSE, IDc=1, gc=2, fc=3,
                    gh = TRUE, N = 31, D=1.702, maxtheta = 6, mintheta = -6, mu=0, sigma=1, counter=1000){

  #message("データチェック")
  ID <- xall[,IDc]
  if(gc == 0){
    group <- rep(1, nrow(xall))
    G <- 1
    x.all <- xall[,fc:ncol(xall)]
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
      weight <- weight/sum(xnodes)
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
    map_apply <- apply(cbind(group, x.all),1,FthetaMAP,a=a,b=b,c=c,mu=mu,sigma=sigma,D=D)
    #message("MAP推定値の計算が終了しました。")
  }

  #----------------------------------------------------------
  #Eatimate MLE
  #----------------------------------------------------------

  if(est == "MLE"){
    # 全受検者のデータをapply関数で与え，MAP推定を実行
    cat("Calculating MLE.\r")
    mle_apply <- apply(cbind(group, x.all),1,FthetaMLE,a=a,b=b,c=c,D=D)
    #message("MLE推定値の計算が終了しました。")
  }


  #----------------------------------------------------------
  #Eatimate WLE
  #----------------------------------------------------------
  if(est == "WLE"){
    # 全受検者のデータをapply関数で与え，MAP推定を実行
    cat("Calculating WLE.\r")
    wle_apply <- apply(x.all,1,FthetaWLE,a=a,b=b,c=c,D=D)
    #message("WLEの計算が終了しました。")
  }

  if(est == "EAP"){

    # EAP
    eap_m <- eap_apply[1,] %>% round(digits = 5)
    # estimate standard error for EAP estimator
    SE    <- eap_apply[2,] %>% round(digits = 5)
    z3 <- apply(cbind(eap_apply[1,], x.all), 1, pfit, a=a,b=b,c=c,D=D) %>% round(digits = 5)
    lol <-  apply(cbind(eap_apply[1,], x.all), 1, lolF, a=a,b=b,c=c,D=D) %>% round(digits = 5)

    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,EAP=eap_m,SE=SE,z3=z3,lol=lol)
    list("res" = result, "EAPmean&sd" = c(mean(eap_m),sd(eap_m)) %>% round(digits = 5))

  }else if(est == "MAP"){
    #estimate standard error for MAP estimator
    SE <- apply(cbind(map_apply, x.all), 1, pitheta_r, a=a,b=b,c=c,sigma=sigma,D=D) %>% round(digits = 5)
    z3 <- apply(cbind(map_apply, x.all), 1, pfit, a=a,b=b,c=c,D=D) %>% round(digits = 5)
    lol <- apply(cbind(map_apply, x.all), 1, lolF, a=a,b=b,c=c,D=D) %>% round(digits = 5)

    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,MAP=map_apply %>% round(digits = 5), SE=SE, z3=z3, lol=lol)
    list("res" = result, "MAPmean&sd" = c(mean(map_apply), sd(map_apply)) %>% round(digits = 5))

  }else if(est == "WLE"){

    z3 <- apply(cbind(wle_apply, x.all), 1, pfit, a=a,b=b,c=c,D=D) %>% round(digits = 5)
    lol <- apply(cbind(wle_apply, x.all), 1, lolF, a=a,b=b,c=c,D=D) %>% round(digits = 5)
    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,WLE=wle_apply,z3=z3,lol=lol)
    list("res" = result, "WLEmean&sd" = c(mean(wle_apply), sd(wle_apply)))

  }else if(est == "MLE"){
    SE <- apply(matrix(mle_apply, ncol = 1), 1, pitheta, a=a, b=b,c=c,D=D) %>% round(digits = 5)

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
    times <- 0

    cat("Sampling Plausible Values based on von Neumann Rejection sampling.\n")

    for(k in 1:n){
      xi <- x.all[k,]
      times <- times + 1

      #すでに計算してあるEAP推定値と，事後分布の分子を取り出して行列として展開。
      eap   <- eap_apply[1,k]
      const <- eap_apply[3,k]

      #乱数発生時の，P(θ)軸の最大値を設定。
      yheight <- Fmaxpdc(xi,map_apply[k],a,b)*1.001

      # 乱数発生時の，θ軸の最大値を設定。
      zmin <- -maxtheta + eap
      zmax <- mintheta + eap

      nofpv <- 0
      times_sub <- 0
      while( nofpv <= nofrands ){

        y <- runif( 1, 0, yheight)
        z <- runif( 1, zmin, zmax)
        times_sub <- times_sub +1
        fg <- apply(xi,1,Ffg,theta=z,a=a,b=b,c=c,D=D)
        fgvalue <- fg/const

        if( y <= fgvalue){
          nofpv <- nofpv + 1
          if( nofpv > nofrands) break
          pv[k,nofpv] <- z
        }
      }

      if(times==counter){
        cat(k ," / ", n," \r")
        #message(k,"人目の受検者の推算値の発生が終了しました。")
        times <- 0
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
      group_var[i,] <- apply(group_pv,2,var)
      group_sd[i,] <- apply(group_pv,2,sd)
    }

    M_M <- apply(group_mean,1,mean)
    M_V <- apply(group_var,1,mean)
    M_SD <- apply(group_sd,1,mean)

    PS <- data.frame(group=c(1:G),N=ng,mean=M_M, variance=M_V, sd=M_SD)

    SE_M <- apply(group_mean,1,MI_SE)
    SE_V <- apply(group_var,1,MI_SE)
    SE_SD <- apply(group_sd,1,MI_SE)

    SE <- data.frame(group=c(1:G),mean=SE_M, variance=SE_V, sd=SE_SD)

    #message("集団統計量の標準誤差の推定が終了しました。")

    #推算値の平均(Right & Wrong)
    pvmeans <- pvmeans_w <- apply(pv,1,mean)
    pvmeans_r <- apply(pv,2,mean)

    #result

    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,EAP=eap_apply[1,],MAP=map_apply,PVmeans_W=pvmeans_w,PV=pv,AREAS=const)
    plausible_values <-data.frame(ID=ID,Group=group,PV=pv)
    if(output==T){
      #message("結果のCSVファイルを書き出しています。")
      cat("\n output result files.")
      write.csv(result,paste0(file,"_result.csv"),quote=F,row.names = F)
      write.csv(plausible_values,paste0(file,"_PVs.csv"),quote=F,row.names = F)
      write.csv(PS,paste0(file,"_PVS population statistics.csv"),quote=F,row.names = F)
      write.csv(SE,paste0(file,"_PVS standard error.csv"),quote=F,row.names = F)
    }
    pv <- data.frame(ID=ID,group=group,SCORE=xscore,EAP=eap_apply[1,],MAP=map_apply,PV=pv)
    list("PVs" = pv, "EAPmean&sd" = c(mean(eap),sd(eap)),
         "MAPmean&sd" = c(mean(map_apply), sd(map_apply)), "PVmean&sd" = c(M_M, M_SD), "PS" = PS, "SE" = SE)
  }
}
