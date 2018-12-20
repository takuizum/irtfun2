
# the function for item fit index

#'Item Fit indexes 'chi^2' and 'g^2' statistics
#'
#'@param x data.frame of item response data.
#'@param para item parameter data.frame estimated by \code{\link{estip}}.
#'@param theta theta parameter vector. This length must be same to \code{x}.
#'@param fc a first column of item response data.frame.
#'@param H the munber of cut point. To know how to cut, see\code{\link[Hmisc]{cut2}}
#'@param p \emph{p} value.
#'@param D a scale constant
#'@param dot_size a point size of ggplot.See\code{\link[ggpplot2]{geom_point}}
#'@param line_sizea line size of ggplot.See\code{\link[ggpplot2]{geom_linea}}
#'@examples
#'res <- ifind(sim_data_1, sim_param$para, sim_eap$res$EAP, fc=2)
#'# result
#'res$X2
#'res$G2
#'res$ggplot
#'@export
#'
ifind <- function(x, para, theta, fc=3, H = 10 , p = 0.05, D = 1.702, dot_size=1, line_size=1){
  #--------------------------------------------------------#
  # estimate item fit index
  # Yen's Q1
  # X2
  # G2
  #--------------------------------------------------------#
  x.all <- x[,fc:ncol(x)]

  #remove a=0 item parameter and response column
  a <- para$a
  x.all <- x.all[, a != 0]
  if(sum((a == 0)*1) != 0) cat("Remove ",sum((a == 0)*1)," item(s) because a is 0.\n")

  # set item parameter data
  para <- para[para$a != 0, ]
  a <- para$a
  b <- para$b
  c <- para$c
  Item <- para$Item

  # Number of Subjects"
  n <- nrow(x.all)
  # number of items
  m <-length(a)
  xscore <- rowSums(x.all,na.rm=T)

  #result matrix
  ppre <- matrix(nrow=m, ncol=H)
  pobs <- matrix(nrow=m, ncol=H)
  X2 <- data.frame(X2=rep(NA,m), p_value=rep(NA,m), result=rep(NA,m))
  rownames(X2) <- Item
  G2 <- data.frame(G2=rep(NA,m), p_value=rep(NA,m), result=rep(NA,m))
  rownames(G2) <- Item
  thetam <- matrix(nrow=m, ncol=H)

  # predict correct probability
  for(j in 1:m){

    key <- !is.na(x.all[,j])
    THETA <- theta[key]

    # set data
    sub <- Hmisc::cut2(THETA, g = H, levels.mean = TRUE) # 各水準に属する受検者数が等しくなるように，レベル分け
    Om <- table(sub)
    thetam[j,] <- tapply(THETA, sub, mean)

    # observed correct probability
    pobs[j,] <- tapply(x.all[key,j], sub, mean, na.rm = T)

    ppre[j,] <- ptheta(thetam[j,], a[j], b[j], c[j], D=D)

    # X2 statistics
    X2[j, 1] <- sum(Om*(pobs[j,]-ppre[j,])^2/(ppre[j,]*(1-ppre[j,])), na.rm=T)
    X2[j, 2] <- stats::pchisq(X2[j, 1], H-2, lower.tail = FALSE)
    X2[j, 3] <- X2[j, 2] <= p

    # G2 statistics
    G2[j, 1] <- 2*sum(Om*(pobs[j,]*log(pobs[j,]/ppre[j,])+(1-pobs[j,])*log((1-pobs[j,])/(1-ppre[j,]))), na.rm = T)
    G2[j, 2] <- stats::pchisq(G2[j, 1], H-2, lower.tail = FALSE)
    G2[j, 3] <- G2[j, 2] <= p

  }
  ppre <- data.frame(Item=Item, ppre)
  pobs <- data.frame(Item=Item, pobs)
  thetam <- data.frame(Item=Item, thetam)
  res <- list(X2=X2, G2=G2, predict=ppre, observed=pobs, cut=thetam)

  # tidyr for ggplot
  fit_d <- data.frame(res$cut %>% tidyr::gather(key=H1, value=theta, -Item),
                      res$predict[,-1] %>% tidyr::gather(key=H2, value=predict),
                      res$observed[,-1] %>% tidyr::gather(key=H3, value=observed))
  # ggplot
  g_fit <- fit_d %>% ggplot(aes(x=theta, group=Item))+
    ggplot2::geom_point(aes(y=observed), size=dot_size)+
    ggplot2::geom_line(aes(y=observed), size=line_size)+
    # ggplot2::geom_smooth(aes(y=predict, colour=Item), size=line_size,
    #                      method = "glm", method.arg=list(family="binomial"), se=FALSE)+
    ggplot2::geom_line(aes(y=predict, colour=Item), size=line_size)+
    ggplot2::labs(x=TeX("$\\theta$"), y=TeX("$P(\\theta)$"))+
    ggplot2::theme(legend.position = 'none')+
    ggplot2::facet_wrap(~Item, ncol=floor(sqrt(m)))

  res <- list(X2=X2, G2=G2, predict=ppre, observed=pobs, cut=thetam, gg_data=fit_d, ggplot=g_fit)

  return(res)
}



#' Item Fit indexes 'chi^2' and 'g^2' statistics that cut scores based on observed score.
#'
#' @param x data.frame of item response data.
#' @param para item parameter data.frame estimated by \code{\link{estip}}.
#' @param theta theta parameter vector. This length must be same to \code{x}.
#' @param fc a first column of item response data.frame.
#' @param Gc the munber of sub groups.
#' @param p \emph{p} value.
#' @param D a scale constant
#@export # comment out

ifind2 <- function(x, para, theta, Gc=2, fc=3, p = 0.05, D=1.702){
  #--------------------------------------------------------#
  # estimate item fit index
  # S-X2
  # S-G2
  #--------------------------------------------------------#

  if(Gc == 0){
    group <- rep(1, nrow(x))
    G <- 1
  }else{
    group <- x[,Gc]
    G <- max(as.numeric(group))
  }
  x.all <- x[,fc:ncol(x)]

  #remove a=0 item parameter and response column
  a <- para$a
  x.all <- x.all[, a != 0]
  if(sum((a == 0)*1) != 0) cat("Remove ",sum((a == 0)*1)," item(s) because a is 0.\n")

  # set item parameter data
  para <- para[para$a != 0, ]
  a <- para$a
  b <- para$b
  c <- para$c
  Item <- para$Item

  # Number of Subjects"
  ni <- nrow(x.all)
  # number of items
  nj <-length(a)
  xscore <- rowSums(x.all,na.rm=T)

  # combile group and x.all
  x.all <- cbind(group, x.all)

  # 復元得点分布(IRT observed score distribution)による項目適合度の計算

  for(g in 1:G){
    g.all <- x.all[group == g,-1]
    ikey <- colSums(g.all, na.rm=T) != 0
    m <- sum(ikey * 1)
    g.all <- g.all[, ikey]
    ag <- a[ikey]
    bg <- b[ikey]
    cg <- c[ikey]
    tkey <- group == g
    THETA <- as.matrix(theta[tkey])

    # 全項目反応パタンを考慮した得点分布＝復元得点分布を計算
    const <- t(apply(THETA, 1, probability, a=ag, b=bg, c=cg, D=D))
    const <- colSums(const, na.rm = T)

    E <- matrix(nrow=m, ncol=m-1)
    # 行番号が，それぞれ取り除いた項目番号に対応している。
    # 列番号が正答数得点である。したがって1～m-1
    for(j in 1:m){
      # 項目jを取り除いた時の復元得点分布を計算
      cat("item ",j,"\r")
      aj <- ag[j]
      bj <- bg[j]
      cj <- cg[j]
      Tj <- ptheta(THETA, aj, bj, cj, D=D)
      proj <- t(apply(THETA, 1, probability, a=ag[-j], b=bg[-j], c=cg[-j], D=D))
      E[j,] <- t(Tj) %*% proj[, -ncol(proj)] # 満点の分は計算しない
    }

    O <- matrix(nrow=m, ncol=m-1)
    Nj <- numeric(m-1)
    score <- rowSums(g.all, na.rm=T)
    dat <- cbind(score, g.all)
    for(s in 1:(m-1)){
      dat2 <- dat[dat[,1]==s,-1]
      Nj[s] <- nrow(dat2)
      O[,s] <- colSums(dat2, na.rm = T)/Nj[s]
    }

    SX2 <- matrix(nrow=m-1, ncol=3)
    SG2 <- matrix(nrow=m-1, ncol=3)

    for(j in 1:(m-1)){
      SX2[j,1] <- sum((Nj[j]*O[j,] - E[j,]/const[j])^2/(E[j,]/const[j]*(1-E[j,]/const[j])), na.rm = T)
      SX2[j, 2] <- stats::pchisq(SX2[j, 1], m-2, lower.tail = FALSE)
      SX2[j, 3] <- SX2[j, 2] <= p

      SG2[j,1] <- 2*sum(Nj*(O[j,]*log(O[j,]/E[j,]/const[j]+(1-O[j,])*log((1-O[j,])/(1-E[j,]/const[j])))), na.rm = T)
      SG2[j, 2] <- stats::pchisq(SG2[j, 1], m-2, lower.tail = FALSE)
      SG2[j, 3] <- SG2[j, 2] <= p
    }

    #----------------------------------------------------------------#

  }

}

#'Item Fit indexes 'In Fit' and 'Out Fit' statistics, which is based on expexted and observed probability.
#'
#'@param x data.frame of item response data.
#'@param para item parameter data.frame estimated by \code{\link{estip}}.
#'@param theta theta parameter vector. This length must be same to \code{x}.
#'@param fc a first column of item response data.frame.
#'@param D a scale constant
#'@examples
#'res <- ifind3(sim_data_1, sim_param$para, sim_eap$res$EAP, fc=2)
#'res
#'@export
#'

ifind3 <- function(x, para, theta, fc=3, D=1.702){
  #--------------------------------------------------------#
  # estimate item fit index
  # In Fit & Out Fit
  #--------------------------------------------------------#

  x.all <- x[,fc:ncol(x)]
  #remove a=0 item parameter and response column
  a <- para$a
  x.all <- x.all[, a != 0]
  message()
  if(sum((a == 0)*1) != 0) cat("Remove ",sum((a == 0)*1)," item(s) because a is 0.\n")

  # set item parameter data
  para <- para[para$a != 0, ]
  a <- para$a
  b <- para$b
  c <- para$c
  Item <- para$Item

  # Number of Subjects"
  n <- nrow(x.all)
  # number of items
  m <-length(a)
  xscore <- rowSums(x.all,na.rm=T)

  InFit <- numeric(m)
  OutFit <- numeric(m)
  Fit <- data.frame(Item=Item, N=rep(NA,m), InFit=rep(NA,m), StdInFit=rep(NA,m),
                    OutFit=rep(NA,m), StdOutFit=rep(NA,m), F0.025p=rep(NA,m), F0.975p=rep(NA,m))

  # predict correct probability
  for(j in 1:m){

    key <- !is.na(x.all[,j]) # 回答している受検者のデータだけ抽出
    u <- x.all[key,j]
    THETA <- theta[key]
    N <- length(u)
    Fit$N[j] <- N
    Fit$F0.025p[j] <- qf(0.025, N-1, Inf)
    Fit$F0.975p[j] <- qf(0.975, N-1, Inf)
    # set data
    p <- ptheta(THETA,a[j],b[j],c[j],D=D)
    w <- p*(1-p)
    mw <- w*(1-3*w)
    q <- sqrt(sum(mw-w^2))/sum(w)
    Fit$OutFit[j] <- sum((u-p)^2/(p*(1-p)))/(N-1)
    Fit$StdOutFit[j] <- (log(Fit$OutFit[j])+Fit$OutFit[j]-1)*(sqrt((N-1)/8))
    Fit$InFit[j] <- sum((u-p)^2)/sum(p*(1-p))
    Fit$StdInFit[j] <- (3/q)*(Fit$InFit[j]^(1/3)-1)+q/3

  }
  cat("A guide: OutFit follows F distribution that DF is N-1.\nStdOutFit follows Standard Normal distribution.\n")
  cat("A guide of InFit is 0.75~1.3, StdInFit is -2.0~2.0.\n")
  return(Fit)
}

r2df <- function(r,x,y){
  r[rownames(r)==x, colnames(r)==y]
}
nfordf <- function(x, X, Y){
  a <- x[,colnames(x)==X]
  b <- x[,colnames(x)==Y]
  na.omit(a+b) %>% length()
}


#'Item Fit indexes 'Q3'statistics. Diagnosis Local Item Independence(LID)
#'
#'@param x data.frame of item response data.
#'@param para item parameter data.frame estimated by \code{\link{estip}}.
#'@param theta theta parameter vector. This length must be same to \code{x}.
#'@param fc a first column of item response data.frame.
#'@param D a scale constant
#'@param use an optional character string giving a method for computing covariances in the presence of missing values. See \code{\link{cor}}
#'@param abs_value an absolute value for Q3.
#'
#'@export
#'
Q3_stat <- function(x, para, theta, fc=3, D = 1.702, use = "pairwise.complete.obs", abs_value=0.2){
  # item response data
  x.all <- x[,fc:ncol(x)]

  #remove a=0 item parameter and response column
  a <- para$a
  x.all <- x.all[, a != 0]
  message()
  if(sum((a == 0)*1) != 0) cat("Remove ",sum((a == 0)*1)," item(s) because a is 0.\n")

  # set item parameter data
  para <- para[para$a != 0, ]
  a <- para$a
  b <- para$b
  c <- para$c
  Item <- para$Item %>% as.character()

  # Number of Subjects"
  n <- nrow(x.all)
  # number of items
  m <-length(a)

  # expected test score
  p <- apply(matrix(theta), 1, ptheta, a=a, b=b, c=c, D=D) %>% t()
  # residual score
  d <- x.all-p

  # correlation of residual matrix
  R <- cor(d, use=use)

  # R matrix for data.frame
  X <- apply(matrix(Item), 1, rep.int, times=m) %>% as.vector()
  Y <- rep(Item, m)
  r <- mapply(r2df, matrix(X), matrix(Y), MoreArgs = list(r=R)) %>% as.vector()
  R_df <- data.frame(X=X[!is.na(r)], Y=Y[!is.na(r)], r=r[!is.na(r)])
  R_df$abs <- R_df$r > abs_value

  #R_df$N <- mapply(nfordf, matrix(R_df$X), matrix(R_df$Y), MoreArgs = list(x=x)) %>% as.vector()
  #R_df$z <- atanh(R_df$r)
  #R_df$p0.025 <- pnorm(0.025, mean=atanh(R_df$r)+5*R_df$r/(2*(R_df$N-2)), sd=sqrt(1/(R_df$N-1)))
  #R_df$p0.975 <- pnorm(0.975, mean=atanh(R_df$r)+5*R_df$r/(2*(R_df$N-2)), sd=sqrt(1/(R_df$N-1)))
  R_df <- R_df[R_df$X!=R_df$Y, ]

  list(df=R_df, matrix=R)
  #R_df
}
