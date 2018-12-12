# general item response model

#' ICC of GIRT model
#'
#' @param theta the location parameter of person
#' @param phi the disprecision parameter of person
#' @param a the discrimination parameter of item
#' @param b the difficulty parameter of item
#' @param D a scale constant
#' @export
#'
gptheta <- function(theta,phi,a,b,D){
  A <- sqrt(1+phi^2*a^2)
  e <- exp(-D*a/A*(theta-b))
  p <- 1/(1+e)
  p
}

# likelihood
gL <- function(u,theta,phi,a,b,D){
  p <- gptheta(theta,phi,a,b,D)
  prod(p^u*(1-p)^(1-u))
}

gLL <- function(u,theta,phi,a,b,D){
  p <- gptheta(theta,phi,a,b,D)
  sum(u*log(p)+(1-u)*log(1-p))
}

#' The density of chi inv dist
#' inverse chi distribution
#'
#' @param phi phi. upper 0.
#' @param v a parameter.
#' @param tau a parameter.
#' @export
#'
dinvchi <- function(phi, v=1, tau=1){
  pow <- v/2
  A <- tau^pow
  B <- 2^(pow-1)*gamma(pow)
  C <- phi^-(v+1)
  D <- exp(-tau/(2*phi^2))
  A/B*C*D
}

sub_rinvx <- function(max, v, tau){
  repeat {
    x <- runif(1,min=0,max=max) # random value
    y <- runif(1) # condition to accept or reject random value
    p <- dinvchi(x, v=v, tau=tau)
    if(p>=y)break
  }
  x
}

#'The random number of chi inv dist
#'
#'@param n the numeber of random variable to generate
#'@param max max ov vector
#'@param v hyper parameter
#'@param tau hyper parameter
#'@export
#'
rinvchi <- function(n, max=5, v=1, tau=1){
  if(is.null(max)) stop("Need a real number to argument 'max' !!")
  if(max <= 0) stop("Need a real number to argument 'max' !!")
  if(n <= 0) stop("Need a positive number to argument 'n' !!")
  x <- matrix(rep(max,n))
  apply(x,1,sub_rinvx,v=v,tau=tau)
}

# M step
# gradient
gr_j <- function(r, N, X, Y, t0, D){
  a <- t0[1]
  b <- t0[2]
  p <- gptheta(X, Y, a, b, D)
  ga <- sum(D*(r-N*p)*(X-b)/sqrt((1+Y^2*a^2)^3))
  gb <- sum(-D*(r-N*p)*a/sqrt(1+Y^2*a^2))
  c(ga,gb)
}

Elnk_j <- function(r,N,t0,X,Y,D){
  a <- t0[1]
  b <- t0[2]
  p <- gptheta(theta=X,phi=Y,a=a,b=b,D=D)
  A <- r*log(p)
  B <- (N-r)*log(1-p)
  sum(A+B)
}

# gradient
grj <- function(r, N, X, Y, a, b, D){
  p <- gptheta(X, Y, a, b, D)
  ga <- sum(D*(r-N*p)*(X-b)/sqrt((1+Y^2*a^2)^3))
  gb <- sum(-D*(r-N*p)*a/sqrt(1+Y^2*a^2))
  c(ga,gb)
}

# Fisher information matrix
Ij <- function(r, N, X, Y, a, b, D){
  p <- gptheta(theta=X, phi=Y, a, b, D)
  ja <- sum(D^2*N*p*(1-p)*(X-b)^2/(1+Y^2*a^2)^3)
  jb <- sum(D^2*N*p*(1-p)*a^2/(1+Y^2*a^2))
  jab <- sum(-D^2*N*p*(1-p)*(X-b)*a/(1+Y^2*a^2)^2)
  matrix(c(ja,jab,jab,jb),nrow = 2)
}

#'Generalized Beta Distribution Beta Distribution in min < x < max
#'
#'@param x A vector consisting of the random variable.
#'@param paramab A vector consisting of a and b parameters.
#'@param rangex A vector consisting of min and max of the random variable.
#'@author Shin-ichi Mayekawa <mayekawa@nifty.com>
#'@export
#'
dgbeta <- function (x, paramab, rangex){
  # copy from lazy.girt package
  # Author: Shin-ichi Mayekawa <mayekawa@nifty.com>
  minscore = rangex[1]
  maxscore = rangex[2]
  a = paramab[1]
  b = paramab[2]
  pdf = stats::dbeta((x - minscore)/(maxscore - minscore), a, b)
  pd = pdf/(maxscore - minscore)^(a + b - 1)
  return(pdf)
}

# functionalize

#'Generalized Item Response Theory parameter estimation
#'
#' @param x DataFrame.
#' @param fc the first column.
#' @param IDc the ID column.
#' @param Ntheta the number of the nodes of theta dist.
#' @param Nphi the number of the nodes of phi dist.
#' @param engine Estep calculation engine.`Cpp` is very faster than `R`.
#' @param method the method of optimiser `optim()` function. Default is "L-BFG-S".
#' @param phi_dist a prior distribution of phi. `invchi` is inverse chi distribution. `lognormal` is log normal distribution.
#' @param v a hyper parameter of invchi for phi
#' @param tau a hyper parameter of invchi for phi
#' @param mu_ph a hyper parameter of lognormal dist for phi
#' @param sigma_ph a hyper parameter of lognormal dist for phi
#' @param min_ph a minimum value of phi distribution
#' @param max_ph a maximum value of phi distribution
#' @param paramab a hyper parameter of generalized beta distribution
#' @param mu_th a hyper parameter of normal dist for theta
#' @param sigma_th a hyper parameter of normal dist for theta
#' @param min_th a minimum value of theta distribution
#' @param max_th a maximum value of theta distribution
#' @param eEM a convergence criterion of item parameters in EM cycle.
#' @param eMLL a convergence criterion of marginal log likelihood in EM cycle.
#' @param maxiter_em the number of iteration of EM cycle.
#'
#' @export
#'
estGip <- function(x, fc=3, IDc=1, Ntheta=31, Nphi=31, engine="Cpp", method="L-BFGS-B",
                   phi_dist = "invchi", v=3, tau=1, mu_ph=0, sigma_ph=0.25, min_ph=0, max_ph=5, paramab=c(1,4),
                   mu_th=0, sigma_th=1, min_th=-4, max_th=4, eEM=0.001, eMLL=0.001, maxiter_em=100){

  if(!(method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent","Fisher_Scoring"))) stop("argument input of `method` is improper string!!")

  X <- as.matrix(x[,fc:ncol(x)])
  nj <- ncol(X)
  nn <- nrow(X)
  nq <- Ntheta
  nr <- Nphi
  if(is.null(IDc)) ID <- 1:nn
  else  ID <- x[,IDc]

  # weight and nodes
  Xq <- seq(min_th, max_th, length.out = nq)
  AX <- dnorm(Xq, mean=mu_th, sd=sigma_th)/sum(dnorm(Xq, mean=mu_th, sd=sigma_th)) # for theta
  Yr <- seq(min_ph, max_ph, length.out = nr)
  if(phi_dist == "invchi"){
    BY <- dinvchi(Yr, v=v, tau=tau)/sum(dinvchi(Yr, v=v, tau=tau)) # for phi
  }else if(phi_dist == "lognormal"){
    BY <- dlnorm(Yr,mu_ph,sigma_ph)/sum(dlnorm(Yr,mu_ph,sigma_ph)) # for phi
  }else if(phi_dist == "uniform"){
    BY <- rep(1/nr,nr)
  }else if(phi_dist == "gbeta"){
    BY <- dgbeta(Yr,paramab = paramab, rangex = c(min_ph,max_ph+1))
    BY <- BY/sum(BY)
  }else {
    stop("argument input of `phi_dist` is improper string!!")
  }
  #
  # initial value
  r <- as.vector(cor(rowSums(X),X))
  pass <- colMeans(X)
  a0 <- 1.702*r/sqrt(1-r^2)
  b0 <- qnorm(pass,0,1, lower.tail = F)/r
  t0 <- data.frame(a=a0,b=b0)

  mll_history <- c(0)

  cat("Estimating Item Parameter!\n")
  t <- 0
  convergence <- T
  while(convergence){
    t <- t + 1
    cat(t,"time EM cycle NOW\n")

    # E step
    knqr <- array(dim = c(nn,nq,nr))
    mll <- 0
    ml <- numeric(nn)
    if(engine=="R"){
      for(r in 1:nr){
        cat(round(100/nr*r, digits = 1),"% / 100%\r")
        for(q in 1:nq){
          l <- apply(X, 1, gL, theta=Xq[q], phi=Yr[r], a=a0, b=b0, D=1.702)*AX[q]*BY[r]
          knqr[,q,r] <- l
          ml <- ml + l
        }
      }
      mll <- sum(log(ml))
    } else if(engine=="Cpp"){
      temp <- Estep_girt(X,a0,b0,Xq,AX,Yr,BY,D=1.702)
      for(i in 1:nn){
        #cat(round(100/nn*i, digits = 1),"% / 100%\r")
        for(q in 1:nq){
          l <- temp$knqr[[i]][[q]]
          knqr[i,q,] <- l
          ml[i] <- ml[i] + sum(l)
        }
      }
      mll <- sum(log(ml))
    }

    cat("-2 Marginal Loglikelihood is",-2*mll,"\n")
    mll_history <- c(mll_history,mll)

    # 規格化
    std_mat <- apply(knqr,1,sum)
    hqr <- sweep(knqr, 1, std_mat, FUN="/")

    rjqr <- array(dim=c(nj,nq,nr))
    Nqr <- array(dim=c(nq,nr))

    for(r in 1:nr){
      for(q in 1:nq){
        rjqr[,q,r] <- t(X)%*%hqr[,q,r]
        Nqr[q,r] <- sum(hqr[,q,r])
      }
    }

    # M step
    # Let Nqr tidyr for FS
    Nqr_dat <- data.frame(Nqr)
    Nqr_dat <- Nqr_dat %>% gather(key=dummy, value=prob)

    t1 <- t0
    for(j in 1:nj){
      cat("Optimising item ",j,"\r")
      # Let Nqr tidyr for FS
      rjqr_dat <- data.frame(rjqr[j,,])
      rjqr_dat <- rjqr_dat %>% gather(key=dummy, value=prob)
      # convert to longer and longer vector
      X_long <- rep(Xq, nr)
      Y_long <- apply(matrix(Yr), 2, rep.int, nq)
      # gradient
      if(method != "Fisher_Scoring"){
        res <- optim(par=c(t0[j,1],t0[j,2]), fn=Elnk_j, gr=gr_j, control = list(fnscale = -1),
                     r=rjqr_dat$prob, N=Nqr_dat$prob, X=X_long, Y=Y_long, D=1.702, method = method)
        t1[j,] <- res$par
      }else{
        # Fisher scoring
        gr <- grj(rjqr_dat$prob, Nqr_dat$prob, X_long, Y_long, t0[j,1], t0[j,2], D=1.702)
        FI <- Ij(rjqr_dat$prob, Nqr_dat$prob, X_long, Y_long, t0[j,1], t0[j,2], D=1.702)
        # solve
        t1[j,] <- t0[j,] + solve(FI)%*%gr
      }
    }

    # convergence check
    a1 <- t1[,1]
    b1 <- t1[,2]
    conv <- c(abs(a0-a1),abs(b0-b1))

    if(all(conv<eEM) || abs(mll_history[t+1]-mll_history[t])<eMLL || maxiter_em == t){
      convergence <- F
      cat("\nConverged!!\n")
      item_para <- as.data.frame(t1)
      break
    }

    t0 <- t1
    a0 <- t0[,1]
    b0 <- t0[,2]

  }

  # Standard Error
  SE <- t1
  Nqr_dat <- data.frame(Nqr)
  Nqr_dat <- Nqr_dat %>% gather(key=dummy, value=prob)

  for(j in 1:nj){
    cat("Optimising item ",j,"\r")
    # Let Nqr tidyr for FS
    rjqr_dat <- data.frame(rjqr[j,,])
    rjqr_dat <- rjqr_dat %>% gather(key=dummy, value=prob)
    # convert to longer and longer vector
    X_long <- rep(Xq, nr)
    Y_long <- apply(matrix(Yr), 2, rep.int, nq)
    # Fisher score matrix
    FI <- Ij(rjqr_dat$prob, Nqr_dat$prob, X_long, Y_long, t1[j,1], t1[j,2], D=1.702)
    # solve
    SE[j,] <- sqrt(diag(solve(FI)))
  }

  cat("Estimating EAP!\n")
  # 個人パラメタの推定
  p1 <- p2 <- numeric(nn)
  for(n in 1:nn){
    cat(round(100/nn*n),"% / 100%\r")
    theta0 <- rep(0,nq)
    phi0 <- rep(0,nr)
    for(r in 1:nr){
      theta0 <- theta0 + hqr[n,,r]
    }
    for(q in 1:nq){
      phi0 <- phi0 + hqr[n,q,]
    }
    p1[n] <- sum(Xq*theta0)
    p2[n] <- sum(Yr*phi0)

  }



  person_para <- data.frame("ID"=ID,"SCORE"=rowSums(X),"theta"=p1, "phi"=p2)

  res <- list(item=item_para, item_SE=SE, person=person_para, hqr=hqr, knqr=knqr)

  return(res)
}
