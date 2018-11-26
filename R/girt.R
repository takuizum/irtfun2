# general item response model

gptheta <- function(theta,phi,a,b,D){
  A <- sqrt(1+phi^2*a^2)
  e <- exp(-D*a/A*(theta-b))
  p <- 1/(1+e)
  p
}
gptheta(1,1,1,0,1)

# likelihood
gL <- function(u,theta,phi,a,b,D){
  p <- gptheta(theta,phi,a,b,D)
  prod(p^u*(1-p)^(1-u))
}

#chi inv dist
#' inverse chi distribution
#'
#' @param phi phi. upper 0.
#' @param v a parameter.
#' @param tau a parameter.
#' @export
#'
invchi <- function(phi, v=3, tau=1){
  pow <- v/2
  A <- tau^pow
  B <- 2^(pow-1)*gamma(pow)
  C <- phi^-(v+1)
  D <- exp(-tau/(2*phi^2))
  A/B*C*D
}

# M step
Elnk_j <- function(rj,N,t0,Xq,Yr,D){
  res <- 0
  a <- t0[1]
  b <- t0[2]
  nq <- length(Xq)
  nr <- length(Yr)
  for(r in 1:nr){
    for(q in 1:nq){
      p <- gptheta(theta=Xq[q],phi=Yr[r],a=a,b=b,D=D)
      A <- rj[q,r]*log(p)
      B <- (N[q,r]-rj[q,r])*log(1-p)
      res <- res + A+B
    }
  }
  res
}

# functionalize

#'Genelar Item Response Theory parameter estimation
#'
#' @param x DataFrame.
#' @param fc the first column.
#' @param IDc the ID column.
#' @param Ntheta the number of the nodes of theta dist.
#' @param Nphi the number of the nodes of phi dist.
#' @param engine Estep calculation engine.`Cpp` is very faster than `R`.
#' @param eEM a convergence criterion of item parameters in EM cycle.
#' @param eMLL a convergence criterion of marginal log likelihood in EM cycle.
#' @param maxiter_em the number of iteration of EM cycle.
#' @export
#'
estGip <- function(x, fc=3, IDc=1, Ntheta=31, Nphi=31, engine="Cpp", eEM=0.001, eMLL=0.001, maxiter_em=100){

  ID <- x[,IDc]
  X <- x[,fc:ncol(x)]
  nj <- ncol(X)
  nn <- nrow(X)
  nq <- Ntheta
  nr <- Nphi

  # weight and nodes
  Xq <- seq(-4,4, length.out = nq)
  AX <- dnorm(Xq)/sum(dnorm(Xq)) # for theta
  Yr <- seq(0.1,5, length.out = nr)
  BY <- invchi(Yr, v=2)/sum(invchi(Yr)) # for phi

  # initial value
  r <- as.vector(cor(rowSums(X),X))
  pass <- colMeans(X)
  a0 <- 1.702*r/sqrt(1-r^2)
  b0 <- qnorm(pass,0,1, lower.tail = F)/r
  t0 <- cbind(a0,b0)

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
    t1 <- t0
    for(j in 1:nj){
      cat("Optimising item ",j,"\r")
      res <- optim(par=c(a0[j],b0[j]), fn=Elnk_j, control = list(fnscale = -1),
                   rj=rjqr[j,,], N=Nqr, Xq=Xq, Yr=Yr, D=1.702)
      t1[j,] <- res$par
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

  cat("Estimating EAP!\n")
  # 個人パラメタの推
  p1 <- p2 <- numeric(nn)
  # EAP
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

  res <- list(item=item_para, person=person_para)

  return(res)
}


#res <- estGip(dat,fc=2, maxiter_em = 2, engine = "R", Ntheta = 10, Nphi=10)
