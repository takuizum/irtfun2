
# M step
# gradient
gr_jm <- function(r, N, X, t0, D){
  if(length(t0)==2){
    a <- t0[1]
    b <- t0[2]
    c <- 0
    p <- ptheta(X, a, b, c, D)
    ga <- sum(D*(X-b)*(r-N*p)*(p-c)/((1-c)*p))
    gb <- sum(-D*a*(r-N*p)*(p-c)/((1-c)*p))
    res <- c(ga,gb)
  }
  if(length(t0)==3){
    a <- t0[1]
    b <- t0[2]
    c <- t0[3]
    p <- ptheta(X, a, b, c, D)
    ga <- sum(D*(X-b)*(r-N*p)*(p-c)/((1-c)*p))
    gb <- sum(-D*a*(r-N*p)*(p-c)/((1-c)*p))
    gc <- sum((r-N*p)/((1-c)*p))
    res <- c(ga,gb,gc)
  }
  res
}

Elnk_jm <- function(r,N,t0,X,D){
  if(length(t0)==2){
    a <- t0[1]
    b <- t0[2]
    c <- 0
    p <- ptheta(X, a, b, c, D)
    A <- suppressWarnings(r*log(p)+(N-r)*log(1-p))
  }
  if(length(t0)==3){
    a <- t0[1]
    b <- t0[2]
    c <- t0[3]
    p <- ptheta(X, a, b, c, D)
    A <- suppressWarnings(r*log(p)+(N-r)*log(1-p))
  }
  sum(A)
}

Elnk_jm_1pl <- function(r,N,t0,X,D,fix_a){
  a <- fix_a
  b <- t0
  c <- 0
  p <- ptheta(X, a, b, c, D)
  A <- r*log(p)
  B <- (N-r)*log(1-p)
  sum(A+B)
}

# gradient for marginal Bayes
gr_jm_b <- function(r, N, X, t0, D, meanlog, sdlog, mean, sd, shape1, shape2){
  if(length(t0)==2){
    a <- t0[1]
    b <- t0[2]
    c <- 0
    p <- ptheta(X, a, b, c, D)
    ga <- sum(D*(X-b)*(r-N*p)*(p-c)/((1-c)*p))-1/a-(log(a)-meanlog)/(a*sdlog^2)
    gb <- sum(-D*a*(r-N*p)*(p-c)/((1-c)*p))-(b-mean)/sd^2
    res <- c(ga,gb)
  }
  if(length(t0)==3){
    a <- t0[1]
    b <- t0[2]
    c <- t0[3]
    p <- ptheta(X, a, b, c, D)
    ga <- sum(D*(X-b)*(r-N*p)*(p-c)/((1-c)*p))-1/a-(log(a)-meanlog)/(a*sdlog^2)
    gb <- sum(-D*a*(r-N*p)*(p-c)/((1-c)*p))-(b-mean)/sd^2
    gc <- sum((r-N*p)/((1-c)*p))+(shape1-2)/c-(shape2-2)/(1-c)
    res <- c(ga,gb,gc)
  }
  res
}

Elnk_jm_b <- function(r,N,t0,X,D,meanlog,sdlog,mean,sd,shape1,shape2){
  if(length(t0)==2){
    a <- t0[1]
    b <- t0[2]
    c <- 0
    p <- ptheta(X, a, b, c, D)
    A <- suppressWarnings(r*log(p)+(N-r)*log(1-p))
    res <- sum(A)+dnorm(b,mean,sd)+dlnorm(a,meanlog,sdlog)
  }
  if(length(t0)==3){
    a <- t0[1]
    b <- t0[2]
    c <- t0[3]
    p <- ptheta(X, a, b, c, D)
    A <- suppressWarnings(r*log(p)+(N-r)*log(1-p))
    res <- sum(A)+dnorm(b,mean,sd)+dlnorm(a,meanlog,sdlog)+dbeta(c,shape1,shape2)
  }
  res
}

Elnk_jm_1pl_b <- function(r,N,t0,X,D,fix_a,mean,sd){
  a <- fix_a
  b <- t0
  c <- 0
  p <- ptheta(X, a, b, c, D)
  A <- r*log(p)
  B <- (N-r)*log(1-p)
  sum(A+B)+dnorm(b,mean,sd)
}

grjm_b <- function(r, N, X, a, b, c, D, model, meanlog, sdlog, mean, sd, shape1, shape2){
  p <- ptheta(X, a, b, c, D)
  if(model=="1PL"){
    gb <- sum(-D*a*(r-N*p)*(p-c)/((1-c)*p))-(b-mean)/sd^2
    res <- gb
  }
  if(model=="2PL"){
    ga <- sum(D*(X-b)*(r-N*p)*(p-c)/((1-c)*p))-1/a-(log(a)-meanlog)/(a*sdlog^2)
    gb <- sum(-D*a*(r-N*p)*(p-c)/((1-c)*p))-(b-mean)/sd^2
    res <- c(ga,gb)
  }
  if(model=="3PL"){
    ga <- sum(D*(X-b)*(r-N*p)*(p-c)/((1-c)*p))-1/a-(log(a)-meanlog)/(a*sdlog^2)
    gb <- sum(-D*a*(r-N*p)*(p-c)/((1-c)*p))-(b-mean)/sd^2
    gc <- sum((r-N*p)/((1-c)*p))+(shape1-2)/c-(shape2-2)/(1-c)
    res <- c(ga,gb,gc)
  }
  res
}

# objective function for regularized marginal maximum likelihood estimation
Elnk_jm_r <- function(r,N,t0,X,D,penalty){
  if(length(t0)==2){
    a <- t0[1]
    b <- t0[2]
    c <- 0
    p <- ptheta(X, a, b, c, D)
    A <- suppressWarnings(r*log(p)+(N-r)*log(1-p))
    res <- sum(A)-penalty
  }
  if(length(t0)==3){
    a <- t0[1]
    b <- t0[2]
    c <- t0[3]
    p <- ptheta(X, a, b, c, D)
    A <- suppressWarnings(r*log(p)+(N-r)*log(1-p))
    res <- sum(A)-penalty
  }
  res
}

#
en_penalty <- function(a, alpha=0.5, lambda=0.1){
  penalty <- alpha*sum(abs(a))+(1-alpha)/2*sum(a^2)
  penalty
}

# gradient for FS
grjm <- function(r, N, X, a, b, c, D, model){
  p <- ptheta(X, a, b, c, D)
  if(model=="1PL"){
    gb <- sum(-D*a*(r-N*p)*(p-c)/((1-c)*p))
    res <- gb
  }
  if(model=="2PL"){
    ga <- sum(D*(X-b)*(r-N*p)*(p-c)/((1-c)*p))
    gb <- sum(-D*a*(r-N*p)*(p-c)/((1-c)*p))
    res <- c(ga,gb)
  }
  if(model=="3PL"){
    ga <- sum(D*(X-b)*(r-N*p)*(p-c)/((1-c)*p))
    gb <- sum(-D*a*(r-N*p)*(p-c)/((1-c)*p))
    gc <- sum((r-N*p)/((1-c)*p))
    res <- c(ga,gb,gc)
  }
  res
}

# Fisher information matrix
Ijm <- function(N, X, a, b, c, D, model){
  p <- ptheta(X, a, b, c, D)
  q <- 1-p
  if(model=="1PL"){
    jb <- sum(N*(p-c)^2*q*a^2*D^2/((1-c)^2*p))
    res <- jb
  }
  if(model=="2PL"){
    ja <- sum(N*(p-c)^2*q*(X-b)^2*D^2/((1-c)^2*p))
    jb <- sum(N*(p-c)^2*q*a^2*D^2/((1-c)^2*p))
    jab <- sum(-N*(p-c)^2*q*a*(X-b)*D^2/((1-c)^2*p))
    res <- matrix(c(ja,jab,jab,jb),nrow = 2)
  }
  if(model=="3PL"){
    ja <- sum(N*(p-c)^2*q*(X-b)^2*D^2/((1-c)^2*p))
    jb <- sum(N*(p-c)^2*q*a^2*D^2/((1-c)^2*p))
    jc <- sum(N*q/(p*(1-c)^2))
    jab <- sum(-N*(p-c)^2*q*a*(X-b)*D^2/((1-c)^2*p))
    jac <- sum(N*(p-c)*q*(X-b)*D/((1-c)^2*p))
    jbc <- sum(-N*(p-c)*q*a*D/((1-c)^2*p))
    res <- matrix(c(ja,jab,jac,
                    jab,jb,jbc,
                    jac,jbc,jc),nrow = 3)
  }
  res
}



#' Unidimensional Item Response Theory parameter estimation
#'
#' @param x DataFrame.
#' @param fc the first column.
#' @param IDc the ID column.
#' @param Gc the grade column.
#' @param bg a mumber of base grade.
#' @param D a scale constant.
#' @param Ntheta the number of the nodes of theta dist.
#' @param method the method of optimiser.  Default is "Fisher_Scoring", but \code{\link[stats]{optim}} function also be able to use.
#' @param model a model vector
#' @param max_func a character of object function. "N" is MML-EM, "B" is marginal Bayes and "R" is Regularized MML.
#' @param mu_th a hyper parameter of normal dist for theta
#' @param sigma_th a hyper parameter of normal dist for theta
#' @param min_th a minimum value of theta distribution
#' @param max_th a maximum value of theta distribution
#' @param eEM a convergence criterion related to item parameters in EM cycle.
#' @param eMLL a convergence criterion related to negative twice log likelihood in EM cycle.
#' @param maxiter_em the number of iteration of EM cycle.
#' @param rm_list a vector of item U want to remove for estimation. NOT list.
#' @param fix_a a fix parameter for slope parameter of 1PLM
#' @param mu_a a hyper parameter for slope parameter prior distribution(lognormal) in marginal Bayes estimation.
#' @param sigma_a a hyper parameter for slope parameter prior distribution(lognormal) in marginal Bayes estimation.
#' @param mu_b a hyper parameter for location parameter prior distribution(normal) in marginal Bayes estimation.
#' @param sigma_b a hyper parameter for location parameter prior distribution(normal) in marginal Bayes estimation.
#' @param mu_c a hyper parameter for lower asymptote parameter prior distribution(beta) in marginal Bayes estimation.
#' @param sigma_c a hyper parameter for lower asymptote parameter prior distribution(beta) in marginal Bayes estimation.
#' @param w_c a weight parameter for lower asymptote parameter prior distribution(beta) in marginal Bayes estimation.
#' @param alpha tuning parameter of elastic net penalty.
#' @param lambda tuning parameter of elastic net penalty.
#' @param th_dist a type of theta dist."normal" or "empirical"
#' @param print How much information you want to display? from 1 to 3. The larger, more information is displayed.
#'
#' @return the output is a list that has item parameters data.frame and these Standard Error.
#' @examples
#' # MMLE
#' res1 <- estip2(x=sim_data_2, Gc=NULL, fc=2, Ntheta=21)
#' # check the parameters
#' res1$para
#' res1$SE
#'
#' # Multigroup MMLE
#' res2 <- estip2(x=sim_dat_st, Gc=2, bg=3, fc=3, Ntheta=10)
#'
#' # Marginal Bayes
#' res3 <- estip2(x=sim_data_2, Gc=NULL, fc=2, Ntheta=21, max_func="B")
#'
#' # Regularized MMLE
#' res4 <- estip2(x=sim_data_2, Gc=NULL, fc=2, Ntheta=21, max_func="R")
#' @export
#'
estip2 <- function(x, fc=3, IDc=1, Gc=NULL, bg=1, Ntheta=31, D=1.702, method="Fisher_Scoring", model="2PL", max_func="N", rm_list=NULL,
                   mu_th=0, sigma_th=1, min_th=-4, max_th=4, eEM=0.001, eMLL=0.001, maxiter_em=100, th_dist="normal",
                   fix_a=1, mu_a=0, sigma_a=1, mu_b=0, sigma_b=1, mu_c=0, sigma_c=1, w_c=1, alpha=0.5, lambda=1,
                   print=0){

  if(!(method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Fisher_Scoring"))) stop("argument input of 'method' is improper string!!")
  if(!(max_func %in% c("B","N","R"))) stop("argument input of 'max_func' is improper string!!")
  if(max_func == "R"){
    if(method == "Fisher_Scoring"){
      cat("'Fisher scoring' cannot be selected as optimization method for Regularized MML !!\nChange 'BFGS'\n")
      method <- "BFGS"
    }
    if(model == "1PL") stop("In 1PLM, Regularized MML cannot be conducted.")
  }
  if(th_dist %in% c("normal", "empirical"))

  # data check
  X <- as.matrix(x[,fc:ncol(x)])
  Item <- colnames(X) %>% as.character()
  nj <- ncol(X)
  ni <- nrow(X)
  nq <- Ntheta
  if(is.null(IDc)) ID <- 1:ni
  else  ID <- x[,IDc]
  if(is.null(Gc)) group <- rep(1,ni)
  else group <- x[,Gc]
  ng <- max(group)
  if(!is.numeric(group)) stop("A group column must be numeric, but has non numeric scalar.")
  if(min(group)==0) stop("A group column must be above 1, but has 0.")
  if(length(model)!=nj) model <- rep(model, nj)
  p_mean <- rep(mu_th, ng)
  p_sd <- rep(sigma_th,ng)

  # design matrix
  ind <- matrix(nrow=ng, ncol=nj)
  for(g in 1:ng){
    ind[g,] <- X[group==g, ] %>% colSums(na.rm = T)
  }
  ind[ind!=0] <- 1
  resp <- X %>% as.matrix()
  resp[!is.na(resp)] <- 1
  resp[is.na(resp)] <- 0

  # weight and nodes
  Xq <- seq(min_th, max_th, length.out = nq)
  # initial weight
  # AX <- matrix(nrow=nq, ncol=ng)
  AX <- dnorm(Xq, mean=mu_th, sd=sigma_th)/sum(dnorm(Xq, mean=mu_th, sd=sigma_th))
  AX <- AX %>% rep.int(times=ng) %>% matrix(ncol=ng)

  # initial value
  r <- as.vector(cor(rowSums(X, na.rm = T),X, use = "pair"))
  pass <- colMeans(X, na.rm = T)
  a0 <- D*r/sqrt(1-r^2)
  b0 <- -log(pass/(1-pass))
  c0 <- rep(0,nj)
  for(j in 1:nj){
    if(model[j]=="3PL")c0 <- 0.2
    if(model[j]=="1PL")a0 <- fix_a
  }
  init <- t0 <- t1 <- data.frame(a=a0,b=b0,c=c0)

  # remove selected item
  if(!is.null(rm_list)){
    rm_ind <- c(1:nj)[Item %in% rm_list]
    model[rm_ind] <- "NONE"
    resp[,rm_ind] <- 0
    ind[,rm_ind] <- 0
    init[rm_ind,] <- t0[rm_ind,] <- t1[rm_ind,] <- 0
    a0 <- init$a
    b0 <- init$b
    c0 <- init$c
    cat("Remove Item ", rm_list, "\n")
  }

  # cat
  cat("The number of subject is " ,ni, ".\nThe number of item is ", nj-length(rm_list),
      ".\nThe number of remove item is ", length(rm_list), ".\n")

  # for imputation
  a1 <- b1 <- c1 <- numeric(nj)

  # reshape beta dist hyper parameter
  shape1 <- w_c*mu_c+1
  shape2 <- w_c*(1-mu_c)+1

  # set mll history vector
  mll_history <- c(0)

  cat("Estimating Item Parameter!\n")

  # stand FLAG
  t <- 0
  convergence <- T
  while(convergence){
    t <- t + 1
    #cat(t,"time EM cycle NOW\n")

    # E step
    Estep <- Estep_irt(xall=X, t0=as.matrix(t0), Xm=Xq, Wm=AX, group=group,
                       ind=ind, resp=resp, D=D, MLL=mll_history)

    mll <- Estep$MLL[t+1]
    if(print == 0) cat(t ,"times -2 Marginal Loglikelihood is",-2*mll,"\r")
    if(print >= 1) cat(t ,"times -2 Marginal Loglikelihood is",-2*mll,"\n")
    mll_history <- Estep$MLL

    Njm <- matrix(0,nrow=nj, ncol=nq)
    rjm <- matrix(0,nrow=nj, ncol=nq)
    for(g in 1:ng){
      # purrr::map(Estep$Njm,function(x){purrr::invoke(.f=rbind, .x=x)})
      Njm <- Njm + purrr::invoke(rbind,Estep$Njm[[g]])
      rjm <- rjm + purrr::invoke(rbind,Estep$rjm[[g]])
    }

    # M step
    a0 <- t0$a
    b0 <- t0$b
    c0 <- t0$c

    # penalty
    penalty <- en_penalty(a0, alpha = alpha, lambda = lambda)
    for(j in 1:nj){
      if(model[j]=="NONE") next
      #cat("Optimising item ",j,"\r")
      # Let Nqr tidyr for FS
      N <- Njm[j,]
      rr <- rjm[j,]
      # normal M step
      if(max_func == "N"){
        if(method != "Fisher_Scoring"){
          #stop("This option dose not work now!!")
          if(model[j]=="1PL"){
            sol <- optimise(Elnk_jm_1pl, interval = c(min_th,max_th)+t0$b[j], maximum = T,
                            r=rr, N=N, X=Xq, D=D, fix_a=fix_a)
            t1[j,2] <- sol$maximum
          }
          if(model[j]=="2PL"){
            sol <- optim(par=c(t0$a[j],t0$b[j]), fn=Elnk_jm, gr=gr_jm, control = list(fnscale = -1),
                         r=rr, N=N, X=Xq, D=D, method = method)
            t1[j,c(1,2)] <- sol$par
          }
          if(model[j]=="3PL"){
            sol <- optim(par=c(t0$a[j],t0$b[j],t0$c[j]), fn=Elnk_jm, gr=gr_jm, control = list(fnscale = -1),
                         r=rr, N=N, X=Xq, D=D, method = method)
            t1[j,c(1,2,3)] <- sol$par
          }
        }else{
          # Fisher scoring
          gr <- grjm(rr, N, Xq, t0$a[j], t0$b[j], t0$c[j], D=D, model=model[j])
          FI <- Ijm(N, Xq, t0$a[j], t0$b[j], t0$c[j], D=D, model=model[j])
          # solve
          sol <- solve(FI)%*%gr
          if(model[j]=="1PL") t1$b[j] <- t0$b[j]+sol
          if(model[j]=="2PL") t1[j,c(1,2)] <- t0[j,c(1,2)]+sol
          if(model[j]=="3PL") t1[j,] <- t0[j,]+sol
        }
      } else if (max_func=="B"){
        if(method != "Fisher_Scoring"){
          #stop("This option dose not work now!!")
          if(model[j]=="1PL"){
            sol <- optimise(Elnk_jm_1pl_b, interval = c(min_th,max_th)+t0$b[j], maximum = T,
                            r=rr, N=N, X=Xq, D=D, fix_a=fix_a, mean=mu_b, sd=sigma_b)
            t1[j,2] <- sol$maximum
          }
          if(model[j]=="2PL"){
            sol <- optim(par=c(t0$a[j],t0$b[j]), fn=Elnk_jm_b, gr=gr_jm_b, control = list(fnscale = -1),
                         r=rr, N=N, X=Xq, D=D, method = method,
                         meanlog=mu_a, sdlog=sigma_a, mean=mu_b, sd=sigma_b, shape1=shape1, shape2=shape2)
            t1[j,c(1,2)] <- sol$par
          }
          if(model[j]=="3PL"){
            sol <- optim(par=c(t0$a[j],t0$b[j],t0$c[j]), fn=Elnk_jm_b, gr=gr_jm_b, control = list(fnscale = -1),
                         r=rr, N=N, X=Xq, D=D, method = method,
                         meanlog=mu_a, sdlog=sigma_a, mean=mu_b, sd=sigma_b, shape1=shape1, shape2=shape2)
            t1[j,c(1,2,3)] <- sol$par
          }
        }else{
          # Fisher scoring
          gr <- grjm_b(rr, N, Xq, t0$a[j], t0$b[j], t0$c[j], D=D, model=model[j],
                       meanlog=mu_a, sdlog=sigma_a, mean=mu_b, sd=sigma_b, shape1=shape1, shape2=shape2)
          FI <- Ijm(N, Xq, t0$a[j], t0$b[j], t0$c[j], D=D, model=model[j])
          # solve
          sol <- solve(FI)%*%gr
          if(model[j]=="1PL") t1$b[j] <- t0$b[j]+sol
          if(model[j]=="2PL") t1[j,c(1,2)] <- t0[j,c(1,2)]+sol
          if(model[j]=="3PL") t1[j,] <- t0[j,]+sol
        }
      }else if(max_func=="R"){
        if(model[j]=="2PL"){
          sol <- optim(par=c(t0$a[j],t0$b[j]), fn=Elnk_jm_r, control = list(fnscale = -1),
                       r=rr, N=N, X=Xq, D=D, method = method, penalty = penalty)
          t1[j,c(1,2)] <- sol$par
        }
        if(model[j]=="3PL"){
          sol <- optim(par=c(t0$a[j],t0$b[j],t0$c[j]), fn=Elnk_jm_r, control = list(fnscale = -1),
                       r=rr, N=N, X=Xq, D=D, method = method, penalty = penalty)
          t1[j,c(1,2,3)] <- sol$par
        }
      }

      # exception handling
      if(t1$c[j] < 0){
        warning(paste0("The model of item ",j," was changed to 2PLM"))
        model[j] <- "2PL"
        t1[j,1] <- init$a[j]
        t1[j,2] <- init$b[j]
        t1[j,3] <- 0
      }
      if(t1$a[j] < 0){
        warning(paste0("In ",t," times iteration, Item ",j," was eliminated. 'a' must be positive, but was negative."))
        model[j] <- "NONE"
        t1[j,] <- t0[j,]
      }
    }

    # calibration
    # calculate mean and sd
    Njm <- matrix(0,nrow=nj, ncol=nq)
    p_mean_t0 <- p_mean
    p_sd_t0 <- p_sd
    for(g in 1:ng){
      gind <- ind[g,]!=0
      # purrr::map(Estep$Njm,function(x){purrr::invoke(.f=rbind, .x=x)})
      Njm <- purrr::invoke(rbind, Estep$Njm[[g]])
      #rjm <- rjm + purrr::invoke(rbind,Estep$rjm[[g]])
      p_mean[g] <- mean(Njm[gind,]%*%matrix(Xq)/rowSums(Njm[gind,]))
      p_sd[g] <- mean(sqrt(Njm[gind,]%*%matrix((Xq-p_mean[g])^2)/rowSums(Njm[gind,])))
    }
    # calculate calibration weight
    A <- sigma_th/p_sd[bg]
    K <- mu_th - A*p_mean[bg]
    # calibrate mean and sd
    p_mean <- A*p_mean+K
    p_sd <- A*p_sd
    for(g in 1:ng){
      gind <- ind[g,]!=0
      Njm <- purrr::invoke(rbind, Estep$Njm[[g]])
      if(th_dist=="normal"){
        # Gaussian dist
        AX[,g] <- dnorm(Xq, mean=p_mean[g], sd=p_sd[g])/sum(dnorm(Xq, mean=p_mean[g], sd=p_sd[g]))
      }
      if(th_dist=="empirical"){
        # empirical dist
        constNjm <- rowSums(Njm[gind,]) %>% rep.int(times=nq) %>% matrix(ncol=nq)
        AX[,g] <- (Njm[gind,]/constNjm) %>% colMeans(na.rm = T)
      }
    }
    # calibrate item parameter
    for(j in 1:nj){
      if(model[j]=="NONE") next
      if(model[j] != "1PL"){
        a1[j] <- t1$a[j]/A
        b1[j] <- t1$b[j]*A+K
        c1[j] <- t1$c[j]
      }else{
        a1[j] <- fix_a
        b1[j] <- t1$b[j]*A+K
        c1[j] <- t1$c[j]
      }
    }
    t1$a <- a1
    t1$b <- b1
    t1$c <- c1

    # convergence check
    conv1 <- c(abs(a0-a1),abs(b0-b1),abs(c0-c1))
    conv4 <- abs(mll_history[t+1]-mll_history[t])
    if(ng > 1){
      conv2 <- c(abs(p_mean - p_mean_t0))
      conv3 <- c(abs(p_sd - p_sd_t0))
    }else{
      conv2 <- conv3 <- 1
    }

    if(all(conv1<eEM) || conv4<eMLL || all(conv2<eEM) || all(conv3<eEM) || maxiter_em == t){
      convergence <- F
      cat("\nConverged!!\n")
      item_para <- as.data.frame(t1)
      break
    }else{
      if(print >= 1){
        cat("Item maximum changed a is", Item[max(abs(a0-a1))==abs(a0-a1)],"=",max(abs(a0-a1)),"\n")
        cat("Item maximum changed b is", Item[max(abs(b0-b1))==abs(b0-b1)],"=",max(abs(b0-b1)),"\n")
        cat("Item maximum changed c is", Item[max(abs(c0-c1))==abs(c0-c1)],"=",max(abs(c0-c1)),"\n")
      }
      t0 <- t1
    }
  }

  # LAST E step
  # E step
  Estep <- Estep_irt(xall=X, t0=as.matrix(t1), Xm=Xq, Wm=AX, group=group,
                     ind=ind, resp=resp, D=D, MLL=mll_history)

  mll <- Estep$MLL[t+1]
  mll_history <- Estep$MLL

  Njm <- matrix(0,nrow=nj, ncol=nq)
  rjm <- matrix(0,nrow=nj, ncol=nq)
  for(g in 1:ng){
    # purrr::map(Estep$Njm,function(x){purrr::invoke(.f=rbind, .x=x)})
    Njm <- Njm + purrr::invoke(rbind,Estep$Njm[[g]])
    rjm <- rjm + purrr::invoke(rbind,Estep$rjm[[g]])
  }

  # Standard Error
  SE <- data.frame(a=rep(0,nj), b=rep(0,nj), c=rep(0,nj))
  for(g in 1:ng){
    # purrr::map(Estep$Njm,function(x){purrr::invoke(.f=rbind, .x=x)})
    Njm <- Njm + purrr::invoke(rbind,Estep$Njm[[g]])
    rjm <- rjm + purrr::invoke(rbind,Estep$rjm[[g]])
  }
  for(j in 1:nj){
    #cat("Optimising item ",j,"\r")
    if(model[j]=="NONE") next
    # Let Nqr tidyr for FS
    N <- Njm[j,]
    rr <- rjm[j,]
    # Fisher score matrix
    FI <- Ijm(N, Xq, t1$a[j], t1$b[j], t1$c[j], D=D, model=model[j])
    # solve
    if(model[j]=="1PL") SE$b[j] <- sqrt(diag(solve(FI)))
    if(model[j]=="2PL") SE[j,c(1,2)] <- sqrt(diag(solve(FI)))
    if(model[j]=="3PL") SE[j,] <- sqrt(diag(solve(FI)))
  }

  # output result
  # Item
  t1 <- data.frame(Item = as.character(Item), t1, model = model)
  # SE
  SE <- data.frame(Item = as.character(Item), SE, model = model)
  # th_dist
  theta_dist <- data.frame(AX)
  colnames(theta_dist) <- as.character(c(1:ng))
  theta_dist$theta <- Xq
  #th_para
  theta_para <- rbind(p_mean,p_sd) %>% data.frame()
  colnames(theta_para) <- as.character(c(1:ng))
  res <- list(para=t1, SE=SE, th_dist=theta_dist, th_para=theta_para, mll_history=mll_history)

  return(res)
}



