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
  C <- phi^(1-v)
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
gr_j_g <- function(r, N, X, Y, t0, D){
  a <- t0[1]
  b <- t0[2]
  p <- gptheta(X, Y, a, b, D)
  ga <- sum(D*(r-N*p)*(X-b)/sqrt((1+Y^2*a^2)^3))
  gb <- sum(-D*(r-N*p)*a/sqrt(1+Y^2*a^2))
  c(ga,gb)
}

Elnk_j_g <- function(r,N,t0,X,Y,D){
  a <- t0[1]
  b <- t0[2]
  p <- gptheta(theta=X,phi=Y,a=a,b=b,D=D)
  A <- r*log(p)
  B <- (N-r)*log(1-p)
  sum(A+B)
}

# gradient
grj_g <- function(r, N, X, Y, a, b, D){
  p <- gptheta(X, Y, a, b, D)
  ga <- sum(D*(r-N*p)*(X-b)/sqrt((1+Y^2*a^2)^3))
  gb <- sum(-D*(r-N*p)*a/sqrt(1+Y^2*a^2))
  c(ga,gb)
}

# Fisher information matrix
Ij_g <- function(r, N, X, Y, a, b, D){
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
#'@author Shin-ichi Mayekawa <mayekawa@@nifty.com>
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
#' @param Gc the grade column.
#' @param bg a mumber of base grade.
#' @param IDc the ID column.
#' @param Ntheta the number of the nodes of theta dist.
#' @param Nphi the number of the nodes of phi dist.
#' @param method the method of optimiser.  Default is "Fisher_Scoring", but \code{\link[stats]{optim}} function also be able to use.
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
#' @param eEM a convergence criterion related to item parameters in EM cycle.
#' @param eMLL a convergence criterion related to negative twice log likelihood in EM cycle.
#' @param eDIST a convergence criterion related to the parameter of theta distribution in EM cycle.
#' @param maxiter_em the number of iteration of EM cycle.
#' @param rm_list a vector of item U want to remove for estimation. NOT list.
#' @param th_dist a type of theta dist."normal" or "empirical"
#' @param print How much information you want to display? from 1 to 3. The larger, more information is displayed.
#' @param esteap logical. If \code{TRUE}, estimate subject theta & phi EAP.
#' @param estdist logical. If \code{TRUE}, estimate population distribution theta & phi.
#'
#' @return the output is a list that has item parameter and person parameter.
#' @examples
#' res <- estGip(x=sim_dat_girt,fc=2, Ntheta=10, Nphi = 5, min_ph = 0.001, max_ph = 2)
#' # check the parameters
#' res$item
#' head(res$person)
#'
#' @export
#'
estGip <- function(x, fc=3, Gc=NULL, bg=1, IDc=1, Ntheta=31, Nphi=10,  method="Fisher_Scoring", th_dist="normal", rm_list=NULL,
                   phi_dist = "invchi", v=3, tau=1, mu_ph=0, sigma_ph=0.25, min_ph=0.01, max_ph=2, paramab=c(1,4),
                   mu_th=0, sigma_th=1, min_th=-4, max_th=4, eEM=0.001, eMLL=1e-6, eDIST=1e-4, maxiter_em=100,
                   print=0, esteap=FALSE, estdist=FALSE){

  if(!(method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent","Fisher_Scoring"))) stop("argument input of `method` is improper string!!")

  # data check
  X <- as.matrix(x[,fc:ncol(x)])
  Item <- colnames(X)
  nj <- ncol(X)
  ni <- nrow(X)
  nq <- Ntheta
  nr <- Nphi
  if(is.null(IDc)) ID <- 1:ni
  else  ID <- x[,IDc]
  if(is.null(Gc)) group <- rep(1,ni)
  else group <- x[,Gc] %>% as.numeric()
  ng <- max(group)
  if(!is.numeric(group)) stop("A group column must be numeric, but has non numeric scalar.")
  if(min(group)==0) stop("A group column must be above 1, but is 0.")
  p_mean <- rep(mu_th, ng)
  p_sd <- rep(sigma_th, ng)
  model <- "GIRT" %>% rep(nj)

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
  AX <- dnorm(Xq, mean=mu_th, sd=sigma_th)/sum(dnorm(Xq, mean=mu_th, sd=sigma_th)) # for theta
  AX <- AX %>% rep.int(times=ng) %>% matrix(ncol=ng)
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
  BY <- BY %>% rep.int(times=ng) %>% matrix(ncol=ng)
  #
  # initial value
  D <- 1.702
  r <- cor(rowSums(X, na.rm = T),X, use = "pair") %>% as.numeric()
  pass <- colMeans(X, na.rm = T)
  a0 <- D*r/sqrt(1-r^2)
  names(a0) <- Item
  b0 <- -log(pass/(1-pass))
  init <- t0 <- t1 <- data.frame(a=a0,b=b0)

  # remove selected item
  if(!is.null(rm_list)){
    rm_ind <- c(1:nj)[Item %in% rm_list]
    model[rm_ind] <- "NONE"
    resp[,rm_ind] <- 0
    ind[,rm_ind] <- 0
    init[rm_ind,] <- t0[rm_ind,] <- t1[rm_ind,] <- 0
    a0 <- init$a
    b0 <- init$b
    cat("Remove Item ", rm_list, " \n")
  }

  # cat
  cat("The number of subject is " ,ni, ".\nThe number of item is ", nj-length(rm_list),
      ".\nThe number of remove item is ", length(rm_list), ".\n")

  # for update imputation
  a1 <- b1 <- numeric(nj)

  # set mll history vector
  mll_history <- c(0)

  # for split unlist vector by group in E step
  g_ind_E <- c(1:ng) %>% matrix() %>% apply(1,rep.int, times=nj*nr*nq) %>% as.data.frame() %>% tidyr::gather(key=dummy, value=value)
  g_ind_E <- g_ind_E$value
  Ngjqr_long <- data.frame(g=g_ind_E, prob=rep(0,ng*nj*nq*nr))
  rgjqr_long <- data.frame(g=g_ind_E, prob=rep(0,ng*nj*nq*nr))
  # empty data.frame
  # group label col
  # item label col
  J <- c(1:nj) %>%
    matrix() %>%
    apply(1,rep.int, times=nr*nq) %>%
    as.data.frame() %>%
    tidyr::gather(key=dummy, value=value)
  J <- J$value# %>% rep.int(times=ng)
  # theta node label col
  Q <- c(1:nq) %>%
    matrix() %>%
    apply(1,rep.int, times=nr) %>%
    as.data.frame() %>%
    tidyr::gather(key=dummy, value=value)
  Q <- Q$value %>%
    rep.int(times=nj)
  # phi label node col
  R <- c(1:nr) %>%
    rep.int(times=nq) %>%
    rep.int(times=nj)
  # data.frame
  Njqr_long <- data.frame(j=J, q=Q, r=R, prob=rep(0,nr*nq*nj))
  rjqr_long <- data.frame(j=J, q=Q, r=R, prob=rep(0,nr*nq*nj))
  #rm(J,Q,R,g_ind_E)

  # convert to longer and longer vector
  X_long <- apply(matrix(Xq), 1, rep.int, nr) %>% as.vector()
  Y_long <- rep(Yr, nq)

  cat("Estimating Item Parameter! \n")

  # stand FLUG
  t <- 0
  convergence <- T
  while(convergence){
    t <- t + 1
    #cat(t,"time EM cycle NOW\n")

    # E step
    Estep <- Estep_girt_mg(X,a0,b0,Xq,AX,Yr,BY,D=1.702,
                           group=group, ind=ind, resp=resp, MLL=mll_history)

    mll <- Estep$MLL[t+1]
    if(print == 0) cat(t ,"times -2 Marginal Loglikelihood is",-2*mll," \r")
    if(print >= 1) cat(t ,"times -2 Marginal Loglikelihood is",-2*mll," \n")
    mll_history <- Estep$MLL

    # unlist
    Ngjqr_long$prob <- Estep$Njqr %>% unlist()
    rgjqr_long$prob <- Estep$rjqr %>% unlist()

    # aggregate
    Njqr_long$prob <- Ngjqr_long %>%
      dplyr::mutate(id=c(1:(nj*nq*nr)) %>% rep(ng)) %>%
      tidyr::spread(key=g, value=prob) %>%
      dplyr::select(-id) %>%
      rowSums()
    rjqr_long$prob <- rgjqr_long %>%
      dplyr::mutate(id=c(1:(nj*nq*nr)) %>% rep(ng)) %>%
      tidyr::spread(key=g, value=prob) %>%
      dplyr::select(-id) %>%
      rowSums()

    # M step
    t1 <- t0
    a0 <- t1$a
    b0 <- t1$b
    for(j in 1:nj){
      if(model[j]=="NONE") next
      Nqr <- Njqr_long$prob[Njqr_long$j==j]
      rqr <- rjqr_long$prob[rjqr_long$j==j]
      # gradient
      if(method != "Fisher_Scoring"){
        res <- optim(par=c(t0[j,1],t0[j,2]), fn=Elnk_j_g, gr=gr_j_g, control = list(fnscale = -1),
                     r=rqr, N=Nqr, X=X_long, Y=Y_long, D=1.702, method = method)
        t1[j,] <- res$par
      }else{
        # Fisher scoring
        gr <- grj_g(rqr, Nqr, X_long, Y_long, t0[j,1], t0[j,2], D=1.702)
        FI <- Ij_g(rqr, Nqr, X_long, Y_long, t0[j,1], t0[j,2], D=1.702)
        # solve
        t1[j,] <- t0[j,] + solve(FI)%*%gr
      }
    }

    # calibration
    # calculate mean and sd
    p_mean_t0 <- p_mean
    p_sd_t0 <- p_sd
    for(g in 1:ng){
      gind <- c(1:nj)[ind[g,]!=0]
      N_of_node <- Ngjqr_long[Ngjqr_long$g==g,] %>%
        dplyr::mutate(J=J, R=R, Q=Q) %>%
        dplyr::filter(J %in% gind)%>%
        dplyr::select(-g) %>%
        tidyr::spread(key=J, value=prob)%>%
        aggregate(by=list(c(1:nq) %>% rep(nr)), FUN=sum) %>%
        dplyr::select(-Group.1, -R, -Q) %>%
        t()
      p_mean[g] <- mean(N_of_node%*%matrix(Xq)/rowSums(N_of_node))
      p_sd[g] <- mean(sqrt(N_of_node%*%matrix((Xq-p_mean[g])^2)/rowSums(N_of_node)))
    }
    # calculate calibration weight
    A <- sigma_th/p_sd[bg]
    K <- mu_th - A*p_mean[bg]
    # calibrate mean and sd
    p_mean <- A*p_mean+K
    p_sd <- A*p_sd
    for(g in 1:ng){
      gind <- c(1:nj)[ind[g,]!=0]
      N_of_node <- Ngjqr_long[Ngjqr_long$g==g,] %>%
        dplyr::mutate(J=J, R=R, Q=Q) %>%
        dplyr::filter(J %in% gind)%>%
        dplyr::select(-g) %>%
        tidyr::spread(key=J, value=prob)%>%
        aggregate(by=list(c(1:nq) %>% rep(nr)), FUN=sum) %>%
        dplyr::select(-Group.1, -R, -Q) %>%
        t()
      if(th_dist=="normal"){
        # Gaussian dist
        AX[,g] <- dnorm(Xq, mean=p_mean[g], sd=p_sd[g])/sum(dnorm(Xq, mean=p_mean[g], sd=p_sd[g]))
      }
      if(th_dist=="empirical"){
        # empirical dist
        constNjq <- rowSums(N_of_node) %>% rep.int(times=nq) %>% matrix(ncol=nq)
        AX[,g] <- (N_of_node/constNjq) %>% colMeans(na.rm = T)
      }
    }
    # calibrate item parameter
    for(j in 1:nj){
      if(model[j]=="NONE") next
      a1[j] <- t1$a[j]/A
      b1[j] <- t1$b[j]*A+K
    }

    t1$a <- a1
    t1$b <- b1

    # change phi dist weight
    # DON'T USE
    #
    # for(g in 1:ng){
    #   gind <- c(1:nj)[ind[g,]!=0]
    #   pN_of_node <- Ngjqr_long[Ngjqr_long$g==g,] %>%
    #     dplyr::mutate(J=J, R=R, Q=Q) %>%
    #     dplyr::filter(J %in% gind)%>%
    #     dplyr::select(-g) %>%
    #     tidyr::spread(key=J, value=prob)%>%
    #     aggregate(by=list(c(1:nr) %>% rep(nq)), FUN=sum) %>%
    #     dplyr::select(-Group.1, -R, -Q) %>%
    #     t()
    #   constpNjr <- rowSums(pN_of_node) %>% rep.int(times=nr) %>% matrix(ncol=nr)
    #   BY[,g] <- (pN_of_node/constpNjr) %>% colMeans(na.rm = T)
    # }

    # convergence check
    conv1 <- c(abs(a0-a1),abs(b0-b1))
    conv4 <- abs(mll_history[t+1]-mll_history[t])
    if(ng > 1){
      conv2 <- c(abs(p_mean - p_mean_t0))
      conv3 <- c(abs(p_sd - p_sd_t0))
    }else{
      conv2 <- conv3 <- 1
    }

    if(all(conv1<eEM)|| conv4<eMLL || all(conv2<eDIST) || all(conv3<eDIST) || maxiter_em == t){
      convergence <- F
      cat("\nConverged!!\n")
      item_para <- as.data.frame(t1)
      break
    }else{
      if(print >= 1){
        cat("Item maximum changed a is", Item[max(abs(a0-a1))==abs(a0-a1)],"=",max(abs(a0-a1))," \n")
        cat("Item maximum changed b is", Item[max(abs(b0-b1))==abs(b0-b1)],"=",max(abs(b0-b1))," \n")
      }
      t0 <- t1
      a0 <- t1$a
      b0 <- t1$b
    }
  }

  # last E step
  Estep <- Estep_girt_mg(X,a1,b1,Xq,AX,Yr,BY,D=1.702,
                         group=group, ind=ind, resp=resp, MLL=mll_history)

  mll <- Estep$MLL[t+1]
  mll_history <- Estep$MLL

  # unlist
  Ngjqr_long$prob <- Estep$Njqr %>% unlist()
  rgjqr_long$prob <- Estep$rjqr %>% unlist()

  # aggregate
  Njqr_long$prob <- Ngjqr_long %>%
    dplyr::mutate(id=c(1:(nj*nq*nr)) %>% rep(ng)) %>%
    tidyr::spread(key=g, value=prob) %>%
    dplyr::select(-id) %>%
    rowSums()
  rjqr_long$prob <- rgjqr_long %>%
    dplyr::mutate(id=c(1:(nj*nq*nr)) %>% rep(ng)) %>%
    tidyr::spread(key=g, value=prob) %>%
    dplyr::select(-id) %>%
    rowSums()

  # Standard Error
  SE <- t1
  for(j in 1:nj){
    #cat("Optimising item ",j,"\r")
    if(model[j]=="NONE") next
    Nqr <- Njqr_long$prob[Njqr_long$j==j]
    rqr <- rjqr_long$prob[rjqr_long$j==j]
    # Fisher score matrix
    FI <- Ij_g(rqr, Nqr, X_long, Y_long, t1[j,1], t1[j,2], D=1.702)
    # solve
    SE[j,] <- sqrt(diag(solve(FI)))
  }

  person_para <- data.frame(ID=ID, GROUP=group, SCORE=rowSums(X,na.rm = T), theta=rep(NA,ni), phi=rep(NA,ni))
  if(esteap==TRUE){
    cat("Estimating EAP! \n")
    # estimate person theta and phi EAP
    for(g in 1:ng){
      gID <- ID[group==g]
      iID <- c(1:ni)[group==g]
      ngi <- iID %>% length() # n of subjects of group g
      hqr <- Estep$hqr[[g]][iID]
      hqr_long <- hqr %>% unlist() %>%
        as.data.frame() %>%
        dplyr::mutate(id = matrix(as.factor(gID)) %>% apply(1, rep.int, times=nq*nr) %>% as.vector()) %>% # add col
        dplyr::mutate(theta = matrix(c(1:nq)) %>% apply(1, rep.int, times=nr) %>% as.vector() %>% rep.int(times=ngi)) %>% # add col
        dplyr::mutate(phi = matrix(c(1:nr)) %>% rep.int(times=ngi*nq) %>% as.vector()) # add col
      hqr_long <- dplyr::rename(.data=hqr_long, prob=.) %>% # rename
        tidyr::spread(key=theta, value=prob) # spread
      for(i in gID){
        i_hqr <- hqr_long %>% dplyr::filter(id==i) %>%
          dplyr::select(-id, -phi) %>%
          as.matrix()
        # theta
        person_para[person_para$ID == i, 4] <- i_hqr %>% colSums() %>% matrix(nrow=1) %*% Xq
        # phi
        person_para[person_para$ID == i, 5] <- i_hqr %>% rowSums() %>% matrix(nrow=1) %*% Yr
      }
    }

  }

  # estimating population distribution
  UX <- (rep(1,nq)/nq) %>% rep.int(times=ng) %>% matrix(ncol=ng)
  UY <- (rep(1,nr)/nr) %>% rep.int(times=ng) %>% matrix(ncol=ng)
  # initialize
  th_mean <- rep(mu_th, ng)
  th_sd <- rep(sigma_th, ng)
  ph_mean <- rep(mu_ph, ng)
  ph_sd <- rep(sigma_ph, ng)
  if(estdist==TRUE){
    cat("Estimating the population theta distribution! \n")
    # theta
    # set mll history vector
    mll_history2 <- c(0)
    t <- 0
    convergence <- T
    while(convergence){
      t <- t + 1
      #cat(t,"time EM cycle NOW\n")

      # E step
      Estep <- Estep_girt_mg(X,t1[,1],t1[,2],Xq,UX,Yr,BY,D=1.702,
                             group=group, ind=ind, resp=resp, MLL=mll_history2)

      mll <- Estep$MLL[t+1]
      if(print == 0) cat(t ,"times -2 Marginal Loglikelihood is",-2*mll," \r")
      if(print >= 1) cat(t ,"times -2 Marginal Loglikelihood is",-2*mll," \n")
      mll_history2 <- Estep$MLL

      # unlist
      Ngjqr_long$prob <- Estep$Njqr %>% unlist()
      rgjqr_long$prob <- Estep$rjqr %>% unlist()

      # aggregate
      Njqr_long$prob <- Ngjqr_long %>%
        dplyr::mutate(id=c(1:(nj*nq*nr)) %>% rep(ng)) %>%
        tidyr::spread(key=g, value=prob) %>%
        dplyr::select(-id) %>%
        rowSums()
      rjqr_long$prob <- rgjqr_long %>%
        dplyr::mutate(id=c(1:(nj*nq*nr)) %>% rep(ng)) %>%
        tidyr::spread(key=g, value=prob) %>%
        dplyr::select(-id) %>%
        rowSums()

      # update the weight of dist
      th_mean_t0 <- th_mean
      th_sd_t0 <- th_sd
      for(g in 1:ng){
        gind <- c(1:nj)[ind[g,]!=0]
        N_of_node <- Ngjqr_long[Ngjqr_long$g==g,] %>%
          dplyr::mutate(J=J, R=R, Q=Q) %>%
          dplyr::filter(J %in% gind)%>%
          dplyr::select(-g) %>%
          tidyr::spread(key=J, value=prob)%>%
          aggregate(by=list(c(1:nq) %>% rep(nr)), FUN=sum) %>%
          dplyr::select(-Group.1, -R, -Q) %>%
          t()
        # empirical dist
        constNjq <- rowSums(N_of_node) %>% rep.int(times=nq) %>% matrix(ncol=nq)
        UX[,g] <- (N_of_node/constNjq) %>% colMeans(na.rm = T)
        # calculate mean and sd
        th_mean[g] <- mean(N_of_node%*%matrix(Xq)/rowSums(N_of_node))
        th_sd[g] <- mean(sqrt(N_of_node%*%matrix((Xq-p_mean[g])^2)/rowSums(N_of_node)))
      }

      # convergence check
      conv2 <- c(abs(th_mean - th_mean_t0))
      conv3 <- c(abs(th_sd - th_sd_t0))

      if(all(conv2<eDIST) || all(conv3<eDIST) || maxiter_em == t){
        convergence <- F
        cat("\nConverged!!\n")
        break
      }
    }

    # phi
    cat("Estimating the population phi distribution!\n")
    # set mll history vector
    mll_history3 <- c(0)
    t <- 0
    convergence <- T
    while(convergence){
      t <- t + 1
      #cat(t,"time EM cycle NOW\n")

      # E step
      Estep <- Estep_girt_mg(X,t1[,1],t1[,2],Xq,AX,Yr,UY,D=1.702,
                             group=group, ind=ind, resp=resp, MLL=mll_history3)

      mll <- Estep$MLL[t+1]
      if(print == 0) cat(t ,"times -2 Marginal Loglikelihood is",-2*mll," \r")
      if(print >= 1) cat(t ,"times -2 Marginal Loglikelihood is",-2*mll," \n")
      mll_history3 <- Estep$MLL

      # unlist
      Ngjqr_long$prob <- Estep$Njqr %>% unlist()
      rgjqr_long$prob <- Estep$rjqr %>% unlist()

      # aggregate
      Njqr_long$prob <- Ngjqr_long %>%
        dplyr::mutate(id=c(1:(nj*nq*nr)) %>% rep(ng)) %>%
        tidyr::spread(key=g, value=prob) %>%
        dplyr::select(-id) %>%
        rowSums()
      rjqr_long$prob <- rgjqr_long %>%
        dplyr::mutate(id=c(1:(nj*nq*nr)) %>% rep(ng)) %>%
        tidyr::spread(key=g, value=prob) %>%
        dplyr::select(-id) %>%
        rowSums()

      # update the weight of dist
      ph_mean_t0 <- ph_mean
      ph_sd_t0 <- ph_sd
      for(g in 1:ng){
        gind <- c(1:nj)[ind[g,]!=0]
        pN_of_node <- Ngjqr_long[Ngjqr_long$g==g,] %>%
          dplyr::mutate(J=J, R=R, Q=Q) %>%
          dplyr::filter(J %in% gind)%>%
          dplyr::select(-g) %>%
          tidyr::spread(key=J, value=prob)%>%
          aggregate(by=list(c(1:nr) %>% rep(nq)), FUN=sum) %>%
          dplyr::select(-Group.1, -R, -Q) %>%
          t()
        constpNjr <- rowSums(pN_of_node) %>% rep.int(times=nr) %>% matrix(ncol=nr)
        UY[,g] <- (pN_of_node/constpNjr) %>% colMeans(na.rm = T)
        # calculate mean and sd
        ph_mean[g] <- mean(pN_of_node%*%matrix(Yr)/rowSums(pN_of_node))
        ph_sd[g] <- mean(sqrt(pN_of_node%*%matrix((Yr-p_mean[g])^2)/rowSums(pN_of_node)))
      }

      # convergence check
      conv2 <- c(abs(ph_mean - ph_mean_t0))
      conv3 <- c(abs(ph_sd - ph_sd_t0))

      if(all(conv2<eDIST) || all(conv3<eDIST) || maxiter_em == t){
        convergence <- F
        cat("\nConverged!! \n")
        break
      }
    }
  }

  # tidyr output data

  theta_dist1 <- data.frame(theta=Xq, AX)
  theta_dist2 <- data.frame(theta=Xq, UX)
  theta_dist1_para <- rbind(p_mean,p_sd) %>% data.frame()
  theta_dist2_para <- rbind(th_mean,th_sd) %>% data.frame()
  phi2_dist1 <- data.frame(phi=Yr, BY)
  phi2_dist2 <- data.frame(phi=Yr, UY)
  phi2_dist2_para <- rbind(ph_mean,ph_sd) %>% data.frame()
  # list
  res <- list(item=item_para, item_SE=SE, person=person_para,
              th_dist1=theta_dist1, phi_dist1=phi2_dist1, th_dist2=theta_dist2, phi_dist2=phi2_dist2,
              th_dist1_para=theta_dist1_para, th_dist2_para=theta_dist2_para, ph_dist1_para=phi2_dist2_para,
              mll_history=mll_history)

  return(res)
}
