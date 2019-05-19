# ML and MAP estimate function using optim function


#----functions----
# theta and phi joint maximum likelihood (optimise) ----
# optim関数を使って，二変数最適化をする。
LLgirt <- function(para, u, a, b, D){
  theta <- para[1]
  phi <- para[2]
  p <- gptheta(theta = theta, phi = phi, a = a, b = b, D = D)
  sum(log(p)*u + log(1-p)*(1-u), na.rm = T)
}

LLBgirt <- function(para, u, a, b, D, mu_th = 0, sigma_th = 1, mu_ph = 1, sigma_ph = 1){
  theta <- para[1]
  phi <- para[2]
  p <- gptheta(theta = theta, phi = phi, a = a, b = b, D = D)
  sum(log(p)*u + log(1-p)*(1-u), na.rm = T) + log(dnorm(theta, mu_th, sigma_th)) + log(dnorm(phi, mu_ph, sigma_ph))
}

LLgrit_gr <- function(para, u, a, b, D){
  theta <- para[1]
  phi <- para[2]
  p <- gptheta(theta = theta, phi = phi, a = a, b = b, D = D)
  A <- D*a/(sqrt(1+phi^2*a^2))
  B <- D*a^3*phi*(theta-b) / (1+phi^2*a^2)^(3/2)
  thetagr <- sum(A*(u-p), na.rm = T)
  phigr <- sum(B*(u-p), na.rm = T)
  c(thetagr, phigr)
}

LLBgrit_gr <- function(para, u, a, b, D, mu_th = 0, sigma_th = 1, mu_ph = 1, sigma_ph = 1){
  theta <- para[1]
  phi <- para[2]
  p <- gptheta(theta = theta, phi = phi, a = a, b = b, D = D)
  A <- D*a/(sqrt(1+phi^2*a^2))
  B <- D*a^3*phi*(theta-b) / (1+phi^2*a^2)^(3/2)
  thetagr <- sum(A*(u-p), na.rm = T) - 1/sigma_th^2 * (theta - mu_th)
  phigr <- sum(B*(u-p), na.rm = T) - 1/sigma_ph^2 * (phi - mu_ph)
  c(thetagr, phigr)
}

# optim(par = c(0, 1), fn = LLgirt, u = c(1,0,1,1), a = c(1,1,1,1), b = c(0,0,0,0), D = 1.0, method = "BFGS", hessian = T, control = list(fnscale = -1))
# optim(par = c(0, 1), fn = LLgirt, gr = LLgrit_gr, u = c(0,0,1,1), a = c(1,1,1,1), b = c(-1,0,1,2), D = 1.0, method = "BFGS", hessian = T, control = list(fnscale = -1))
# optim(par = c(0, 1), fn = LLgirt, gr = LLgrit_gr, u = c(1,1,1,0,1), a = c(1,1,1,1,1), b = c(-2,-1,0,1,2), D = 1.0, method = "BFGS", hessian = T, control = list(fnscale = -1))
# optim(par = c(0, 1), fn = LLBgirt, gr = LLBgrit_gr, u = c(1,0,1,0,0), a = c(1,1,1,1,1), b = c(-2,-1,0,1,2), D = 1.0, method = "BFGS", hessian = T, control = list(fnscale = -1))
# optim(par = c(0, 1), fn = LLgirt, u = c(1,0,1,1,1), a = c(1,1,1,1,1), b = c(-2,-1,0,1,2), D = 1.0, hessian = T, control = list(fnscale = -1))



# first partial derivative
fpdG <- function(u, theta, phi, a, b, D){
  z <- D * a / (sqrt(1 + phi^2*a^2)) * (theta - b)
  p <- 1/(1 + exp(-z))
  sum(a^2 * phi /(1 + phi^2*a^2) * z * (u - p), na.rm = T)
}
# second partial derivative
spdG <- function(u, theta, phi, a, b, D){
  z <- D * a / (sqrt(1 + phi^2*a^2)) * (theta - b)
  p <- 1/(1 + exp(-z))
  sum(a^2 / (1 + phi^2 * a^2)^2 * z * (u - p), na.rm = T)
}
# Fisher information
fiG <- function(u, theta, phi, a, b, D){
  z <- D * a / (sqrt(1 + phi^2*a^2)) * (theta - b)
  p <- 1/(1 + exp(-z))
  D^2 * sum(a^6*phi^2*(theta - b)^2 / (1 + phi^2 * a^2)^3 * p * (1 - p), na.rm = T)
}

# visualization
LLG_visual <- function(u, a, b, theta, phi, D){
  apply(matrix(phi), 1, LLG, u = u, theta = theta, a = a, b = b, D = D)
}
BLLG_visual <- function(u, a, b, theta, phi, D, mode){
  apply(matrix(phi), 1, BLLG, u = u, theta = theta, a = a, b = b, D = D, mode)
}
# test
LLG(dat[1, -1], fit$person$theta[1], 1, fit$item$a, fit$item$b, D = 1.702)
gr1 <- fpdG(dat[1, -1], fit$person$theta[1], 1, fit$item$a, fit$item$b, D = 1.702)
He1 <- spdG(dat[1, -1], fit$person$theta[1], 1, fit$item$a, fit$item$b, D = 1.702)
Fi1 <- fiG(dat[1, -1], fit$person$theta[1], 1, fit$item$a, fit$item$b, D = 1.702)

1 - gr1/He1
1 + gr1/Fi1

LLG_visual(dat[1, -1], fit$item$a, fit$item$b, fit$person$theta[1], c(0.1, 1, 2, 3), D = 1.702)


# ECI 関数----
ECI <- function(dat, a, b, D = 1.702){
  theta <- dat[1]
  u <- dat[-1]
  # T(¥theta)
  z <- D*a*(theta - b)
  t <- mean(1/(1+exp(-z)))
  p <- 1/(1+exp(-z))
  numerator <- sum((p-u)*(p-t))
  denominator <- sqrt(sum(p*(1-p)*(p-t)^2))
  # numerator/denominator
  sum(abs((p-u)*(p-t)), na.rm = T)
}

# eci1 <- apply(dat2, 1, ECI, a = fit$item$a, b = fit$item$b)
#
# hist(fit$person$phi)
# hist(fitindex)
# hist(eci1)
# plot(phi, eci1)
# phi_map <- apply(cbind(dat2$theta, eci1, dat2[,c(-1,-2)]), 1, BLLG_apply2, a = fit$item$a, b = fit$item$b, D = 1.702, mode = 1, sigma = 1)
# plot(eci1, phi_map)
# cor(phi_map, phi)
# plot(fit$person$phi, phi_map)

# theta and phi optimise function-----
mapG_apply <- function(xi, a, b, D, method, mu_th, sigma_th, mu_ph, sigma_ph, Hessian = F){
  # initial value
  mm <- sum(!is.na(xi))
  if(sum(xi, na.rm = TRUE) == 0){
    t0 <- log(0.5)
  }else if(sum((xi==1)*1,na.rm=T) == length(na.omit(xi))){
    t0 <- log(mm-0.5)
  }else{
    t0 <- log(sum(xi, na.rm = TRUE)/(mm-sum(xi, na.rm = TRUE)))
  }
  eci <- ECI(c(t0, xi), a = a, b = b, D = D)
  # cat(eci ," | ")
  # optimise
  opt <- optim(c(t0,eci), fn = LLBgirt, gr = LLBgrit_gr, method = method, u = xi, a = a, b = b, D = D,
               mu_th = mu_th, sigma_th = sigma_th, mu_ph = mu_ph, sigma_ph = sigma_ph, hessian = Hessian)
  if(Hessian){
    res <- diag(opt$hessian)
  } else {
    res <- opt$par
  }
  res
}

mlG_apply <- function(xi, a, b, D, method, Hessian = F){
  # initial value
  mm <- sum(!is.na(xi))
  if(sum(xi, na.rm = TRUE) == 0){
    t0 <- log(0.5)
  }else if(sum((xi==1)*1,na.rm=T) == length(na.omit(xi))){
    t0 <- log(mm-0.5)
  }else{
    t0 <- log(sum(xi, na.rm = TRUE)/(mm-sum(xi, na.rm = TRUE)))
  }
  eci <- ECI(c(t0, xi), a = a, b = b, D = D)
  # optimise
  opt <- optim(c(t0,eci), fn = LLgirt, gr = LLgrit_gr, method = method, u = xi, a = a, b = b, D = D, hessian = Hessian)
  if(Hessian){
    res <- diag(opt$hessian)
  } else {
    res <- opt$par
  }
  res
}

estGtheta <- function(xall, param, IDc = 1, fc = 2, gc = 2, est = "MAP", method = "BFGS",
                      mu_th = 0, sigma_th = 1, mu_ph = 1, sigma_ph = 1, D = 1.702){
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
  if(is.data.frame(param) && param %>% colnames() %in% c("a", "b") %>% sum() == 2){
    a <- param$a
    x.all <- x.all[, a != 0]

    # set item parameter data
    param <- param[param$a != 0, ]
    a <- param$a
    b <- param$b
  } else {
    a <- param[,1]
    x.all <- x.all[, a != 0]

    # set item parameter data
    param <- param[param[,1] != 0, ]
    a <- param[,1]
    b <- param[,2]
  }

  # Number of Subjects"
  n <- nrow(x.all)
  ng <- numeric(G)
  # number of items
  m <-length(a)
  xscore <- rowSums(x.all,na.rm=T)
  cat("n of subject is ",n,".\n")
  cat("n of item is ",m,".\n")
  cat(est," ESTIMATION IS RUNNING NOW!\n")

  groupitem <- numeric(G)
  for(g in 1:G){
    key <- group==g
    temp <- x.all[key,]
    temp <- temp %>% colSums(na.rm = T)
    groupitem[g] <- temp[temp!=0] %>% length()
  }

  if(est == "MAP"){
    cat("MAP estimation is buggynow. \n")
    t(apply(x.all, 1, mapG_apply, a = a, b = b, method = method, mu_th = mu_th, sigma_th = sigma_th, mu_ph = mu_ph, sigma_ph = sigma_ph, Hessian = F, D = D))
  } else if(est == "MLE"){
    t(apply(x.all, 1, mlG_apply, a = a, b = b, method = method, Hessian = F, D = D))
  }
}


map2 <- estGtheta(dat, param = fit$item, fc = 2, gc = 0, est = "MAP", method = "BFGS")

# for graph
LLgirt2 <- function(theta, phi, u, a, b, D){
  p <- gptheta(theta = theta, phi = phi, a = a, b = b, D = D)
  sum(log(p)*u + log(1-p)*(1-u))
}
LLgirt_apply <- function(theta, phi, u, a, b, D){
  mapply(FUN = LLgirt2, theta, phi, MoreArgs = list(u = u, a = a, b = b, D = D), SIMPLIFY = T) %>% as.vector
}
tibble(theta = seq(-4, 4, length.out = 31) %>% rep(31), phi = apply(matrix(seq(0, 4, length.out = 31)), 1, rep, 31) %>% as.vector) %>%
  mutate(LL = LLgirt_apply(theta = theta, phi = phi, u = dat[3,-1], a = fit$item$a, b = fit$item$b, D = 1.702)) %>%
  ggplot(aes(x = theta, y = phi, z = LL)) + geom_contour(bins = 200)

test %>% ggplot(aes(x = theta, y = phi, z = LL)) + geom_contour(bins = 200)
test %>% ggplot(aes(x = theta, y = phi, z = LL)) + geom_raster(aes(fill = LL), hjust = 0, vjust = 0, interpolate = F)
