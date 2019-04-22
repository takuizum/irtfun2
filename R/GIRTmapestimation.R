# ML and MAP estimate function using optim function
library(tidyverse);library(irtfun2)
#----data gen-----
set.seed(0204)
theta <- rnorm(3000)
phi <- rinvchi(3000, max = 2)
a <- rlnorm(30, sdlog = 0.25)
b <- rnorm(30)
dat <- sim_gen(theta=theta, phi=phi, a=a, b=b)
fit <- estGip(dat, fc = 2, esteap = T)
# for apply
dat2 <- bind_cols(theta = fit$person$theta, phi = fit$person$phi, dat[,-1])
#----prior setting----
tibble(phi = 0.0001:4) %>% ggplot(aes(x = phi)) + stat_function(fun = dinvchi, args = list(v = 2, tau = 3))
tibble(phi = 0.0001:4) %>% ggplot(aes(x = phi)) + stat_function(fun = dlnorm, args = list(meanlog = 0.5, sdlog = 1))
tibble(phi = 0.0001:4) %>% ggplot(aes(x = phi)) + stat_function(fun = dnorm, args = list(mean = 0.5, sd = 1))


#----functions----

# Maximum Likelihood
LLG <- function(u, theta, phi, a, b, D){
  z <- D * a / (sqrt(1 + phi^2*a^2)) * (theta - b)
  p <- 1/(1 + exp(-z))
  sum(u * log(p) + (1-u) * log(1-p), na.rm = T)
}
# Maximau a posteriori(prior is turncated normal distribution)
BLLG <- function(u, theta, phi, a, b, D, mode){
  z <- D * a / (sqrt(1 + phi^2*a^2)) * (theta - b)
  p <- 1/(1 + exp(-z))
  sum(u * log(p) + (1-u) * log(1-p), na.rm = T) + log(dnorm(phi, mean = mode))
}

# for apply version
LLG_apply <- function(dat, a, b, D){
  theta <- dat[1]
  # phi <- dat[2]
  u <- dat[c(-1,-2)]
  opt <- optimise(LLG, interval = c(0.001, 5), u = u, theta = theta, a = a, b = b, D = D, maximum = T)
  opt$maximum
}
# for apply version
BLLG_apply <- function(dat, a, b, D, mode){
  theta <- dat[1]
  # phi <- dat[2]
  u <- dat[c(-1,-2)]
  opt <- optimise(BLLG, interval = c(0.001, 5), u = u, theta = theta, a = a, b = b, D = D, mode = mode, maximum = T)
  opt$maximum
}
# phi penalized
BLLG_apply2 <- function(dat, a, b, D){
  theta <- dat[1]
  mode <- abs(dat[2])
  u <- dat[c(-1,-2)]
  opt <- optimise(BLLG, interval = c(0.001, 5), u = u, theta = theta, a = a, b = b, D = D, mode = mode, maximum = T)
  opt$maximum
}

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

#----optimize function----
optimise(LLG, interval = c(0.001, 5), u = dat[1, -1], theta = fit$person$theta[1], a = fit$item$a, b = fit$item$b, D = 1.702, maximum = T)
phi[1]

# est_phi <- numeric(3000)
# for(i in 1:3000){
#   cat(i, "\r")
#   res <- optimise(LLG, interval = c(0.01, 5), u = as.vector(dat[i, -1]), theta = fit$person$theta[i], a = fit$item$a, b = fit$item$b, D = 1.702, maximum = T)
#   est_phi[i] <- res$maximum
# }

# apply version
phi_ml <- apply(dat2, 1, LLG_apply, a = fit$item$a, b = fit$item$b, D = 1.702)
phi_map <- apply(dat2, 1, BLLG_apply, a = fit$item$a, b = fit$item$b, D = 1.702, mode = 1)

cor(phi, phi_ml)
cor(phi, phi_map)
plot(phi_ml, phi_map)
plot(phi, phi_map)
plot(phi, phi_ml)


# Newton-Raphton
t0 <- 1
s <- 0 # counter
while(TRUE){
  s <- s + 1
  cat(s, "\n")
  gr1 <- fpdG(dat[3, -1], fit$person$theta[3], t0, fit$item$a, fit$item$b, D = 1.702)
  # Hi1 <- fiG(dat[3, -1], fit$person$theta[3], t0, fit$item$a, fit$item$b, D = 1.702)
  Fi1 <- fiG(dat[3, -1], fit$person$theta[3], t0, fit$item$a, fit$item$b, D = 1.702)
  # t1 <- t0 - gr1/Hi1
  t1 <- t0 + gr1/Fi1
  if(abs(t1-t0) < 0.0001) break
  t0 <- t1
}

#----Vizualization a maximum point of objective function----
phi_vec <- seq(0.001, 4, length.out = 101)
# ML function
# plot(y = LLG(dat[1, -1], fit$person$theta[1], phi_vec[1], fit$item$a, fit$item$b, D = 1.702), x = phi_vec[1], xlim = c(0, 4), ylim = c(-20, 0), type = "o", ylab = "", xlab = "")
# for(j in 1:100){
#   par(new = T)
#   plot(y = LLG(dat[1, -1], fit$person$theta[1], phi_vec[j+1], fit$item$a, fit$item$b, D = 1.702), x = phi_vec[j+1], xlim = c(0, 4), ylim = c(-20, 0), type = "o", ylab = "", xlab = "")
# }

subject <- 10
# ML
tibble(phi = 0.001:4) %>% ggplot(aes(x = phi)) +
  stat_function(fun = LLG_visual, args = list(u = dat[subject, -1], a = fit$item$a, b = fit$item$b, theta = fit$person$theta[subject], D = 1.702))

# log posterior distribution function
# plot(y = BLLG(dat[1, -1], fit$person$theta[1], phi_vec[1], fit$item$a, fit$item$b, D = 1.702, mode = 1), x = phi_vec[1], xlim = c(0, 4), ylim = c(-20, -15), type = "o", ylab = "", xlab = "")
# for(j in 1:100){
#   par(new = T)
#   plot(y = BLLG(dat[1, -1], fit$person$theta[1], phi_vec[j+1], fit$item$a, fit$item$b, D = 1.702, mode = 1), x = phi_vec[j+1], xlim = c(0, 4), ylim = c(-20, -15), type = "o", ylab = "", xlab = "")
# }
# MAP
tibble(phi = 0.001:4) %>% ggplot(aes(x = phi)) +
  stat_function(fun = BLLG_visual,
                args = list(u = dat[subject, -1], a = fit$item$a, b = fit$item$b, theta = fit$person$theta[subject], D = 1.702, mode = 1))


#----person fit index----
# E_3(theta) : expectation
E3 <- function(theta, a, b, D){
  z <- D * a * (theta - b)
  p <- 1/(1 + exp(-z))
  sum(p * log(p) + (1-p) * log(1-p),na.rm = T)
}

# sigma3
S3 <- function(theta, a, b, D){
  z <- D * a * (theta - b)
  p <- 1/(1 + exp(-z))
  sqrt(sum(p * (1 - p) * (log(p/(1-p)))^2,na.rm = T))
}

# person fit index: z_3 statistics
pfit <- function(dat,a,b,D){
  #decompose dat
  theta <- dat[1]
  xi <- dat[-1]
  (LLG(xi, theta, phi = 0,  a, b, D) - E3(theta, a, b, D)) / S3(theta, a, b, D)
}

fitindex <- apply(dat2[,-2], 1, pfit, a = fit$item$a, b = fit$item$b, D = 1.702)
# MAP
phi_map <- apply(cbind(dat2$theta, fitindex, dat2[,c(-1,-2)]), 1, BLLG_apply2, a = fit$item$a, b = fit$item$b, D = 1.702)
phi[1]
plot(phi_map, phi)
cor(phi, phi_map)
plot(phi_map, fitindex)
plot(phi, fitindex)

# ECI 関数

