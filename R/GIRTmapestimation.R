# ML and MAP estimate function using optim function

LLG <- function(u, theta, phi, a, b, D){
  z <- D * a / (sqrt(1 + phi^2*a^2)) * (theta - b)
  p <- 1/(1 + exp(-z))
  sum(u * log(p) + (1-u) * log(1-p), na.rm = T)# + log(dlnorm(phi))
}

fpdG <- function(u, theta, phi, a, b, D){
  z <- D * a / (sqrt(1 + phi^2*a^2)) * (theta - b)
  p <- 1/(1 + exp(-z))
  sum(a^2 * phi /(1 + phi^2*a^2) * z * (u - p), na.rm = T)
}

spdG <- function(u, theta, phi, a, b, D){
  z <- D * a / (sqrt(1 + phi^2*a^2)) * (theta - b)
  p <- 1/(1 + exp(-z))
  sum(a^2 / (1 + phi^2 * a^2)^2 * z * (u - p), na.rm = T)
}

fiG <- function(u, theta, phi, a, b, D){
  z <- D * a / (sqrt(1 + phi^2*a^2)) * (theta - b)
  p <- 1/(1 + exp(-z))
  D^2 * sum(a^6*phi^2*(theta - b)^2 / (1 + phi^2 * a^2)^3 * p * (1 - p), na.rm = T)
}

LLG(dat[1, -1], fit$person$theta[1], 1, fit$item$a, fit$item$b, D = 1.702)
gr1 <- fpdG(dat[1, -1], fit$person$theta[1], 1, fit$item$a, fit$item$b, D = 1.702)
He1 <- spdG(dat[1, -1], fit$person$theta[1], 1, fit$item$a, fit$item$b, D = 1.702)
Fi1 <- fiG(dat[1, -1], fit$person$theta[1], 1, fit$item$a, fit$item$b, D = 1.702)

1 - gr1/He1
1 + gr1/Fi1

# optimize function
optimise(LLG, interval = c(0.001, 5), u = dat[1, -1], theta = fit$person$theta[1], a = fit$item$a, b = fit$item$b, D = 1.702, maximum = T)
phi[1]

est_phi <- numeric(3000)
for(i in 1:3000){
  cat(i, "\r")
  res <- optimise(LLG, interval = c(0.01, 5), u = dat[i, -1], theta = fit$person$theta[i], a = fit$item$a, b = fit$item$b, D = 1.702, maximum = T)
  est_phi[i] <- res$maximum
}

cor(phi, est_phi)
plot(phi, est_phi)


phi_vec <- seq(0.001, 4, length.out = 101)
plot(y = LLG(dat[1, -1], fit$person$theta[1], phi_vec[1], fit$item$a, fit$item$b, D = 1.702), x = phi_vec[1], xlim = c(0, 4), ylim = c(-20, -18), type = "o", ylab = "", xlab = "")
for(j in 1:100){
  par(new = T)
  plot(y = LLG(dat[1, -1], fit$person$theta[1], phi_vec[j+1], fit$item$a, fit$item$b, D = 1.702), x = phi_vec[j+1], xlim = c(0, 4), ylim = c(-20, -18), type = "o", ylab = "", xlab = "")
}

# Newton-Raphton
t0 <- 0.5
s <- 0 # counter
while(TRUE){
  s <- s + 1
  cat(s, "\r")
  gr1 <- fpdG(dat[1, -1], fit$person$theta[1], t0, fit$item$a, fit$item$b, D = 1.702)
  Fi1 <- spdG(dat[1, -1], fit$person$theta[1], t0, fit$item$a, fit$item$b, D = 1.702)
  t1 <- t0 + gr1/Fi1
  if(abs(t1-t0) < 0.0001) break
  t0 <- t1
}
