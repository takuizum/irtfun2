# ML and MAP estimate function using optim function

LLG <- function(u, theta, phi, a, b, D){
  p <- gptheta(theta = theta, phi = phi, a = a, b = b, D = D)
  sum(u * log(p) + (1-u) * log(p), na.rm = T)
}

fpdG <- function(u, theta, phi, a, b, D){
  p <- gptheta(theta = theta, phi = phi, a = a, b = b, D = D)
  z <- D * a / (sqrt(1 + phi^2*a^2))
  sum(a^2 * phi /(1 + phi^2*a^2) * z * (u - p), na.rm = T)
}

spdG <- function(u, theta, phi, a, b, D){
  p <- gptheta(theta = theta, phi = phi, a = a, b = b, D = D)
  z <- D * a / (sqrt(1 + phi^2*a^2))
  sum(a^2 / (1 + phi^2 * a^2)^2 * z * (u - p), na.rm = T)
}

LLG(dat1[1, -1], fit1$person$theta[1], 1, fit1$item$a, fit1$item$b, D = 1.702)
gr1 <- fpdG(dat1[1, -1], fit1$person$theta[1], 1, fit1$item$a, fit1$item$b, D = 1.702)
He1 <- spdG(dat1[1, -1], fit1$person$theta[1], 1, fit1$item$a, fit1$item$b, D = 1.702)
