# IRT true score equating function
#ptheta <- function(theta,a,b,c,D){
#  c + (1-c)/(1+exp(-D*a*(theta-b)))
#}

func1 <- function(tau,theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  tau - sum(p,na.rm = T)
}

func2 <- function(theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  -sum(D*a*(1-p)*(p-c)/(1-c), na.rm=T)
}

ty <- function(theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  sum(p,na.rm=T)
}

#'Newton-Raphton method in IRT true score equating
#'Reference: Kolen & Brennan (2014, p.194)
#'@param tau true score, integer, vector.
#'@param t0 initial value of theta in NR.
#'@param a slope parameter vector.
#'@param b location parameter vector.
#'@param c asymptote parameter vector.
#'@param D factor constant.
#'@export

tse <- function(tau,t0,a,b,c,D){

  # Newton-raphton
  conv <- 0
  iter <- 0
  while(conv == 0){
    iter <- iter + 1
    t1 <- t0 - func1(tau,t0,a,b,c,D)/func2(t0,a,b,c,D)
    if(abs(t1-t0) < 0.001) conv <- 1
    t0 <- t1
  }
  #res <- data.frame(tau=tau,theta=t0,iter=iter)
  res <- c(tau,t0,iter)
  return(res)
}

# Hypothetical data
#a <- c(0.6,1.2,1.0,1.4,1.0)
#b <- c(-1.7,-1.0,0.8,1.3,1.4)
#c <- c(0.2,0.2,0.25,0.25,0.2)
#sum(c)
#tau <- 2
#t0 <- -2
#D <- 1.7
#tse(tau,t0=-2,a,b,c,D=1.7)

#tau <- c(2,3,4) %>% as.matrix() # sum(c)<tau<length(a)
#apply(tau,1,tse,t0=0,a=a,b=b,c=c,D=1.702) %>% t()

#'IRT true score equating
#'Reference: Kolen & Brennan (2014, p.194)
#'@param paraX df of itemparameter on Form X.
#'@param paraY df of itemparameter on Form Y.
#'@param D factor constant.
#'@export
irt_ytx <- function(paraX,paraY,D=1.702){

  # paraX
  keyx <- paraX$a != 0
  ax <- paraX$a[keyx]
  bx <- paraX$b[keyx]
  cx <- paraX$c[keyx]

  # paraY
  keyy <- paraY$a != 0
  ay <- paraY$a[keyy]
  by <- paraY$b[keyy]
  cy <- paraY$c[keyy]

  tau_x <- floor(sum(cx,na.rm = T)+1):(length(ax)-1) %>% matrix()
  theta <- apply(tau_x,1,tse,a=ax,b=bx,c=cx,D=D,t0=0) %>% t()

  tau_y <- apply(theta[,2] %>% as.matrix,1,ty,a=ay,b=by,c=cy,D=D)

  tau_x <- c(0,tau_x,length(ax))
  tau_y <- c(0,tau_y,length(ay))
  theta <- c(NA,theta[,2],NA)

  res <- data.frame(theta = theta,tau_X=tau_x, tau_Y=tau_y)
  res
}



