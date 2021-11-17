# Group_09, Chang LU, Hongxuan ZHU, Zheyi SHEN
# https://github.com/Zoey1102/BFGS_Optimizer

## Practical 4 Group 9 ####
## This practical is aimed to write a R function, bfgs, implementing the BFGS quasi-Newton minimization method.

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100) {
  " theta: vector of initial optimization parameters
        f: objective function to minimize
      tol: convergence tolerance
   fscale: estimate of f at the optimum
    maxit: maximum number of iterations 
   return: a list of value at th minimum, including
           - objective function, parameters, iterations, gradients, approximate Hessian matrix"
  f <- f(theta,...) ## return the function value
  p <- length(theta)
  B <- diag(p) ## initialize approximate Hessian matrix as a diagonal matrix
  c2 <- 0.9; maxit.s <- 20
  
  fdg <- function(theta,f,...) {
    g <- rep(0,p) ## empty vector to store the gradient value
    eps <- 1e-7
    for (i in 1:p) {
      pos = rep(0,length(theta)); pos[i] = 1
      theta_add <- theta + pos*eps
      theta_substract <- theta - pos*eps
      g[i] <- (f(theta_add)-f(theta_substract))/(2*esp)
    }
  }
  
  for (i in 1:maxit) {
    g0 <- fdg(theta,f,...)
    if (max(abs(g0))<(abs(f)+fscale)*tol) break
    s <- -bk %*% g0
    for (k in 1:maxit.s) {
      if (t(fdg(theta+s,f,...))*s>=c2*t(g0)*s) {
        theta <- theta+s
        f <- f(theta+s)
        g <- fdg(theta+s)
        y <- g-g0
        ro <- 1/(t(s)*y)
        B <- (I(p)-ro*s%*%t(y))%*%B%*%(I(p)-ro*s%*%t(y))+ro*s*t(s)
      } else s <- s/2
    }
    if (k==max.half) {
      warning("function not a normal convergence")
      break
    }
  }
  if (i==maxit) warning("iteration limit reached")
  list(f=f,theta=theta,iter=i,g=g,H=H)
}