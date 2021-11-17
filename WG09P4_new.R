# Group_09, Chang LU, Hongxuan ZHU, Zheyi SHEN
# https://github.com/Zoey1102/BFGS_Optimizer

## Practical 4 Group 9 ####
## This practical is aimed to write a R function, bfgs, implementing the BFGS quasi-Newton minimization method.

fdg <- function(theta,f,...) {
  if(is.null(attr(f(theta,...), 'gradient'))){
    p <- length(theta); g <- rep(0,p) ## empty vector to store the gradient value
    eps <- 1e-7
    for (i in 1:p) {
      th_add <- th_sub <- theta
      th_add[i] <- th_add[i] + eps; th_sub[i] <- th_sub[i] - eps
      g[i] <- (f(th_add,...)-f(th_sub,...))/(2*eps)
    }
  }
  else{
    g = attr(f(theta,...),'gradient')
  }
  g
}

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100) {
  " theta: vector of initial optimization parameters
        f: objective function to minimize
      tol: convergence tolerance
   fscale: estimate of f at the optimum
    maxit: maximum number of iterations 
   return: a list of value at th minimum, including
           - objective function, parameters, iterations, gradients, approximate Hessian matrix"
  
  g0 <- fdg(theta,f,...)
  f0 <- f(theta,...)
  if (abs(f0) == Inf | sum(abs(g0)) == Inf) stop('the initiail input is infinite')
  
  p <- length(theta)
  B <- diag(p) ## initialize approximate Hessian matrix as a diagonal matrix
  c1 <- 0.01; c2 <- 0.9; maxit.s <- 20
  
  for (i in 1:maxit) {
    g0 <- fdg(theta,f,...)
    f0 <- f(theta,...) ## return the function value
    if (max(abs(g0))<(abs(f0)+fscale)*tol) break
    s <- -B %*% g0
    for (k in 1:maxit.s) {
      theta1 <- theta+s
      g1 <- fdg(theta1,f,...)
      f1 <- f(theta1,...)
      if (f1<=f0+c1*g0%*%s & g1%*%s>=c2*g0%*%s) {
        y <- g1-g0
        ro <- 1/as.numeric((t(s)%*%y))
        B <- (diag(p)-ro*s%*%t(y))%*%B%*%(diag(p)-ro*y%*%t(s))+ro*s%*%t(s)
        theta <- theta1
        H <- 0.5*(t(B)+B)
        break
      } else s <- s/2
    }
    if (k==maxit.s) {
      warning("function not a normal convergence")
      break
    }
  }
  if (i==maxit) warning("iteration limit reached")
  list(f=f1,theta=theta,iter=i,g=g1,H=H)
}
