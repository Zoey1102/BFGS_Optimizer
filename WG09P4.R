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
  
  fn <- function(theta) f(theta,...) ## return the function value
  m = 1 ## number of iterations
  hep = 1e-07 ## a tiny value for calculating gradient
  gr = rep(0,length(theta)) ## empty vector to store the gradient value
  bk = diag(length(theta)) ## initialize approximate Hessian matrix as a diagonal matrix
  while(m <= maxit){
    # calculate gradient value
    for(i in 1:length(theta)){
      pos = rep(0,length(theta))
      pos[i] = 1
      theta_add = theta + pos*hep
      theta_substract = theta - pos*hep
      gr[i] = (fn(theta_add) - fn(theta_substract))/(2 * hep)
    }
    # figure out the step size 'a'
    dk = -gr %*% solve(bk)
    a = 10
    lameda = 0.8
    c = 0.5
    b = TRUE
    while(b){
      if(fn(theta + a*dk)<fn(theta) + a * c * gr %*% t(dk)) b = FALSE
      else a = a * lameda
    }
    # renew theta
    theta_new = theta + a * dk
    # renew the approximate Hessian matrix
    sk = theta_new - theta
    gr_new = rep(0,length(theta))
    for(i in 1:length(theta)){
      pos = rep(0,length(theta))
      pos[i] = 1
      theta_add = theta_new + pos*hep
      theta_substract = theta_new - pos*hep
      gr_new[i] = (fn(theta_add) - fn(theta_substract))/(2 * hep)
    }
    yk = gr_new - gr
    # if sum(abs(yk)) < tol, which means this function converge, then loop break
    if(sum(abs(yk)) < tol) break
    # renew bk
    else{
      if(yk %*% t(sk) > 0){
        bk = bk - (bk %*% t(sk) %*% sk %*% bk)/as.double(sk %*% bk %*% t(sk)) +
          (t(t(yk)) %*% yk)/as.double(yk %*% t(sk))
      }
      else bk = bk
    }
    # renew theta,m
    theta = theta_new
    m = m +1
  }
  gr = gr_new
  f = fn(theta)
  list(f = f,theta = theta,iter = m,g = gr,h = bk)
}
