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
  
  fn <- function(theta) f(theta,...)
  if(abs(fn(theta)) == Inf) warning('the initiail input is infinite')
  m = 1 ## number of iterations
  bk = diag(length(theta)) ## initialize approximate Hessian matrix as a diagonal matrix
  if(is.null(attr(f(theta,...), 'gradient'))){
    hep = 1e-07 ## a tiny value for calculating gradient
    gr = rep(0,length(theta)) ## empty vector to store the gradient value
    while(m <= maxit){
      # calculate gradient value
      for(i in 1:length(theta)){
        pos = rep(0,length(theta))
        pos[i] = 1
        theta_add = theta + pos*hep
        theta_substract = theta - pos*hep
        gr[i] = (fn(theta_add) - fn(theta_substract))/(2 * hep)
      }
      if(max(abs(gr)) < (abs(fn(theta))+fscale)*tol) break
      else{
        # figure out the step size 'a'
        dk = -gr %*% solve(bk)
        a = 10
        lameda = 0.8
        c1 = 0.5
        c2 = 0.9
        b = TRUE
        while(b){
          theta_dk = theta + a * dk
          gr_dk = rep(0,length(theta))
          for(i in 1:length(theta)){
            pos = rep(0,length(theta))
            pos[i] = 1
            theta_add = theta_dk + pos*hep
            theta_substract = theta_dk - pos*hep
            gr_dk[i] = (fn(theta_add) - fn(theta_substract))/(2 * hep)
          }
          if(fn(theta + a*dk) <= fn(theta) + a * c1 * gr %*% t(dk)
             & gr_dk %*% t(dk) >= c2 * gr %*% t(dk)) b = FALSE
          else a = a * lameda
        }
        # renew theta
        theta_new = theta + a * dk
        if(fn(theta) - fn(theta_new) < 0) warning('The optimization process cannot reduce the value of the function')
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
        # renew bk
        if(yk %*% t(sk) > 0){
          bk = bk - (bk %*% t(sk) %*% sk %*% bk)/as.double(sk %*% bk %*% t(sk)) +
            (t(t(yk)) %*% yk)/as.double(yk %*% t(sk))
        }
        else bk = bk
      }
      m = m + 1
      theta = theta_new
    }
  }
  # input 'f' have attributes 'gradients'
  else{
    while(m < maxit){
      gr = attr(f(theta,...),'gradient')
      if(max(abs(gr)) < (abs(fn(theta))+fscale)*tol) break
      else{
        # figure out the step size 'a'
        dk = -gr %*% solve(bk)
        a = 10
        lameda = 0.8
        c1 = 0.5
        c2 = 0.9
        b = TRUE
        while(b){
          theta_dk = theta + a * dk
          gr_dk = attr(f(theta_dk,...),'gradient')
          if(fn(theta + a*dk) < fn(theta) + a * c1 * gr %*% t(dk)
             & gr_dk %*% t(dk) >= c2 * gr %*% t(dk)) b = FALSE
          else a = a * lameda
        }
        # renew theta
        theta_new = theta + a * dk
        if(fn(theta) - fn(theta_new) < 0) warning('The optimization process cannot reduce the value of the function')
        # renew the approximate Hessian matrix
        sk = theta_new - theta
        gr_new = attr(f(theta_new,...),'gradient')
        yk = gr_new - gr
        # if sum(abs(yk)) < tol, which means this function converge, then loop break
        # renew bk
        if(yk %*% t(sk) > 0){
          bk = bk - (bk %*% t(sk) %*% sk %*% bk)/as.double(sk %*% bk %*% t(sk)) +
            (t(t(yk)) %*% yk)/as.double(yk %*% t(sk))
        }
        else bk = bk
      }
      m = m + 1
      theta = theta_new
    }
  }
  gr = gr_new
  if(m >= maxit & max(abs(gr)) > (abs(fn(theta))+fscale)*tol) warning('maxit is reached without convergence')
  fn = fn(theta)
  list(f = fn,theta = theta,iter = m,g = gr,h = bk)
}
