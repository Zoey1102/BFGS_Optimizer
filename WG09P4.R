# Group_09, Chang LU, Hongxuan ZHU, Zheyi SHEN
# https://github.com/Zoey1102/BFGS_Optimizer

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
  # accept the arguments and return the function value
  fn <- function(theta) f(theta,...)
  # initiate the number of iterations
  m = 1
  # set a tiny value 'hep' for figuring gradient
  hep = 1e-07
  # set an empty vector 'gr' to store the gradient value
  g = rep(0,length(theta))
  # set a diagonal matrix as the initial approximate Hessian matrix
  bk = diag(length(theta))
  while(m <= maxit){
    # calculate gradient value
    for(i in 1:length(theta)){
      pos = rep(0,length(theta))
      pos[i] = 1
      theta_add = theta + pos*hep
      theta_substract = theta - pos*hep
      g[i] = (fn(theta_add) - fn(theta_substract))/(2 * hep)
    }
    # when the gradient value meet the request, loop end
    if(sum(abs(g)) < tol) break
    # when fail to meet the request, figure out the step size 'a'
    else{
      # figure out update direction 'dk'
      dk = -g %*% solve(bk)
      a = 10
      lameda = 0.8
      c = 0.5
      b = TRUE
      while(b){
        if(fn(theta + a*dk)<fn(theta) + a * c * g %*% t(dk)) b = FALSE
        else a = a * lameda
      }
      # renew theta
      theta_new = theta + a * dk
      # renew the approximate Hessian matrix
      sk = theta_new - theta_new
      g_new = rep(0,length(theta))
      for(i in 1:length(theta)){
        pos = rep(0,length(theta))
        pos[i] = 1
        theta_add = theta_new + pos*hep
        theta_substract = theta_new - pos*hep
        g_new[i] = (fn(theta_add) - fn(theta_substract))/(2 * hep)
      }
      yk = g_new - g
      if(yk %*% t(sk) > 0){
        bk = bk - (bk %*% t(sk) %*% sk %*% bk)/(sk %*% bk %*% t(sk)) +
          (t(yk) %*% yk)/(yk %*% t(sk))
      }
      else bk = bk
      theta = theta_new
    }
    m = m +1
  }
  f = fn(theta)
  g = g_new
  list(f = f,theta = theta,iter = m,g = g,h = bk)
}
