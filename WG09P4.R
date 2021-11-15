# Group_09, Chang LU, Hongxuan ZHU, Zheyi SHEN
# https://github.com/Zoey1102/BFGS_Optimizer

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
  fn <- function(theta) f(theta,...)
  m = 1
  hep = 1e-07
  gr = rep(0.0000000,length(theta))
  bk = diag(length(theta))
  while(m <= maxit){
    for(i in 1:length(theta)){
      pos = rep(0,length(theta))
      pos[i] = 1
      theta_add = theta + pos*hep
      theta_substract = theta - pos*hep
      gr[i] = (fn(theta_add) - fn(theta_substract))/(2 * hep)
    }
    if(sum(abs(gr)) < tol) break
    else{
      dk = -gr %*% solve(bk)
      a = 10
      lameda = 0.8
      c = 0.5
      b = TRUE
      while(b){
        if(fn(theta + a*dk)<fn(theta) + a * c * gr %*% t(dk)) b = FALSE
        else a = a * lameda
      }
      theta_new = theta + a * dk
      sk = theta_new - theta_new
      gr_new = rep(0,length(theta))
      for(i in 1:length(theta)){
        pos = rep(0,length(theta))
        pos[i] = 1
        theta_add = theta_new + pos*hep
        theta_substract = theta_new - pos*hep
        gr_new[i] = (fn(theta_add) - fn(theta_substract))/(2 * hep)
      }
      yk = gr_new - gr
      if(yk %*% t(sk) > 0){
        bk = bk - (bk %*% t(sk) %*% sk %*% bk)/(sk %*% bk %*% t(sk)) +
          (t(yk) %*% yk)/(yk %*% t(sk))
      }
      else bk = bk
      theta = theta_new
    }
    m = m +1
  }
  theta
}
