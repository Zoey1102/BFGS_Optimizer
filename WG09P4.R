# Group_09, Chang LU, Hongxuan ZHU, Zheyi SHEN
# https://github.com/Zoey1102/BFGS_Optimizer

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
  # A BFGS method optimiser.
  # theta is a vector of initial values for the optimisation parameters.
  # f is the objective function to minimise. Called as f(theta, gradientLogical,...) 
  #   - theta is the vector of optimisation parameters.
  #   - gradientLogical is a logical indicating if gradients of the objectives w.r.t the parameters should be computed.
  #   - ... where all other arguments are passed
  #   - if gradientLogical = TRUE, f will have an attribute "gradient" (vector)
  # tol is the convergence tolerance.
  # fscale is a tough estimate of the magnitude of f at the minimum. Used in convergence testing.
  # maxit is the maximum number of BFGS iterations to try before giving up.
  f0 <- f(theta,...)
  
  
  
  list(f=as.numeric(f0),theta=theta,iter=i,g=g,H=)
}
