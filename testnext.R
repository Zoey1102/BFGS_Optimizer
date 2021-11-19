# Group_09, Chang LU, Hongxuan ZHU, Zheyi SHEN
# https://github.com/Zoey1102/BFGS_Optimizer

## Practical 4 Group 9 ####
## This practical is aimed to write a R function, bfgs, implementing the BFGS quasi-Newton minimization method.


## 1. Define the function to find gradients (either extract or compute)
fdg <- function(theta,f,...) {   
    ## Extract the gradient if it's available in f
    ## ... or compute it using definite difference
    if(is.null(attr(f(theta,...), 'gradient'))){ ## If gradient is not provided in f
        p <- length(theta); g <- rep(0,p) ## empty vector to store the gradient value
        eps <- 1e-7     ## finite difference interval
        for (i in 1:p) {    ## loop over parameters
            th_add <- theta; th_add[i] <- th_add[i] + eps ## increase theta[i] by eps
            f1 <- f(th_add,...) ## compute resulting f
            g[i] <- (f1 - f0)/eps   ## approximate df/dthta[i]
        }
    }
    else{
        g = attr(f(theta,...),"gradient") # Extract the gradient as it's available
    }
    g   ## Return the gradient
}



## 2. The BFGS quasi-Newton optimizer
bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100) {
    " theta: vector of initial optimization parameters
        f: objective function to minimize
      tol: convergence tolerance
   fscale: estimate of f at the optimum
    maxit: maximum number of iterations 
   return: a list of value at th minimum, including
           - objective function, parameters, iterations, gradients, approximate Hessian matrix
    "
    
    ## Bail out if the objective is infinite at the initial theta    
    f0 <- f(theta,...)    
    if (is.infinite(f0)) stop("infinite objective at the initial theta")
    
    ## Bail out if any derivatives are infinite at initial theta
    g0 <- fdg(theta,f,...)   
    if (sum(is.infinite(g0))) stop("infinite derivatives at the initial theta")
    
    
    ## Step trials   
    p <- length(theta)   ## take the dimension for convenience
    B <- diag(p) ## initialize approximate Hessian matrix using a p-dimensional Identity
    
    c1 <- 0.01 ## coefficient for sufficing a decrease
    c2 <- 0.9 ## coefficient for sufficing a positive definite matrix B
    maxit.s <- 20   ## maximum allowed times of halving the step in each iteration
    
    for (i in 1:maxit) {  ## loop until reaching the maximum number of iterations 
        g0 <- fdg(theta,f,...) ## call the fdg() to get the gradients
        f0 <- f(theta,...) ## return the initial value of objective function
        
        ##  --- theta is satisfactory if convergence is achieved ---
        if (max(abs(g0))<(abs(f0)+fscale)*tol) break   
        
        s <- -B %*% g0      ## ... or we will propose our 1st step by setting s = -B * g0, a p*1 vector
        for (k in 1:maxit.s) {  ## loop until reaching the maximum times of halving the step
            theta1 <- theta+s   ##  Take the proposed step
            g1 <- fdg(theta1,f,...) ## ... compute the resulting gradients
            f1 <- f(theta1,...)     ## ... and the resulting objective value
            
            
            ## Check Wolfe(1): if it's an improve/decrease on objectives
            if (f1<=f0+c1*g0%*%s ) {      ## if f is reduced
                y <- g1-g0    ## ... then compute the B matrix
                ro <- 1/as.numeric((t(s)%*%y))
                B <- (diag(p)-ro*s%*%t(y))%*%B%*%(diag(p)-ro*y%*%t(s))+ro*s%*%t(s)
                
                theta <- theta1     ## update the theta for next iteration
                H <- 0.5*(t(B)+B)   ## compute the approximate Hessian
                break
            } else s <- 0.5*s  ## halve the proposed step to reduce f
        }   ## end of the halve step loop
        
        if (k==maxit.s) {  ## If however "small" a step cannot reduce the objective function f, 
            warning("function not a normal convergence") ## ... report that the convergence criteria not met
            break
        }
        
        
        ## Check Wolfe(2): if matrix B stays positive definite
        if (g1%*%s>=c2*g0%*%s) {     ## if B is positive definite
            y <- g1-g0    ## ... then compute the B matrix
            ro <- 1/as.numeric((t(s)%*%y))
            B <- (diag(p)-ro*s%*%t(y))%*%B%*%(diag(p)-ro*y%*%t(s))+ro*s%*%t(s)
            
            theta <- theta1     ## update the theta for next iteration
            H <- 0.5*(t(B)+B)   ## compute the approximate Hessian
            break
        } else s <- 0.5*s  ## halve the proposed step to reduce f
    }   ## end of the halve step loop
    
    
    
    
    
    if (i==maxit) warning("iteration limit reached")
    list(f=f1,theta=theta,iter=i,g=g1,H=H)
}
