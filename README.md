# BFGS_Optimizer
R function bfgs( ) implementing the BFGS quasi-Newton minimization method

Course name: MATH11176 Statistical Programming, University of Edinburgh

Group member: Zheyi SHEN, Chang LU, Hongxuan ZHU


# BFGS Overview
Working out Hessians is a bit tedious. Hugging a purpose of getting away with just gradients, the obvious way is to try to build an optimization method on the first order Taylor expansion, rather than the second order one. 
However, there are two problems with this. One is that the first order approximation does not have a minimum. The other is that the first order approximation gets worse as we approach the minimum of the objective, since the gradient is tending to zero, leaving ever diminishing grounds for having neglected the second order terms in the expansion. 
These theoretical problems matter in practice. 
The ‘steepest descent’ optimization method that results from the first order approximation is often extremely slow to converge near the optimum.

All is not lost however. 

A better method using only first derivatives is to combine gradient information over several steps to build up an approximation to the Hessian matrix. The fact that using any positive definite matrix in place of the Hessian in the Newton step still gives us a descent direction implies that it is not a big problem our Hessian will only be approximate. 

The way to build up the approximation is to insist that at each optimization step, the quadratic model implied by the approximate Hessian must be made to match the gradient evaluated at the start and the end of the step. 

This project is just the implementation of the most common BFGS22 version of the algorithm. It updates the inverse approximate Hessian, thereby avoiding the need for a matrix solve.
