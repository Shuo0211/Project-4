## Members of group: Yan Chen s2318048, Shuo Wang s2439249, Yuxin Yang s2163400
## Address of github repo:https://github.com/Shuo0211/Project-4.git
## Contribution: @Shuo Wang and @Yuxin Yang take part in steps 1-4 given in the
## overview below, @Yan Chen takes part in steps 4-6 and the description of 
## the overview.
## Overview:
## The whole file is about an R function, newt, implementing Newton’s method for
## minimization of functions. In the program, we are mainly divided into six steps. 
## Step 1: evaluate function value, gradient and hessian matrix
## Step 2: test whether k-step theta is minimum and terminate if it is.
## Step 3: if k-step hessian matrix is not positive definite, perturb it.
## Step 4: solve H*Delta = -gradient for the search direction Delta.
## Step 5: if function value of next step is not less than this step,repeatedly 
##        halve Delta until it is.
## Step 6: set k+1-step theta = k-step theta + Delta, 
##        increment k by 1 and return to Step 1
## We multiply a small multiple (1e-10) of the norm of the Hessian matrix to
## perturb it.
## And the project issue errors or warnings by try() function when the objective
## or derivatives are not finite at the initial theta, the step fails to reduce 
## the objective, maxit is reached without convergence and the Hessian 
## is not positive definite at convergence.

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  ## Calculate f, gradient, Hessian at the initial theta
  f <- func(theta, ...) # the objective function to minimize at the initial theta
  gradient <- grad(theta, ...)# the gradient at the initial theta
  H <- hess(theta,...) # the Hessian matrix at the initial theta
  
  ## check Hessian and recalculate hessian by finite differencing if has no 
  ## Hessian attribute
  if (any(is.null(H))) {
    ## create a new matrix for Hessian
    H <- matrix(0, nrow = length(theta), ncol = length(theta))
    ## loop over parameters 
    for (i in 1:length(theta)){
      ## increase theta[i] by eps
      th1 <- theta; th1[i] <- th1[i]+eps
      H <- grad(th1,...) 
      ## approximate second derivs
      H[i,] <- (H - gradient)/eps
    }
    ## make Hessian matrix symmetric
    H <- (t(H) + H)/2
  }
  
  ## check objective and derivatives are whether nor not finite
  if(is.finite(f)==FALSE || any(is.finite(gradient))==FALSE || any(is.finite(H))==FALSE){
    # if they are not finite, stop and print the error information
    stop("The objective or derivatives are not finite at the initial theta.")
  }
  ## start iteration
  iter <- 0
  ## Check if the maximum number of loops has been reached.
  while(iter<=maxit){
    ## check convergence
    while(max(abs(gradient)) > tol*(abs(f)+fscale)){
      per_iter <- 0
      ## perturb Hessian if Hessian is not positive definite
      while(isTRUE(inherits(try(chol(H), silent = TRUE), "try-error"))){
        I <- diag(nrow = nrow(H)) # create a identity matrix
        ## multiply a small multiple (1e-10) of the norm of the Hessian matrix
        H <- H + I*norm(H)*eps*10^per_iter 
        ## Increment the number of iteration by one
        per_iter <- per_iter + 1 
      }
      ## calculate Delta = inv(H)*(-gradient)
      Delta <- t(chol2inv(chol(H))%*%(-gradient))
      
      ## try the max.half step halving
      half <- 0
      ## check if the minimum value of the function at the next step is 
      ## less than the previous step.
      while (func(theta + Delta, ...)[1] > f[1]){
        ## if the minimum value at the next step is greater than the previous step
        if(half < max.half){
          ## step halving
          Delta <- Delta/2
          ## Increment the number of iteration by one.
          half <- half + 1
        }else{
          ## If the step exceeds the max.half steps, the minimum value of the
          ## function in the next step is still greater than the previous step.
          ## we stop and print the error information.
          stop("The step fails to reduce the objective despite trying",max.half,"halvings.")
        }
      }
      
      ## renew theta
      theta <- theta + Delta
      ## renew f
      f <- func(theta, ...)
      ## renew gradient
      gradient <- grad(theta, ...)
      ## renew Hessian
      H <- hess(theta, ...)
      
      ## check renew Hessian and recalculate Hessian if has no Hessian attribute
      if (any(is.null(H))) {
        #create a new matrix for Hessian
        H <- matrix(0, nrow = length(theta), ncol = length(theta))
        ## loop over parameters
        for (i in 1:length(theta)){
          ## increase theta[i] by eps
          th1 <- theta; th1[i] <- th1[i]+eps
          H <- grad(th1,...)
          ## approximate second derivs
          H[i,] <- (H - gradient)/eps
        }
        ## make Hessian matrix symmetric
        H <- 0.5 * (t(H) + H)
      }
      ## update iteration
      iter <- iter + 1
    }
    
    ## check convergence while outside while loop
    ## check whether or not Hessian is positive definite at convergence
    if (isTRUE(class(try(chol(H),silent=TRUE)) == "try-error")){
      warning("The Hessian is not positive definite at convergence.")  
    }else{
      ## calculate inverse Hessian
      Hi <- chol2inv(chol(H))
      newton_list <- list(f=f, theta=theta, iter=iter, g=gradient, Hi=Hi)
      ## if convergence exists then return newton_list
      return(newton_list)
    } 
  }
  ## try the maximum iteration
  warning("Maximum 100 iterations is reached without convergence.")
}## newt
