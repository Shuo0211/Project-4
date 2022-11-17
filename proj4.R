#Name:
#Address of github repo:
#Contribution:
#Overview:
#The whole file is about an R function, newt, implementing Newton’s method for
#minimization of functions. In the program, we are mainly divided into six steps. 
#Step 1: evaluate function value, gradient and hessian matrix
#Step 2: test whether k-step theta is minimum and terminate if it is.
#Step 3: if k-step hessian matrix is not positive definite, perturb it.
#Step 4: solve H*Delta = -gradient for the search direction Delta.
#Step 5: if function value of next step is not less than this step,repeatedly 
#  halve Delta until it is.
#Step 6: set k+1-step theta = k-step theta + Delta, 
#  increment k by 1 and return to Step 1
#We multiply a small multiple (1e-10) of the norm of the Hessian matrix to perturb it
#And the project issue errors or warnings by try() function when the objective or 
#derivatives are not finite at the initial theta, the step fails to reduce 
#the objective, maxit is reached without convergence and the Hessian 
#is not positive definite at convergence.

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  #calculate f, gradient, Hessian at the initial theta
  f <- func(theta, ...)
  gradient <- grad(theta, ...)
  H <- hess(theta,...)
  
  #check Hessian and recalculate hessian by finite differencing if has no Hessian attribute
  if (any(is.null(H))) {
    #create a new matrix for Hessian
    H <- matrix(0, nrow = length(theta), ncol = length(theta))
    for (i in 1:length(theta)){
      th1 <- theta; th1[i] <- th1[i]+eps
      H <- grad(th1,...)
      H[i,] <- (H - gradient)/eps
    }
    #make Hessian matrix symmetric
    H <- 0.5 * (t(H) + H)
  }
  
  #check objective and derivatives are whether nor not finite
  if(is.finite(f)==FALSE || any(is.finite(gradient))==FALSE || any(is.finite(H))==FALSE){
    stop("The objective or derivatives are not finite at the initial theta.")
  }
  
  iter <- 0
  while(iter<=maxit){
    #convergence check
    while(max(abs(gradient)) > tol*(abs(f)+fscale)){
      per_iter <- 0
      #perturb Hessian if Hessian is not positive definite
      while(isTRUE(inherits(try(chol(H), silent = TRUE), "try-error"))){
        I <- diag(nrow = nrow(H))
        H <- H + I*norm(H)*eps*10^per_iter
        per_iter <- per_iter + 1
      }
      #calculate Delta
      Delta <- t(chol2inv(chol(H))%*%(-gradient))
      
      #try the max.half step halving
      half <- 0
      while (func(theta + Delta, ...)[1] > f[1]){
        if(half < max.half){
          Delta <- Delta/2
          half <- half + 1
        }else{
          stop("The step fails to reduce the objective despite trying 20 halvings.")
        }
      }
      
      #renew theta, f, gradient and Hessian
      theta <- theta + Delta
      f <- func(theta, ...)
      gradient <- grad(theta, ...)
      H <- hess(theta, ...)
      
      #check renew Hessian and recalculate Hessian if has no Hessian attribute
      if (any(is.null(H))) {
        #create a new matrix for Hessian
        H <- matrix(0, nrow = length(theta), ncol = length(theta))
        for (i in 1:length(theta)){
          th1 <- theta; th1[i] <- th1[i]+eps
          H <- grad(th1,...)
          H[i,] <- (H - gradient)/eps
        }
        #make Hessian matrix symmetric
        H <- 0.5 * (t(H) + H)
      }
      iter <- iter + 1
    }
    
    #check convergence while outside while loop
    #check whether or not Hessian is positive definite at convergence
    if (isTRUE(class(try(chol(H),silent=TRUE)) == "try-error")){
      warning("The Hessian is not positive definite at convergence.")  
    }else{
      #inverse Hessian
      Hi <- chol2inv(chol(H))
      newton_list <- list(f=f, theta=theta, iter=iter, g=gradient, Hi=Hi)
      return(newton_list)
    } 
  }
  #try the maximum iteration
  warning("Maximum 100 iterations is reached without convergence.")
}

