newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  f <- func(theta, ...)
  gradient <- grad(theta, ...)
  H <- hess(theta,...)
  if (any(is.null(H))) {
    H <- matrix(0L, nrow = length(theta), ncol = length(theta))
    for (i in 1:length(theta)){
      th1 <- theta; th1[i] <- th1[i]+eps
      H <- grad(th1,...)
      H[,i] <- (H - gradient)/eps
    }
    H <- 0.5 * (t(H) + H)
  }
  
  if(is.finite(f)==FALSE || any(is.finite(gradient))==FALSE || any(is.finite(H))==FALSE){
    stop("The objective or derivatives are not finite at the initial theta.")
  }
  iter <- 0
  while(iter<maxit){
    while(max(abs(gradient)) > tol*(abs(f)+fscale)){
      per_iter <- 0
      #perturb Hessian if Hessian is not positive definite
      while(isTRUE(inherits(try(chol(H), silent = TRUE), "try-error"))){
        I <- diag(nrow = nrow(H))
        H <- H + I*norm(H)*eps*10^per_iter
        
        #I <- diag(abs(max(H))*10^per_iter*eps, nrow=nrow(H), ncol=ncol(H))
        #H <- H + I
        per_iter <- per_iter + 1
      }
      # Calculate forward step towards optimum (minimum)
      #Delta <- chol2inv(chol(H), chol2inv(t(chol(H)), -gradient))
      Delta <- backsolve(chol(H), forwardsolve(t(chol(H)), -gradient))
      
      #try the max.half step halvings
      half <- 0
      while(is.finite(func(theta+Delta),...)==FALSE || any(is.finite(grad(theta+Delta),...))==FALSE ||any(is.finite(hess(theta+Delta),...))==FALSE){
        if(half < max.half){
          Delta <- Delta/2
          half <- half + 1
        }else{
          stop("The step fails to reduce the objective despite trying 20 steps halvings.")
        }
      }
      #renew theta, f, gradient and hessian
      theta <- theta + Delta
      f <- func(theta, ...)
      gradient <- grad(theta, ...)
      H <- hess(theta, ...)
      
      #check renew Hessian
      if (any(is.null(H))) {
        H <- matrix(0L, nrow = length(theta), ncol = length(theta))
        for (i in 1:length(theta)){
          th1 <- theta; th1[i] <- th1[i]+eps
          H <- grad(th1,...)
          H[,i] <- (H - gradient)/eps
        }
        H <- 0.5 * (t(H) + H)
      }
      iter <- iter + 1
    }
    if (max(abs(gradient)) < (abs(f)+fscale)*tol){
      if (isTRUE(class(try(chol(H),silent=TRUE)) == "try-error")){
        print("The Hessian is not positive definite at convergence")
      }else{
        Hi <- chol2inv(chol(H))
        newton_list1 <- list(f=f, theta=theta, iter=iter, g=gradient, Hi=Hi)
        return(newton_list1)
      } 
    }
  }
  stop("Maximum 100 iterations is reached without convergence.")
}

#test function
rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}


gb <-  function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}

hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}


Output1 <- newt(c(1, 1),rb,gb,hb)

