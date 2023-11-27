stage2 <- function(betahat,G,corY,n){ 
  n1 <- dim(G)[1]
  if(abs(det(t(betahat)%*% t(G) %*% G %*% betahat/n1)) < 1e-10)
  {
    return(0)
  }
  criterion <- abs(det(t(betahat)%*% t(G) %*% G %*% betahat/n1))
  theta = solve(t(betahat)%*% t(G) %*% G %*%betahat/n1, tol = 0) %*% t(betahat) %*% corY
  sigma2 = 
    (1 - 2*t(corY) %*% betahat %*% theta + crossprod(G %*% betahat %*% theta)/n1)
  sigma2 = as.numeric(sigma2)
  vartheta = solve(crossprod(G%*%betahat)/n1*n, tol = 0)*sigma2
  return(list(theta = theta,
              sigma2 = sigma2,
              vartheta = vartheta, 
              criterion = criterion))
}


stage2_ind <- function(betahat, Z1, Z2, Y){
  n2 <- length(Y)
  n1 <- dim(Z1)[1]
  if (abs(det(t(betahat) %*% t(Z1) %*% Z1 %*% betahat / n1)) < 1e-10) {
    return(0)
  }
  criterion <- abs(det(t(betahat) %*% t(Z1) %*% Z1 %*% betahat / n1))
  corY <- cor(Z2, Y)
  theta <- solve(t(betahat) %*% t(Z2) %*% Z2 %*% betahat / n2, tol = 0) %*% t(betahat) %*% corY
  sigma2 <-  (1 - 2 * t(corY) %*% betahat %*% theta + crossprod(Z2 %*% betahat %*% theta)/n2)
  sigma2 <- as.numeric(sigma2) 
  vartheta = solve(crossprod(Z2 %*% betahat), tol = 0) * sigma2
  return(list(theta = theta,
              sigma2 = sigma2,
              vartheta = vartheta, 
              criterion = criterion))
}
