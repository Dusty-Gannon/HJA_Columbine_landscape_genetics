### Functions used in edgeweight model fit and model checking ###

# DIC function using posterior draws
DIC_wishart <- function(S, estims, loglik){
  with(estims, {
    -2*(dWishart(S, df=kappa, Sigma = Sigma_hat, log = T) - 
          var(loglik))
  })
}



# Build contrast matrix 
contrast_mat <- function(n){
  L <- diag(x=1, nrow = n-1, ncol = n)
  Ljp1 <- diag(x=1,nrow = n-1, ncol = n-1)
  Ljp1 <- cbind(rep(0,n-1),Ljp1)
  L-Ljp1
}