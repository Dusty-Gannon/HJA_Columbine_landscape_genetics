###
# Functions to compute or create components 
# of the model for the sake of prediction and figures
###

#compute distance matrix among observations in matrix m
dist_matrix <- function(m){
  n <- nrow(m)
  ones <- rep(1,n)
  Gmat <- m%*%t(m) 
  diag_G <- diag(Gmat)
  return(ones%*%t(diag_G) -2*Gmat + diag_G%*%t(ones))
}




# construct the weights matrix given array X,
#  vector beta, matrix Z, and vector alpha
make_W <- function(size, X, beta, Z, alpha, tau){
  W_upper <- matrix(data=0, nrow = size, ncol = size)
  for(i in 1:size){
    for(j in i:size){
      W_upper[i,j] <- exp(X[i,j,]%*%beta + Z[i,]%*%alpha*tau)
    }
  }
  W_upper <- `diag<-`(W_upper, 0)
  return(W_upper + t(W_upper))
}




# create symmetric matrices with {i,j} pair means
sym_matrix <- function(v){
  m <- matrix(data=0, nrow = length(v), ncol = length(v))
  for (i in 1:length(v)) {
    for(j in i:length(v)){
      m[i,j] <- 0.5*(v[i] + v[j])
    }
  }
  m_zdiag <- m
  m_zdiag <- `diag<-`(m,0)
  return(m + t(m_zdiag))
}
