#' @title datgen
#' @description Function to generate data from multivariate normal with different mean.
#' @param n number of samples
#' @param m number of cords
#' @param m0 number of non sparse elements
#' @param sigeff magnitude of signal
#' @param Sigma Covariance matrix
#' @return X data matrix of size nxm
#' @details This function generates data from multivariate normal distribution with given covariance matrix. The mean values are either zero or constant sigeff, randomly permuted among the coordinates.
#' @examples 
#' \dontrun{
#' n = 50;m = 250;m0 = 20;
#' sigeff = 1;
#' Sigma <- 0.25*diag(m)
#' X <- datgen(n,m,m0,sigeff,Sigma = Sigma)
#' }
#' @export

datgen <- function(n,m,m0,sigeff,Sigma){
  nsindex <- sample(m,m0)
  meanvec <- numeric(m)
  meanvec[nsindex] <- sigeff
  e1.tmp <- eigen(Sigma)
  Sig.sqrt <- e1.tmp$vectors%*%tcrossprod(diag(sqrt(e1.tmp$values)),e1.tmp$vectors)
  X <- matrix(rnorm(n*m,mean=0,sd = 1),n,m)
  X <- X%*%Sig.sqrt
  X <- scale(X,center = meanvec)
  return(X)
}