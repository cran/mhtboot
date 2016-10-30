#' @title Multiple hypothesis testing based on p value distribution for one sample test
#' @description Implements multiple hypothesis testing based on bootstrap distribution of p values.
#' @param X matrix of data
#' @param B bootstrap sample size, default is 100
#' @param test one sample test. by default t.test(), user can provide own function, must return p values in $p.value
#' @param nbx size of the bootstrap sample
#' @param ncpus number of cpu to use
#' @param rbuff right buffer for change detection
#' @param h window size for change detection
#' @param qi the quantile to use for change detection
#' @details This function takes the dataset and produces the bootstrap distribtution of the transformed and ordered p values using the user given parameters. Then detects the change in the bootstrap distribution using the corner detection method. This method requires the user to specify the quantile to use for change detection. The change point is an estimate of the location of change from alternative to null and used to get the coordinates of the true signals.
#' @return list with two elements. cutoff: the location of corner, signal: the index of the detected coordinates.
#' @examples 
#' n = 50;m = 100;m0 = 20;
#' sigeff = 1;
#' Sigma <- 0.25*diag(m)
#' X <- datgen(n,m,m0,sigeff,Sigma = Sigma)
#' out1 <- mht.1sample(X,B=100,ncpus = 1)
#' out1$cutoff
#' out1$signal
#' @export

mht.1sample <- function(X,B=100,test=t.test,nbx=NROW(X),
                        ncpus=8,rbuff=25,h=30,qi=0.9){
  porder <- pboot.1sample(X = X,B = B,test = test,nbx=nbx,ncpus = ncpus)
  out.elbow <- qelbow(porder = porder,rbuff = rbuff,h=h,qi=qi)
  p.orig <- apply(X,2,function(x) test(x)$p.value)
  cutoff <- min(out.elbow) #TODO apply other logic for cutoff selection here
  tmp <- sort.int(p.orig,index.return = TRUE)
  signifcoef <- which(tmp$ix < cutoff + 1)
  return(list(cutoff=cutoff,signal=signifcoef))
}