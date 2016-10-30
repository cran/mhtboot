#' @title Generate Bootstrap Distribution of p values for one sample tests.
#' @param X matrix of data, each row is an independent observation nxm
#' @param B bootstrap sample size
#' @param test function for testing. default is t.test(). Must return a data frame with p value in $p.value.
#' @param nbx Sample size for the bootstrap samples. Default is NROW(X), which is same as the original data sample size.
#' @param ncpus Number of cpus to use for bootstrap. We use parallel:multicore() to parallelize the bootstrap. For windows, use ncpus = 1, for any other machine, you can use the maximum permissible number for your system.
#' @description Performs bootstrap to generate empirical distribution of order statistics of p values
#' @details We generate the bootstrap distribution of the order statistics of the p values. We are performing one sample test on each coordinate of the original dataset. The bootstrap used here is standard version with default bootstrap sample size being equal to data sample size. The default one sample test is t.test(), however the user can provide their own test functions. The only requirement is that it must return p values in $p.value column of the output.
#' The bootstrap is parallelized using multicore from the library parallel. Windows machines at this point does not support using multiple cores, so the ncpus option should be equal to 1 for windows. For other systems, it can be higher to speed up the process.
#' We also use a transofrmation of the p values, by default the transformation is -log(1-p). But the user can provide their own transformation function. They should be monotonically increasing functions.
#' @return matrix of dimension Bxm. (Where m coordinates), each row indicates transformed p values for that bootstrap sample.
#' @examples 
#' \dontrun{
#' n = 50;m = 250;m0 = 20;
#' sigeff = 1;
#' Sigma <- 0.25*diag(m)
#' X <- datgen(n,m,m0,sigeff,Sigma = Sigma)
#' porder <- pboot.1sample(X=X,B=100,ncpus = 1)
#' plotpboot(porder)
#' }
#' @export


pboot.1sample <- function(X,B=100,test=t.test,nbx=NROW(X),ncpus = 8){
  datX <- X
  nx <- NROW(X)
  tt <- function(i){
    bx <- sample.int(nx,size=nbx,replace = TRUE)
    X.b <- datX[bx,]
    pvec <- apply(X.b,2,function(x,nbx) test(x[1:nbx])$p.value,nbx=nbx)
    return(pvec)
  }
  res <- parallel::mclapply(X=seq_len(B),tt,mc.cores=ncpus)
  pmat <- do.call(rbind,res)
  porder <- apply(pmat, 1, sort)
  t(porder)
}


#' @title Generate p value distributions and estimate of sample correlation matrix using bootstrap.
#' @param X data matrix
#' @param B Bootstrap size
#' @param test test to perform
#' @param nbx bootstrap sample size, by default same as the data sample size
#' @param ncpus number of cpus to use
#' @param sout if correlation matrix is needed or not
#' @description If the user chooses to keep sout as TRUE, then this function generates bootstrap distribution of p values and returns the mean of the correlation matrices of all the bootstrap samples generated.
#' @return a list with a matrix containing the p value distributions, and another matrix of correlation matrix.
#' @examples 
#' \dontrun{
#' n = 50;m = 250;m0 = 20;
#' sigeff = 1;
#' Sigma <- 0.25*diag(m)
#' X <- datgen(n,m,m0,sigeff,Sigma = Sigma)
#' porder <- pboot.1sample.s(X=X,B=100,sout = TRUE,ncpus = 1)
#' plotpboot(porder)
#' } 
#' @export
pboot.1sample.s <- function(X,B=100,test=t.test,nbx=NROW(X),
                            ncpus = 8,sout=FALSE){
  datX <- X
  nx <- NROW(X)
  if(sout){
    tt <- function(i){
      bx <- sample.int(nx,size=nbx,replace = TRUE)
      X.b <- datX[bx,]
      pvec <- apply(X.b,2,function(x,nbx) test(x[1:nbx])$p.value,nbx=nbx)
      S.b <- cor(X.b)
      return(list(pvec=pvec,S=S.b))
    }
    res <- parallel::mclapply(X=seq_len(B),tt,mc.cores=ncpus)
    pmat <- t(mapply(function(x) x$pvec,res))
    porder <- apply(pmat, 1, sort)
    S <- Reduce("+",lapply(res,function(x)x$S))
    return(list(porder = t(porder), S.boot = S))
  } else {
    tt <- function(i){
      bx <- sample.int(nx,size=nbx,replace = TRUE)
      X.b <- datX[bx,]
      pvec <- apply(X.b,2,function(x,nbx) test(x[1:nbx])$p.value,nbx=nbx)
      return(pvec)
    }
    res <- parallel::mclapply(X=seq_len(B),tt,mc.cores=ncpus)
    pmat <- do.call(rbind,res)
    porder <- apply(pmat, 1, sort)
    return(list(porder = t(porder), S.boot = NULL))
  }
}

#' @title Transformation of order statistics of the p value distributions
#' @param porder matrix of p value order statistics, rows indicate replicates
#' @param trans one of ("default","normal","none") indicating trnasformation of -log(1-p), which is by default. Or inverse normal cdf transformation or no transformation.
#' @description This function applys transformation on the bootstrap distribution of order statistics of p values.
#' @details The transformation of p values must be monotonically increasing. The user can use their own transofrmation, however, this function supports only the commonly used transformations. These are -log(1-p) transformation, inverse normal cdf and identiy transformation.
#' @return matrix with transformed distribution.
#' @examples 
#' \dontrun{
#' X <- datgen(n=100,m=80,m0=20,sigeff=1,Sigma = 0.25*diag(80))
#' porder <- pboot.1sample(X=X,B=100,ncpus = 1)
#' porder.tr <- ptrans(porder,trans="normal")
#' plotpboot(porder.tr)
#' }
#' @export
ptrans <- function(porder,trans="default"){
  logtr <- function(mat){
    -log(1-mat)
  }
  
  phitr <- function(mat){
    qnorm(mat)
  }
  switch (trans,
    default = logtr(porder),
    normal = phitr(porder),
    none = porder
  )
}

#' @title Generate bootstrap distribution of p values based on user given two sample tests.
#' @param X matrix of data, each row is an independent observation nxm
#' @param Y matrix of data, sample 2. each row is an independent observation nxm.
#' @param B bootstrap sample size
#' @param test function for testing. default is t.test(). Must return a data frame with p value in $p.value.
#' @param nbx Sample size for the bootstrap samples. Default is NROW(X), which is same as the original data sample size.
#' @param nby Sample size for the bootstrap samples for second dataset. Default is NROW(X), which is same as the original data sample size.
#' @param ncpus Number of cpus to use for bootstrap. We use parallel:multicore() to parallelize the bootstrap. For windows, use ncpus = 1, for any other machine, you can use the maximum permissible number for your system.
#' @description Performs bootstrap to generate empirical distribution of order statistics of p values from two sample tests.
#' @details We generate the bootstrap distribution of the order statistics of the p values. We are performing one sample test on each coordinate of the original dataset. The bootstrap used here is standard version with default bootstrap sample size being equal to data sample size. The default one sample test is t.test(), however the user can provide their own test functions. The only requirement is that it must return p values in $p.value column of the output.
#' The bootstrap is parallelized using multicore from the library parallel. Windows machines at this point does not support using multiple cores, so the ncpus option should be equal to 1 for windows. For other systems, it can be higher to speed up the process.
#' We also use a transofrmation of the p values, by default the transformation is -log(1-p). But the user can provide their own transformation function. They should be monotonically increasing functions.
#' @return matrix of dimension Bxm. (Where m coordinates), each row indicates transformed p values for that bootstrap sample.
#' @examples 
#' \dontrun{
#' n = 50;m = 250;m0 = 20;
#' sigeff = 1;
#' Sigma <- 0.25*diag(m)
#' X <- datgen(n,m,m0,sigeff,Sigma = Sigma)
#' porder <- pboot.1sample(X=X,B=100,ncpus = 1)
#' plotpboot(porder) 
#' }
#'@export
pboot.2sample <- function(X,Y,B=100,test=t.test,nbx=NROW(X),nby=NROW(Y),
                          ncpus = 8){
  datX <- X
  datY <- Y
  nx <- NROW(X)
  ny <- NROW(Y)
  # dat <- rbind(datX,datY)
  tt <- function(i){
    bx <- sample.int(nx,size=nbx,replace = TRUE)
    X.b <- datX[bx,]
    by <- sample.int(ny,size=nby,replace = TRUE)
    Y.b <- datX[by,]
    dat <- rbind(X.b,Y.b)
    pvec <- apply(dat,2,function(x,nbx,nby) test(x[1:nbx],x[(nbx+1):nby])$p.value,nbx=nbx,nby=nby)
    -log(1-pvec)
  }
  res <- parallel::mclapply(X=seq_len(B),tt,mc.cores=ncpus)
  pmat <- do.call(rbind,res)
  porder <- apply(pmat, 1, sort)
  t(porder)
}

#' Plot function for pboot
#' @title Quantile plots for p value distributions.
#' @param porder Matrix feeds from pboot. This is a matrix of p values from the bootstrap sampls. Of size Bxm, each row for one bootstrap. The columns indicate the coordinates for testing. 
#' @description Produces density plots of quantiles of transformed order statistics of p values
#' @details This function plots the order statistics of the quantiles of the transformed p values. As the distribution of the statistic changes as the number of coordinates increase, it should show a change in the curve. 
#' 
#' This function uses ggplot2 and reshape library to manipulate data. The final object returned is a ggplot2 image that can be fed into ggsave or any other supported functions.
#' @return ggplot2 object contatining the plot.
#' @examples 
#' \dontrun{
#' n = 50;m = 250;m0 = 20;
#' sigeff = 1;
#' Sigma <- 0.25*diag(m)
#' X <- datgen(n,m,m0,sigeff,Sigma = Sigma)
#' Y <- datgen(n,m,m0,sigeff,Sigma = Sigma)
#' porder <- pboot.2sample(X=X,Y = Y, B=100,ncpus = 1)
#' plotpboot(porder)
#' }
#' @export
plotpboot <- function(porder){
  niter <- nrow(porder)
  m <- NCOL(porder)
  porder.df <- data.frame(porder)
  porder.df$ordstat <- factor(paste("z",1:niter,sep=""))
  tmp3 <- apply(porder.df[,1:m],2,function(x) quantile(x,probs=c(1:9)/10))
  tmp3 <- t(tmp3)
  cpoint <- which(tmp3[,5]>0.05)[1] - 1
  tmpsmall <- tmp3[1:cpoint,]
  tmpsmall.df <- data.frame(tmpsmall)
  colnames(tmpsmall.df) <- colnames(tmpsmall)
  tmpsmall.df$index <- 1:cpoint
  tmpsmall.melt <- reshape2::melt(tmpsmall.df,variable.name = "percentile",id.vars = "index")
  colnames(tmpsmall.melt) <- c("index", "percentile", "value")
  gp4 <- ggplot2::ggplot(tmpsmall.melt,ggplot2::aes_string(x='index',y='value',group='percentile',color='percentile'))+
    ggplot2::geom_line() + ggplot2::geom_hline(aes(yintercept=0.05))+ggplot2::theme_bw()
  gp4
}

#' @title Plot area under p value cdf below a cutoff.
#' @description Function to plot the area under the cdf below a certain cutoff.
#' @details The alpha parameter specifies the cutoff, the plot is the ecdf under alpha. So the right tail of the ecdf would have probability alpha.
#' @param porder the feed from porder.1sample or porder.2sample. matrix of size Bxm. of ordered transformed p values.
#' @param alpha the cutoff of ecd, by default 0.005.
#' @examples 
#' \dontrun{
#' n = 50;m = 250;m0 = 20;
#' sigeff = 1;
#' Sigma <- 0.25*diag(m)
#' X <- datgen(n,m,m0,sigeff,Sigma = Sigma)
#' porder <- pboot.1sample(X, B = 100, ncpus = 1)
#' hitplots(porder)
#' }
#' @export
hitplots <- function(porder,alpha=0.005){
  pcdfs <- apply(porder,2,function(x){y<-ecdf(x);y(alpha)})
  plot(pcdfs[pcdfs>0.05],type="l")
  cp.hit <- which(pcdfs<1)[1]
  abline(v=cp.hit,col="red")
  return(cp.hit)
}
