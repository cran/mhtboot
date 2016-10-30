#'@title Finding corner of a vector of ordered transformed p values
#'@param zvec vector of ordered transformed p values
#'@param rbuff scaler, by default 25. Controls the right buffer.
#'@param h scaler, default 30. Controls the window size.
#'@description Finds corner of a vector of ordered transformed p values.
#'@details The corner point of ordered p values indicate the point where the change from the alternative to null happens. So, by detecting that point we get an estimate of the number of true alternatives. 
#'
#'This function uses two methods for corner detection. One method is by transforming the vectors by taking their first difference and centering them around a theoretical mean for null case. The other method is by detecting the maximum change in gradient at each point. These methods will be denoted by dav and dlm respectively.
#'@return vector with two elements, containing estimates of the index of corner 
#'
#'        $dav: by average method.
#'        $dlm: by maximum gradient method.
#' @examples
#' \dontrun{         
#' n = 50;m = 250;m0 = 20;
#' sigeff = 1;
#' Sigma <- 0.25*diag(m)
#' X <- datgen(n,m,m0,sigeff,Sigma = Sigma)
#' porder <- pboot.1sample(X=X,B=500,ncpus = 1)
#' out <- elbow(zvec = porder[1,])
#' out
#' }
#'@export        
elbow <- function(zvec,rbuff=25,h=30){
  m <- length(zvec)
  cp <- which(zvec>0.05)[1] - 1
  emeans <- 1/(m+1-(1:cp))
  testdat <- data.frame(y=zvec[2:(cp+1)] - zvec[1:cp] - emeans,x=1:cp)
  smcut <- numeric(2)
  smcut[1] <- which.max(difave(testdat$y,rbuff = rbuff,h=h))
  smcut[2] <- which.max(diflm(testdat$x,testdat$y))
  #smcut[3] <- which.max(lcp(testdat))
  names(smcut) <- c("dav","dlm")
  smcut
}

#' @title Finding corner of a quantile of ordered transformed p values
#' @description Given a matrix of empirical distribution of ordered transformed p values, this function finds the corner point for a particular quantile.
#' @param porder matrix, usually feed from pboot functions. Bxm matrix of ordered p values, where B is the replication size and m is dimension.
#' @param rbuff right buffer, scaler, control for elbow()
#' @param h window size, default 30.
#' @param qi number between 0 and 1, quantile of the distribution. default 0.9.
#' @details In the distribution of the transformed ordered p values, we choose a particular quantile given by the user. We estimate the change point, which will be an estimate of the number of true alternatives corresponding to that quantile of the p values. As the values of the quantile increases, the estimates can only increasing, because we are dealing with ordered p values.
#' @return vector with two elements. estimates of the corner point by two methods.
#' @examples     
#' \dontrun{    
#' n = 50;m = 250;m0 = 20;
#' sigeff = 1;
#' Sigma <- 0.25*diag(m)
#' X <- datgen(n,m,m0,sigeff,Sigma = Sigma)
#' porder <- pboot.1sample(X=X,B=100,ncpus = 1)
#' out <- qelbow(porder = porder)
#' out
#' }
#' @export
qelbow <- function(porder,rbuff=25,h=30,qi=0.9){
  B <- nrow(porder)
  porder.df <- data.frame(porder)
  porder.df$ordstat <- factor(paste("z",1:B,sep=""))
  m.tmp <- ncol(porder.df)
  tmp2 <- apply(porder.df[,-m.tmp],2,function(x) quantile(x,probs=qi))
  smcut <- elbow(zvec = tmp2,rbuff = rbuff,h=h)
  return(smcut)
}

#' @title plotchange 
#' @param zvec vector of transformed order statistic of p values
#' @param rbuff right buffer
#' @param h window size
#' @param ... any graphical parameters passed to the plot function
#' @description Plot the change function that is maximized to find the change point.
#' @details Currently there are two types of change functions supported. The difference between first difference series and the difference in gradients at each point. Both of these functions should have a theoretical maximum at the change point. We plot these two series side by side along with indicating the change point.
#' @return Nothing
#' @examples   
#' \dontrun{      
#' n = 50;m = 250;m0 = 20;
#' sigeff = 1;
#' Sigma <- 0.25*diag(m)
#' X <- datgen(n,m,m0,sigeff,Sigma = Sigma)
#' porder <- pboot.1sample(X=X,B=100,ncpus = 1)
#' plotchange(porder[1,])
#' }
#' @export

plotchange <- function(zvec,rbuff = 25, h=30,...){
  m <- length(zvec)
  cp <- which(zvec>0.05)[1] - 1
  emeans <- 1/(m+1-(1:cp))
  testdat <- data.frame(y=zvec[2:(cp+1)] - zvec[1:cp] - emeans,x=1:cp)
  davseries <- difave(testdat$y,rbuff = rbuff,h=h)
  dlmseries <- diflm(testdat$x,testdat$y)
  par(mfrow=c(1,2))
  plot(davseries,type = "l",ylab = "dav")
  abline(v=which.max(davseries),col="red")
  plot(dlmseries,type = "l",ylab = "dlm")
  abline(v=which.max(dlmseries),col="red")
}

llm <- function(x,y,i){
  x <- x[1:(i-1)]
  y <- y[1:(i-1)]
  m <- lm(y~x)
  m
}

rlm <- function(x,y,i){
  x <- x[-(1:(i-1))]
  y <- y[-(1:(i-1))]
  m <- lm(y~x)
  m
}

diflm <- function(x,y,rbuff=25){
  n <- length(x)
  if(n < 2*rbuff) rbuff <- 0
  n <- n-rbuff
  out <- numeric(n)
  #   for(i in 1:n){
  #     out[i] <- predict(rlm(x,y,i),new=data.frame(x=x[i])) - predict(llm(x,y,i),new=data.frame(x=x[i]))
  #     #out[i] <- rlm(x,y,i)$coefficients[2]-llm(x,y,i)$coefficients[2]
  #   }
  out2 <- sapply(3:n, function(i) {newdat <- data.frame(x=x[i]);predict(rlm(x,y,i),new=newdat) - predict(llm(x,y,i),new=newdat)})
  out <- c(0,0,out2)
}



lave <- function(x,i,h=30){
  if(i > h){
    ilo <- (i-h)
  } else{
    ilo <- 1
  }
  mean(x[ilo:(i-1)])
}

rave <- function(x,i,h=30){
  n <- length(x)
  if(i < n+1-h){
    ihi <- i+h
  }else{
    ihi <- n
  }
  mean(x[i:ihi])
}

difave <- function(x,rbuff=25,h=30){
  n <- length(x)
  if(n < 2*rbuff) rbuff <- 0
  n <- n-rbuff
  out <- numeric(n)
  out <- sapply(1:n, function(i,h) rave(x,i,h)-lave(x,i,h),h=h)
  out
}

# lcp <- function(dat,h=10,deg=1,rbuff=25){
#   n <- nrow(dat)
#   if(n < 2*rbuff) rbuff <- 0
#   n <- n - rbuff
#   #dat <- dat[1:n,]
#   xev <- dat$x[1:n]
#   fitl <- locfit(y ~ left(x,h=h,deg=deg),ev=xev,data=dat)
#   fitr <- locfit(y ~ right(x,h=h,deg=deg),ev=xev,data=dat)
#   (predict(fitl,where = "ev") - predict(fitr,where = "ev"))^2
# }


