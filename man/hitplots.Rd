% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootfunctions.R
\name{hitplots}
\alias{hitplots}
\title{Plot area under p value cdf below a cutoff.}
\usage{
hitplots(porder, alpha = 0.005)
}
\arguments{
\item{porder}{the feed from porder.1sample or porder.2sample. matrix of size Bxm. of ordered transformed p values.}

\item{alpha}{the cutoff of ecd, by default 0.005.}
}
\description{
Function to plot the area under the cdf below a certain cutoff.
}
\details{
The alpha parameter specifies the cutoff, the plot is the ecdf under alpha. So the right tail of the ecdf would have probability alpha.
}
\examples{
\dontrun{
n = 50;m = 250;m0 = 20;
sigeff = 1;
Sigma <- 0.25*diag(m)
X <- datgen(n,m,m0,sigeff,Sigma = Sigma)
porder <- pboot.1sample(X, B = 100, ncpus = 1)
hitplots(porder)
}
}

