% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootfunctions.R
\name{pboot.1sample.s}
\alias{pboot.1sample.s}
\title{Generate p value distributions and estimate of sample correlation matrix using bootstrap.}
\usage{
pboot.1sample.s(X, B = 100, test = t.test, nbx = NROW(X), ncpus = 8,
  sout = FALSE)
}
\arguments{
\item{X}{data matrix}

\item{B}{Bootstrap size}

\item{test}{test to perform}

\item{nbx}{bootstrap sample size, by default same as the data sample size}

\item{ncpus}{number of cpus to use}

\item{sout}{if correlation matrix is needed or not}
}
\value{
a list with a matrix containing the p value distributions, and another matrix of correlation matrix.
}
\description{
If the user chooses to keep sout as TRUE, then this function generates bootstrap distribution of p values and returns the mean of the correlation matrices of all the bootstrap samples generated.
}
\examples{
\dontrun{
n = 50;m = 250;m0 = 20;
sigeff = 1;
Sigma <- 0.25*diag(m)
X <- datgen(n,m,m0,sigeff,Sigma = Sigma)
porder <- pboot.1sample.s(X=X,B=100,sout = TRUE,ncpus = 1)
plotpboot(porder)
} 
}

