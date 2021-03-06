#' mhtboot: A package for multiple hypothesis testing using bootstrap distribution of p values.
#'
#' The mhtboot package provides three categories of important functions:
#' pboot, elbow and mht.
#' 
#' @section pboot functions:
#' pboot functions provide bootstrap distribution of p values. The pvalues are ordered and transformed. Currently the default transformation is fn(p) = -log(1-p) and in future some more transformations would be provided. There are support for two type of tests. One sample and two sample tests. The corresponding two functions are pboot.1sample and pboot.2sample. The test function by default is taken to be t.test(), while the user can provide their own test function. Both of these functions are parallelized using multicore for better performance.
#' @section elbow functions:
#' The purpose of elbow functions is to detect the change in distribution of the ordered transfromed p values. The basic function for detecting this change is elbow(), which takes in a particular p value curve and estimates the change point. We also provide a function to process the bootstrap distribution of p values and generate the estimate of the change point corresponding to a quantile of the empirical distribution.
#' @section mht:
#' The general function implementing the proceedure for multiple hypothesis testing based on bootstrap distribution of the p values. All the controls associated with pboot functions and elbow functions are transferred in mht functions too. There are two functions corresponding to one sample and two sample tests. These functions are mht.1sample and mht.2sample.
#'
#' @docType package
#' @name mhtboot
NULL