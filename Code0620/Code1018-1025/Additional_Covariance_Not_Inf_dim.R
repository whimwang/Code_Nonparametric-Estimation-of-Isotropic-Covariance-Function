#' ---
#' title: "Additional_Covariance_Functions.Rmd"
#' author: "Yiming Wang"
#' date: "10/17/2019"
#' output: html_document
#' ---
#' 
#' 
#' Additional Isotropic Covariance Function: Not applicable in theory.
#' 
#' 
#' 
#' 
#' Circular Covariance Function
#' Applied in $R^2$ dimension.
#' @param a scalar, distance in $R^2$ or higher dimension
#' @param A scalar, covariance parameter
#' @param A scalar, range parameter
#' @return same dimension as rho
#' @example
#' 
#' plot of circular covariance function
#' 
#' rho_seq = seq(from=0,to=2,by=0.01)
#' plot(rho_seq, sapply(rho_seq, Cov_Circular, theta_1 = 1, theta_2 = 1),
#'      ylab="covariance value",xlab="distance", main="Circular")
#' 
## ------------------------------------------------------------------------
Cov_Circular <- function(rho, theta_1, theta_2){
  return(geoR::cov.spatial(rho,cov.model = 'circular',cov.pars = c(theta_1,theta_2)))
  
  #temp = rho/theta_2
  #return(ifelse(temp>1, 0, (acos(temp)-temp*sqrt(1-temp^2))*(2/pi)))
}



#' 
#' 
#' Spherical Covariance Function
#' Applied in $R^3$ or higher dimension
#' @param a scaler, distance
#' @param a scaler, covariance parameter
#' @param a scaler, range parameter
#' @return a scaler
#' @example
#' 
#' plot of spherical covariance function
#' 
#' rho_seq = seq(from=0, to=2, by=0.01)
#' plot(rho_seq, sapply(rho_seq, Cov_Spherical, theta_1=1, theta_2=1),
#'      ylab="Covariance Value", xlab="distance", main="Spherical")
#' 
#' 
## ------------------------------------------------------------------------
Cov_Spherical <- function(rho, theta_1, theta_2){
  return(geoR::cov.spatial(rho,cov.model = 'spherical',cov.pars = c(theta_1,theta_2)))
  #temp = rho/theta_2
  #return(ifelse(temp>1, 0, 1- temp*1.5 + 0.5*temp^3))
}



#' 
#' 
#' 
#' Wave Covariance Function
#' Applied in $R^3$ or higher dimension
#' @param a scalar, distance
#' @param a scaler, covariance parameter
#' @param a scaler, range parameter
#' @return a scalar
#' @example
#' 
#' Plot of Wave Covariance Function
#' 
#' rho_seq = seq(from=0, to=5, by=0.01)
#' plot(rho_seq, sapply(rho_seq, Cov_Wave, theta_1=1, theta_2=1),
#'      ylab="Covariance Value", xlab="distance", main="Wave",col="red")
#' abline(h=0)
#' 
## ------------------------------------------------------------------------
Cov_Wave <- function(rho, theta_1, theta_2){
  return(geoR::cov.spatial(rho,cov.model = 'wave',cov.pars = c(theta_1,theta_2)))
  #if(rho == 0){return(theta_1)}
  #temp = rho/theta_2
  #return(sin(temp)/temp)
}



#' 
