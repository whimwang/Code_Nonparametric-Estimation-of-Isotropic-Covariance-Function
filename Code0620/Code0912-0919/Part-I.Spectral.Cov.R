#' ---
#' title: "Data_Generation_new.log.Rmd"
#' author: "Whim Wang"
#' date: "9/13/2019"
#' output: html_document
#' ---
library("geoR") # for matern function

#'Additional Covariance Function Not Applicable in Inf Dim
#source("\\\\wolftech.ad.ncsu.edu/cos/stat/Redirect/ywang225/Desktop/Exp_lambda_basis/2019September/Code1018-1025/Additional_Covariance_Not_Inf_dim.R")

#' 
#' 
#' #### Part I - Covariance function
#' #### Function 
#' 
#' #### Description -- This function computes the Exponential covariance and Gaussian covariancefunction 
#' 
#' #### Input -- rho = distance; distance in any dimension
#' ####          theta_1 
#' ####          theta_2
#' 
#' #### Exponential
#' $$C(\rho)=\theta_1 exp(-\frac{\rho}{\theta_2})$$
#' 
#' #### Gaussian
#' $$C(\rho)=\theta_1 exp(-(\frac{\rho}{\theta_2})^2)$$
#' 
#' #### Rational Quadratic
#' 
#' $$C(\rho)=\theta_1 (1-\frac{\rho^2}{\rho^2+\theta_2})=\theta_1\frac{1}{1+(\frac{h}{\theta_2})^2} $$
## ------------------------------------------------------------------------

#Exp
Cov_Exp<-function(rho,theta_1=1,theta_2=2){
  return(geoR::cov.spatial(rho,cov.model = 'exponential',cov.pars = c(theta_1,theta_2)))
  #return(theta_1*exp(-rho/theta_2))
}

#Gaussian
Cov_Gaussian<-function(rho,theta_1=1,theta_2=2){
  return(geoR::cov.spatial(rho,cov.model = 'gaussian',cov.pars = c(theta_1,theta_2)))
  #return(theta_1*exp(-(rho/theta_2)^2))
}


## RQ
Cov_RQ<-function(rho,theta_1=1,theta_2=2){
  return( theta_1/(1+(rho/theta_2)^2))
}


## Matern
Cov_Matern <- function(rho, theta_1,theta_2,theta_3){
  return(theta_1 * geoR::matern(u=rho, phi=theta_2, kappa = theta_3))
}

## GenCauchy
## special case: consider theta_3 = 0.5;theta_4=2
Cov_GenCauchy <- function(rho, theta_1,theta_2,theta_3,theta_4){
  return(cov.spatial(rho,cov.model = 'gencauchy',cov.pars = c(theta_1,theta_2),kappa=c(theta_3,theta_4)))
}


Cov_Cauchy <- function(rho,theta_1,theta_2){
  return(geoR::cov.spatial(rho,cov.model = 'cauchy',cov.pars = c(theta_1,theta_2,1/2)))
}


Cov_LinearMatern <- function(rho,theta_1,theta_2,theta_3,theta_4,theta_5){
  return(theta_1 * (0.5 * geoR::matern(u=rho, phi=theta_2, kappa = theta_3) + 
                      0.5 * geoR::matern(u=rho,phi=theta_4,kappa=theta_5)))
}



#' Use Basis function to Generate a sum of Basis function
#' @param distance, rho
#' @param theta_1 covariance value at rho=0
#' @param theta_2 range parameter
#' @param w_tilde, coeffcient of basis functions, sum equal to 1
#' @return value of covariance

Cov_sum_basis <- function(rho, theta_1, theta_2, w_tilde){
  if(sum(w_tilde)!=1){
    stop("sum of w_tilde are not equal to 1")
  }
  m = length(w_tilde)
  basis_value = sapply(seq(from=1,to=m,by=1), function(k)(return(exp(Log.Basis.Cor(rho/theta_2, k=k,m=m)))))
  return(as.numeric(w_tilde %*% basis_value) * theta_1)
}








#' 
#' 
#' 
#' #### Description -- This function return a final covariance function.
## ------------------------------------------------------------------------
# define Cov_function

Cov_function <- function(rho,option = "Exp",par.cov){
  if(option == "Exp"){
    return(Cov_Exp(rho,theta_1 = par.cov[1],theta_2 = par.cov[2]))
  }
  if(option == "Gaussian"){
    return(Cov_Gaussian(rho,theta_1 = par.cov[1],theta_2 = par.cov[2]))
  }
  if(option == "RQ"){
    return(Cov_RQ(rho,theta_1 = par.cov[1],theta_2 = par.cov[2]))
  }
  if(option == "Matern"){
    return(Cov_Matern(rho,theta_1 = par.cov[1],theta_2 = par.cov[2], theta_3 = par.cov[3]))
  }
  if(option == 'GenCauchy'){
    return(Cov_GenCauchy(rho,theta_1 = par.cov[1], theta_2 = par.cov[2], theta_3 = par.cov[3], theta_4 = par.cov[4]))
  }
  #======== sum of basis covariance function Applicable in Inf Dim===#
  if(option == "Basis"){
    #return(Cov_sum_basis(rho,theta_1 = par.cov[1], theta_2 = par.cov[2], w_tilde = par.cov[-(1:2)]))
    return(Cor.Approx.2d.Log(rho,w_tilde=par.cov))
  }
  if(option == 'Cauchy'){
    return(Cov_Cauchy(rho,theta_1 = par.cov[1], theta_2=par.cov[2]))
  }
  if(option == 'LinearMatern'){
    return(Cov_LinearMatern(rho,theta_1=par.cov[1],theta_2=par.cov[2],theta_3=par.cov[3],
                            theta_4=par.cov[4],theta_5=par.cov[5]))
  }
  
 
  
  
  #===== Additional Covariance Function Not Applicable in Inf Dim ========#
  if(option == "Wave"){
    return(Cov_Wave(rho, theta_1 = par.cov[1], theta_2 = par.cov[2]))
  }
  
  if(option == "Circular"){
    return(Cov_Circular(rho, theta_1 = par.cov[1], theta_2 = par.cov[2]))
  }
  
  if(option == "Spherical"){
    return(Cov_Spherical(rho, theta_1 = par.cov[1], theta_2 = par.cov[2]))
  }
  stop('cov of option is not defined')
}

#' 
#' 
#' 
#' 

#' 
#' 
#' 
#' Spectral Spherical of dimension=2
#' @param theta_1
#' @param theta_2
#' @return
#' $$f(w)=2\theta_1\frac{w^2e^{-\theta_2w}}{1+w^2}$$
#' 
## ------------------------------------------------------------------------------------------------
Spectral_Cauchy<-function(w,theta_1,theta_2){
  2 * theta_1 * w^2 * exp(-theta_2 * w)/(1+w^2)
}

#' 
#' Spectral matern of dim=2
#' @param theta_1
#' @param theta_2  range
#' @param smootheness
#' 
#' $$f(w)=4\theta_1\frac{\theta_3 w^3}{\theta_2^{2\theta_3}(w^2+\frac{1}{\theta_2^2})^{\theta_3+1}(1+w^2)}$$
#' 
## ------------------------------------------------------------------------------------------------
Spectral_Matern<-function(w,theta_1,theta_2,theta_3){
  tmp1 = (w^2+(1/theta_2)^2)^(theta_3+1)
  return(2 * theta_1 * (2 * theta_3 * w^3)/(theta_2^(2*theta_3))/tmp1/(1+w^2))
}

Spectral_LinearMatern<-function(w,theta_1,theta_2,theta_3,theta_4,theta_5){
  return(0.5 * Spectral_Matern(w, theta_1, theta_2, theta_3) + 
           0.5 *Spectral_Matern(w, theta_1, theta_4, theta_5))
}

#' 
#' Spectral Transmatern of dim=2
#' @param theta_1
#' @param theta_2  range
#' @param smootheness
#' 
#' 
## ------------------------------------------------------------------------------------------------
Spectral_TransMatern<-function(w,theta_1,theta_2,theta_3){
  return(Spectral_Matern(w,theta_1,theta_2=theta_2/2/sqrt(theta_3),theta_3))
}

#' 
#' Spectral LinearTransmatern of dim=2
#' @param theta_1
#' @param theta_2 range
#' @param theta_3 smoothness
#' @param theta_4
#' @param theta_5
## ------------------------------------------------------------------------------------------------
Spectral_LinearTransMatern<-function(w,theta_1,theta_2,theta_3,theta_4,theta_5){
  return(0.5 * Spectral_TransMatern(w, theta_1, theta_2, theta_3) + 
           0.5 *Spectral_TransMatern(w, theta_1, theta_4, theta_5))
}

#' 
#' 
#' Spectral Wave of dim=2
#' @param theta_1
#' @param theta_2 range parameter
#' 
#' $$f(w)=2\frac{\theta_1\theta_2^2 w^3}{(1-\theta_2^2w^2)^{\frac{1}{2}}(1+w^2)}$$
## ------------------------------------------------------------------------------------------------
Spectral_Wave<-function(w,theta_1,theta_2){
  
  ifelse(w>=1/theta_2,0, 2 * theta_1 * (theta_2^2) * (w^3)/(1+w^2)/((1-theta_2^2 * w^2)^(0.5)))
  
}



#' 
#' 
## ------------------------------------------------------------------------------------------------
Variog_Spherical <- function(rho, theta_1, theta_2){
  return(2 * theta_1 - 2*geoR::cov.spatial(rho,cov.model = 'spherical',cov.pars = c(theta_1,theta_2)))
}

Variog_Matern <- function(rho, theta_1,theta_2,theta_3){
  return(2 * theta_1 - 2*theta_1 * geoR::matern(u=rho, phi=theta_2, kappa = theta_3))
}

Variog_Exp <- function(rho, theta_1,theta_2,theta_3){
  return(2 * theta_1 - 2*geoR::cov.spatial(rho,cov.model = 'exponential',cov.pars = c(theta_1,theta_2)))
}

Variog_TransMatern <- function(rho, theta_1,theta_2,theta_3){
  return(2 * theta_1 - 2*theta_1 * geoR::matern(u=rho, phi=theta_2/2/sqrt(theta_3), kappa = theta_3))
}

Variog_Gaussian<-function(rho,theta_1,theta_2){
  return(2 * theta_1 - 2 * geoR::cov.spatial(rho,cov.model = 'gaussian',cov.pars = c(theta_1,theta_2)))
}

Variog_LinearMatern <- function(rho,theta_1,theta_2,theta_3,theta_4,theta_5){
  return(2*Cov_LinearMatern(0,theta_1,theta_2,theta_3,theta_4,theta_5) - 2*Cov_LinearMatern(rho,theta_1,theta_2,theta_3,theta_4,theta_5))
}

Variog_LinearTransMatern <- function(rho,theta_1,theta_2,theta_3,theta_4,theta_5){
  return(2*Cov_LinearTransMatern(0,theta_1,theta_2,theta_3,theta_4,theta_5) - 2*Cov_LinearTransMatern(rho,theta_1,theta_2,theta_3,theta_4,theta_5))
}


Variog_Cauchy <- function(rho,theta_1,theta_2){
  return(2*theta_1 -  2*(geoR::cov.spatial(rho,cov.model = 'cauchy',cov.pars = c(theta_1,theta_2,1/2))))
}

Variog_GenCauchy <- function(rho,theta_1,theta_2,theta_3,theta_4){
  return(2*theta_1 -  2*(geoR::cov.spatial(rho,cov.model = 'gencauchy',cov.pars = c(theta_1,theta_2),kappa=c(theta_3,theta_4))))
}


Variog_Wave <- function(rho,theta_1,theta_2){
  return(2*theta_1 - 2*Cov_Wave(rho,theta_1,theta_2))
}


Variog_RQ <- function(rho,theta_1,theta_2){
  return(2*theta_1 - 2*Cov_RQ(rho,theta_1,theta_2))
}





Variog_function <- function(rho,option = "Exp",par.cov){
  if(option == "Exp"){
    return(Variog_Exp(rho,theta_1 = par.cov[1],theta_2 = par.cov[2]))
  }
  if(option == "Gaussian"){
    return(Variog_Gaussian(rho,theta_1 = par.cov[1],theta_2 = par.cov[2]))
  }
  if(option == "RQ"){
    return(Variog_RQ(rho,theta_1 = par.cov[1],theta_2 = par.cov[2]))
  }
  if(option == "Matern"){
    return(Variog_Matern(rho,theta_1 = par.cov[1],theta_2 = par.cov[2], theta_3 = par.cov[3]))
  }
  if(option == 'GenCauchy'){
    return(Variog_GenCauchy(rho,theta_1 = par.cov[1], theta_2 = par.cov[2], theta_3 = par.cov[3], theta_4 = par.cov[4]))
  }
  #======== sum of basis covariance function Applicable in Inf Dim===#
  if(option == 'Cauchy'){
    return(Variog_Cauchy(rho,theta_1 = par.cov[1], theta_2=par.cov[2]))
  }
  if(option == 'LinearMatern'){
    return(Variog_LinearMatern(rho,theta_1=par.cov[1],theta_2=par.cov[2],theta_3=par.cov[3],
                            theta_4=par.cov[4],theta_5=par.cov[5]))
  }
  
  #===== Additional Covariance Function Not Applicable in Inf Dim ========#
  if(option == "Wave"){
    return(Variog_Wave(rho, theta_1 = par.cov[1], theta_2 = par.cov[2]))
  }
  if(option == "Spherical"){
    return(Variog_Spherical(rho, theta_1 = par.cov[1], theta_2 = par.cov[2]))
  }
  
  stop('variog of option is not defined')
  
}


Spectral_function <- function(rho,option = "Exp",par.cov){
  if(option == "Exp"){
    return(Spectral_Exp(rho,theta_1 = par.cov[1],theta_2 = par.cov[2]))
  }
  if(option == "Gaussian"){
    return(NULL)
  }
  
  if(option == "Matern"){
    return(Spectral_Matern(rho,theta_1 = par.cov[1],theta_2 = par.cov[2], theta_3 = par.cov[3]))
  }
  if(option == 'GenCauchy'){
    return(NULL)
  }
  #======== sum of basis covariance function Applicable in Inf Dim===#
  if(option == 'Cauchy'){
    return(Spectral_Cauchy(rho,theta_1 = par.cov[1], theta_2=par.cov[2]))
  }
  if(option == 'LinearMatern'){
    return(Spectral_LinearMatern(rho,theta_1=par.cov[1],theta_2=par.cov[2],theta_3=par.cov[3],
                               theta_4=par.cov[4],theta_5=par.cov[5]))
  }
  
  #===== Additional Covariance Function Not Applicable in Inf Dim ========#
  if(option == "Wave"){
    return(Spectral_Wave(rho, theta_1 = par.cov[1], theta_2 = par.cov[2]))
  }
  if(option == "Spherical"){
    return(NULL)
  }
  
  stop('spectral of option is not defined')
  
}






#' 
#' 
#' 
#' 