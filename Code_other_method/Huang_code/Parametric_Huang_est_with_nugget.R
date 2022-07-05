#' ---
#' title: "Parametric_Huang_with_nugget"
#' author: "Whim Wang"
#' date: "10/5/2020"
#' output: html_document
#' ---
#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Gaussian_Variog_WeightedMSE_with_nugget <- function(theta,h,Variog_hat,weight){
  nugget = theta[1]
  theta_1 = theta[2]
  theta_2 = theta[3]
  return(sum((Variog_hat - 2*nugget - sapply(h,Variog_Gaussian,
                                  theta_1 = theta_1, theta_2 = theta_2))^2 * weight))
}

Matern_Variog_WeightedMSE_with_nugget <- function(theta,h,Variog_hat,weight,theta_3=NULL){
  nugget = theta[1]
  theta_1 = theta[2]
  theta_2 = theta[3]
  theta_3 = theta_3
  return(sum((Variog_hat -  2*nugget - sapply(h,Variog_Matern,
                                  theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3))^2 * weight))
}


Cauchy_Variog_WeightedMSE_with_nugget <- function(theta,h,Variog_hat,weight){
  nugget = theta[1]
  theta_1 = theta[2]
  theta_2 = theta[3]
  return(sum((Variog_hat -  2*nugget - sapply(h,Variog_Cauchy,
                                  theta_1 = theta_1, theta_2 = theta_2))^2 * weight))
}


GenCauchy_Variog_WeightedMSE_with_nugget <- function(theta,h,Variog_hat,weight){
  nugget = theta[1]
  theta_1 = theta[2]
  theta_2 = theta[3]
  theta_3 = 0.5
  theta_4 = 2
  return(sum((Variog_hat - 2*nugget- sapply(h,Variog_GenCauchy,
                                  theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3, theta_4 = theta_4))^2 * weight))
}



#' 
#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Matern_Parametric_est_with_nugget_Huang_model <- function(h,Variog_hat,weight){
  theta_3.vector <- seq(from=1,to=10,by=1)
  value.vector <- numeric(length(theta_3.vector))
  for(i in 1:length(theta_3.vector)){
    tmp <- try(optim(c(1,1,0.8),Matern_Variog_WeightedMSE_with_nugget,
                     h = h,
                     Variog_hat = Variog_hat,
                     weight = weight,
                     theta_3 = theta_3.vector[i],
                     lower=c(1e-6,1e-6),method="L-BFGS-B"),
               silent = TRUE)
    print(tmp$par)
    if('try-error' %in% class(tmp)){
      value.vector[i] = NA
    }else{
      value.vector[i] = tmp$value
    }
  }
  index = which.min(value.vector)
  optout_Matern <- try(optim(c(1,1,0.8),Matern_Variog_WeightedMSE_with_nugget,
                             h = h,
                             Variog_hat = Variog_hat,
                             weight = weight,
                             theta_3 = theta_3.vector[index],
                             lower=c(1e-6,1e-6),method="L-BFGS-B"),
                       silent = TRUE)
  optout_Matern$par = c(optout_Matern$par,theta_3.vector[index])
  return(optout_Matern)
}




## Gaussian #####
Gaussian_Parametric_est_with_nugget_Huang_model <- function(h,Variog_hat,weight){
  optout_Gaussian <- try(optim(c(1,1,3),Gaussian_Variog_WeightedMSE_with_nugget,
                               h = h,
                               Variog_hat = Variog_hat,
                               weight = weight,
                               lower=c(1e-6,1e-6),method="L-BFGS-B"),
                         silent = TRUE)
  if('try-error' %in% class(optout_Gaussian)){
    optout_Gaussian = list(par=NULL)
  }
  return(optout_Gaussian)
}


## Cauchy ###
Cauchy_Parametric_est_with_nugget_Huang_model <- function(h,Variog_hat,weight){
  optout_Cauchy <- try(optim(c(1,1,0.8),Cauchy_Variog_WeightedMSE_with_nugget,
                               h = h,
                               Variog_hat = Variog_hat,
                               weight = weight,
                               lower=c(1e-6,1e-6),method="L-BFGS-B"),
                         silent = TRUE)
  if('try-error' %in% class(optout_Cauchy)){
    optout_Cauchy = list(par=NULL)
  }
  return(optout_Cauchy)
}

## GenCauchy ###
GenCauchy_Parametric_est_with_nugget_Huang_model <- function(h,Variog_hat,weight){
  optout_GenCauchy <- try(optim(c(1,1,0.3),GenCauchy_Variog_WeightedMSE_with_nugget,
                             h = h,
                             Variog_hat = Variog_hat,
                             weight = weight,
                             lower=c(1e-6,1e-6),method="L-BFGS-B"),
                       silent = TRUE)
  if('try-error' %in% class(optout_GenCauchy)){
    optout_GenCauchy = list(par=NULL)
  }else{
    optout_GenCauchy$par = c(optout_GenCauchy$par,0.5,2)
  }
  return(optout_GenCauchy)
}


#' 
## ----eval=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------
## knitr::purl('Parametric_Huang_est_with_nugget.Rmd',documentation = 2L)

