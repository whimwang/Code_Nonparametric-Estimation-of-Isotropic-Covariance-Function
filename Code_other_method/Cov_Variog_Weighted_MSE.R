


#' 
#' 
#' 
#' Functions calculates weighted mean square given type of covariance functions
#' $\sum w_i||\hat{C}(h_i)-C_{\theta}(h_i) ||^2$
#' @param theta, parameters for covariance function
#' @param h, vector
#' @param Cov_hat, empirical average of $C(h)$
#' @param weight
#' @return value of weighted mean square
#' 
#' 
## ------------------------------------------------------------------------------------------------
Spherical_WeightedMSE <- function(theta,h,Cov_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  return(sum((Cov_hat - sapply(h,Cov_Spherical,theta_1 = theta_1, theta_2 = theta_2))^2 * weight))
}

Gaussian_WeightedMSE <- function(theta,h,Cov_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  return(sum((Cov_hat - sapply(h,Cov_Gaussian,theta_1 = theta_1, theta_2 = theta_2))^2 * weight))
}

Matern_WeightedMSE <- function(theta,h,Cov_hat,weight,theta_3=NULL){
  theta_1 = theta[1]
  theta_2 = theta[2]
  theta_3 = theta_3
  return(sum((Cov_hat - sapply(h,Cov_Matern,theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3))^2 * weight))
}

TransMatern_WeightedMSE <- function(theta,h,Cov_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  theta_3 = theta[3]
  return(sum((Cov_hat - sapply(h,Cov_TransMatern,theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3))^2 * weight))
}

Cauchy_WeightedMSE <- function(theta,h,Cov_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  return(sum((Cov_hat - sapply(h,Cov_Cauchy,theta_1 = theta_1, theta_2 = theta_2))^2 * weight))
}

GenCauchy_WeightedMSE <- function(theta,h,Cov_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  theta_3 = 0.5
  theta_4 = 2
  return(sum((Cov_hat - sapply(h,Cov_GenCauchy,theta_1 = theta_1, theta_2 = theta_2,theta_3=theta_3, theta_4=theta_4))^2 * weight))
}

Wave_WeightedMSE <- function(theta,h,Cov_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  return(sum((Cov_hat - sapply(h,Cov_Wave,theta_1 = theta_1, theta_2 = theta_2))^2 * weight))
}



#' 
#' 
#' 
#' 
#' 
#' 
#' @param beta length p
#' @param y length n
#' @param X matrix n by p
#' @param weight length n
#' @return $$\sum_i w(y-X\beta)^2$$
#' 
## ------------------------------------------------------------------------------------------------
Weighted_MSE <- function(beta,y,X,weight){
  sum((y-X %*% beta)^2 * weight)
}

#' 
#' 
#' 
#' 
#' 
#' Functions calculates weighted mean square given type of covariance functions
#' $$\sum w_i||2\hat{\gamma}(h_i)-2\gamma(h_i) ||^2$$
#' @param theta, parameters for covariance function
#' @param h, vector
#' @param Cov_hat, empirical average of $C(h)$
#' @param weight
#' @return value of weighted mean square
#' 
#' 
## ------------------------------------------------------------------------------------------------
Spherical_Variog_WeightedMSE <- function(theta,h,Variog_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  return(sum((Variog_hat - sapply(h,Variog_Spherical,
                                  theta_1 = theta_1, theta_2 = theta_2))^2 * weight))
}

Gaussian_Variog_WeightedMSE <- function(theta,h,Variog_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  return(sum((Variog_hat - sapply(h,Variog_Gaussian,
                                  theta_1 = theta_1, theta_2 = theta_2))^2 * weight))
}

Matern_Variog_WeightedMSE <- function(theta,h,Variog_hat,weight,theta_3=NULL){
  theta_1 = theta[1]
  theta_2 = theta[2]
  theta_3 = theta_3
  return(sum((Variog_hat - sapply(h,Variog_Matern,
                                  theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3))^2 * weight))
}

TransMatern_Variog_WeightedMSE <- function(theta,h,Variog_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  theta_3 = theta[3]
  return(sum((Variog_hat - sapply(h,Variog_TransMatern,
                                  theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3))^2 * weight))
}

Cauchy_Variog_WeightedMSE <- function(theta,h,Variog_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  return(sum((Variog_hat - sapply(h,Variog_Cauchy,
                                  theta_1 = theta_1, theta_2 = theta_2))^2 * weight))
}


Wave_Variog_WeightedMSE <- function(theta,h,Variog_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  return(sum((Variog_hat - sapply(h,Variog_Wave,theta_1 = theta_1, theta_2 = theta_2))^2 * weight))
}


GenCauchy_Variog_WeightedMSE <- function(theta,h,Variog_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  theta_3 = 0.5
  theta_4 = 2
  return(sum((Variog_hat - sapply(h,Variog_GenCauchy,
                                  theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3, theta_4 = theta_4))^2 * weight))
}


#' 
#' 
#' 
#' Functions calculates weighted ratio of variog function
#' $$\sum w_i||2\hat{\gamma}(h_i)/2\gamma_\theta(h_i)-1 ||^2$$
#' @param theta, parameters for covariance function
#' @param h, vector
#' @param Cov_hat, empirical average of $C(h)$
#' @param weight
#' @return value of weighted mean square
#' 
#' 
## ------------------------------------------------------------------------------------------------
Spherical_VariogRatio_WeightedMSE <- function(theta,h,Variog_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  return(sum((Variog_hat/sapply(h,Variog_Spherical,
                                theta_1 = theta_1, theta_2 = theta_2)-1)^2 * weight))
}

Gaussian_VariogRatio_WeightedMSE <- function(theta,h,Variog_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  return(sum((Variog_hat - sapply(h,Variog_Gaussian,
                                  theta_1 = theta_1, theta_2 = theta_2)-1)^2 * weight))
}

Matern_VariogRatio_WeightedMSE <- function(theta,h,Variog_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  theta_3 = theta[3]
  return(sum((Variog_hat - sapply(h,Variog_Matern,
                                  theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3)-1)^2 * weight))
}

TransMatern_VariogRatio_WeightedMSE <- function(theta,h,Variog_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  theta_3 = theta[3]
  return(sum((Variog_hat - sapply(h,Variog_TransMatern,
                                  theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3)-1)^2 * weight))
}

Cauchy_VariogRatio_WeightedMSE <- function(theta,h,Variog_hat,weight){
  theta_1 = theta[1]
  theta_2 = theta[2]
  return(sum((Variog_hat - sapply(h,Variog_Cauchy,
                                  theta_1 = theta_1, theta_2 = theta_2)-1)^2 * weight))
}



#' 
#' 
#' 