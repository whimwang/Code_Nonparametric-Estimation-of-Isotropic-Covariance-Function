#' ---
#' title: "Data_Generation_new.log.Rmd"
#' author: "Whim Wang"
#' date: "9/13/2019"
#' output: html_document
#' ---
#' 
#' 
#' 
#' 
#' #######################################################################################
#' 
#' 
#' #### Part II - Data Generation
#' 
#' 
## ------------------------------------------------------------------------
library("mvtnorm")

#' 
#' 
#' 
#' #### Description -- This function generate random coordinate.
#' 
#' #### Input - d = dimension
#' ####         n = number of locations
#' ####         l_min = lower bound of coordinates
#' ####         l_max = upper bound of coordinates
#' ####         seed 
#' 
#' #### Return - a n by d matrix
## ------------------------------------------------------------------------
Coordinates_Generation <- function(d,n,l_min,l_max,seed=0){
  set.seed(seed)
  coord_matrix = matrix(runif(n*d, min = l_min, max = l_max),nrow=n, ncol = d)

  return(coord_matrix)
}


#' 
#' 
#' 
#' 
#' 
#' #### Description -- This function generate Distance matrix based on coordinate matrix
#' #### Input - coord_matrix = a n by d matrix
#' #### Return - rho_matrix = a n by n matrix
#' 
#' 
## ------------------------------------------------------------------------
Distance <- function(coord_matrix){
  n = nrow(coord_matrix); d = ncol(coord_matrix)
  rho_matrix = matrix(0,n,n)

  #upper 
  for(i in 1:(n-1)){
    for(j in (i+1):n){
       rho_matrix[i,j] = norm(coord_matrix[i,]-coord_matrix[j,],type="2")
    }
  }

  #diag are zero
  
  #lower copy from upper
  #rho_matrix = (rho_matrix + t(rho_matrix))/2
  rho_matrix = rho_matrix + t(rho_matrix) #updated on 0921

  #diag is 0
  return(rho_matrix)
}


#' 
#' Function: Standardize Distance
#' Description: 
#' @param rho_matrix: a distance matrix or vector
#' @param value: a value which is used for transformation if specified; if it is not specified, will divide by max of rho_matrix
#' 
#' @return a list: transformed_rho_matrix and value which is ratio of rho_matrix and transformed_rho_matrix
## ------------------------------------------------------------------------
Transform.Distance <- function(rho_matrix, ratio = NULL){
  if(is.null(ratio)){#null
    ratio = max(rho_matrix)
  }
  return(list(transformed_rho_matrix = rho_matrix/ratio,
              ratio = ratio))
}

#' 


#' 
#' 
#' 
#' 
#' ### Description -- This function generate data Y and covariance matrix Sigma
#' 
#' #### Input -- Cov_function = the function returns covariance given distance rho
#' ####          rho_matrix = a n by n dist matrix, or say rho matrix
#' ####          r = number of replicates
#' 
#' #### Return -- a list: Y = generated observation
#' ####                   Sigma_matrix = True Sigma (Covariance) matrix based on rho_matrix and Cov_function
#' 
## ------------------------------------------------------------------------
library("MASS")
Data_Generation <- function(Cov_function,rho_matrix,r,option = "Exp",
                            par.cov = c(1,2),seed=0){
  set.seed(seed)
  n = nrow(rho_matrix)

  #Get Sigma matrix (covariance matrix)
  Sigma_matrix = matrix(sapply(rho_matrix,Cov_function,option,par.cov),nrow=n,ncol=n)
  
  
  #eigen.decompose = eigen(Sigma_matrix)
  #eigen.approx = sapply(eigen.decompose$values, function(x){ifelse(x<10^(-3),0,x)})
  #Sigma_R.approx = eigen.decompose$vectors %*% diag(eigen.approx) %*% t(eigen.decompose$vectors)
  
  
  # r by n matrix
  #Y = rmvnorm(r,rep(0,n),Sigma_matrix) #r*n matrix, use chol to decompose Sigma_matrix, cause error when n=500, not stable compared to eigen function
  Y = mvrnorm(r,rep(0,n),Sigma_matrix) # use eigen function in R
  #Y = mvrnorm(r,rep(0,n),Sigma_R.approx)
  

  return(list(Y=Y,Sigma_matrix=Sigma_matrix))
}




#' 
#' 
#' 
#' 
#' 
