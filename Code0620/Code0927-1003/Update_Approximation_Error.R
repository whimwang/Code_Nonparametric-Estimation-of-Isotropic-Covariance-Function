#' ---
#' title: "RQ.Norm.Rmd"
#' author: "Whim Wang"
#' date: "9/26/2019"
#' output: html_document
#' ---
#' 
#' 
## ------------------------------------------------------------------------


##source("\\\\wolftech.ad.ncsu.edu/cos/stat/Redirect/ywang225/Desktop/Exp_lambda_basis/2019September/Code0912-0919/Part-III.Objective.functions.R")


#' 
#' 
#' 
#' 
#' 
#' 
#' #### Function
#' #### RQ.Diff.Matrix.before.integration.func
#' #### Description -- This function calculates the function before integration
#' 
#' $$
#' [\frac{Beta(k+\rho^2)}{Beta(k,m-k+1)}-R(\rho)]\times[\frac{Beta(j+\rho^2)}{Beta(j,m-j+1)}-R(\rho)]
#' $$
#' #### Input --  rho
#' ###            k,j,m
#' ###            method = "Lgamma" or "Lbeta" for Log.Basis.Cor approximation
#' 
## ------------------------------------------------------------------------

Diff.Matrix.before.integration.func <- function(rho,k,j, m, option, theta_2){
  
  Cor_seq = Cov_function(rho, option = option, par.cov = c(1,theta_2))
  A = exp(Log.Basis.Cor(rho, k, m)) - Cor_seq
  B = exp(Log.Basis.Cor(rho, j, m)) - Cor_seq
  return(A*B)
}

#' 
#' 
#' #### Function:  RQ.Diff.Matrix.ele.func
#' #### Description -- This function calculates the integration of function.
#' 
#' $$
#' D[k,j] =\int_0^\infty [\frac{Beta(k+\rho^2)}{Beta(k,m-k+1)}-R(\rho)]\times[\frac{Beta(j+\rho^2)}{Beta(j,m-j+1)}-R(\rho)]d\rho
#' $$
#' #### Input --  rho
#' ###            k,j,m
#' ###            method = "Lgamma" or "Lbeta" for Log.Basis.Cor approximation
#' 
## ------------------------------------------------------------------------
Diff.Matrix.ele.func <- function(k,j, m,option,theta_2){
  value = integrate(Diff.Matrix.before.integration.func,lower=0,upper = Inf,
                    k,j, m, option,theta_2)$value
  return(value)
}

#' 
#' 
#' 
#' 
#' #### Function:  RQ.Diff.Matrix.func
#' #### Description -- This function calculates a (m,m) Difference Matrix
#' 
#' #### Input -- m
#' ####          method="Lgamma" or "Lbeta"
#' ####          theta_2
#' 
#' #### Return-- m by m Diff matrix
## ------------------------------------------------------------------------
Diff.Matrix.func <- function(m,option,theta_2){
  
  Diff.matrix = matrix(0,nrow=m,ncol=m)
  #diag and upper
  for(k in seq(from=1,to=m,by=1)){
    #j th row; k th col
    Diff.matrix[seq(from=1,to=k,by=1),k] = sapply(seq(from=1,to=k,by=1), Diff.Matrix.ele.func,k,m, option,theta_2)
  }
  
  Diff.matrix = Diff.matrix + t(Diff.matrix)
  diag(Diff.matrix) = diag(Diff.matrix)/2
  return(Diff.matrix)
}

#' 
#' 
#' 
#' Test RQ.Diff.Matrix.func
## ----eval=FALSE----------------------------------------------------------
## m = 50; theta_2 = 2;
## 
## RQ.Diff.matrix <- RQ.Diff.Matrix.func(m,theta_2)
## 
## eigen(RQ.Diff.matrix)$values

#' 
#' 
#' 
#' 
##' #### Function
#' #### Description -- This function computes using log.Basis.Cor
#' $$\sum\tilde{w}_k\frac{Beta(k+\rho^2,m-k+1)}{Beta(k,m-k+1)}$$
#' Note this function is different from Cor.Approx.2d.Log defined before.
#' #### Input -- rho: distance;  w_tilde; method="Lgamma" or "Lbeta
## ------------------------------------------------------------------------
Cor.Approx.2d.Log.Check.Diff.Norm <- function(rho,w_tilde){
  m = length(w_tilde)
  
  Log.basis_func <- function(k){ 
    return(Log.Basis.Cor(rho=rho,k,m))
  }
  
  Log.basis = sapply(seq(from=1,to=m,by=1),Log.basis_func) 
  # Log.basis = sapply(1:5,Log.basis_func) can not be used here
  #return(sum(exp(Log.basis + log(w_tilde))))
  return(sum(exp(Log.basis) * w_tilde))
}


#' 
#' 
#' 
#' 
#' 
#' 
#' #### Function
#' #### Description -- This function plot the estimated correlation function and true correlation function
## ------------------------------------------------------------------------


Plot.Result.Check.Diff.Norm <- function(result.list,option, theta_2){
  
  rho_seq = seq(from=0, to=6, by=0.1)
  
  # true cor
  true.cor = sapply(rho_seq, Cov_function, option, par.cov=c(1,theta_2))
  plot(rho_seq, true.cor, pch=20,ylim=c(0,1),main= paste0(option, " theta2:",paste0(theta_2,collapse = ' ')))
  lines(rho_seq, true.cor, lty="dashed")
  
  # est cor
  # and plot
  m_seq = numeric(length(result.list))
  for(i in 1:length(result.list)){
    #for(i in 4){
    w_tilde = result.list[[i]]$est.par
    m_seq[i] = result.list[[i]]$m
    est.cor = sapply(rho_seq, Cor.Approx.2d.Log.Check.Diff.Norm, w_tilde)
    
    lines(rho_seq, est.cor,col=i+1)
  }
  
  #true again, for clarity: last layer is true
  lines(rho_seq, true.cor, lty="dashed")
  
  legend("topright",legend = paste0("m:",m_seq),col=1:length(m_seq)+1,
         lty=rep(1,length(m_seq)))
  
}



##L2 norm of True covariance

L_2norm <- function(option,par.cov){
  square_function <- function(rho,option,par.cov){
    Cov_function(rho,option,par.cov) * Cov_function(rho,option,par.cov)
  }
  temp = integrate(square_function, lower = 0, upper = Inf,option=option,par.cov=par.cov)$value
  return(sqrt(temp))
}



