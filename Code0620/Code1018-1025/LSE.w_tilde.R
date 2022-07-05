#' ---
#' title: "Initial_Solnp.Rmd"
#' author: "Yiming Wang"
#' date: "10/17/2019"
#' output: html_document
#' ---
#' 
## ------------------------------------------------------------------------
library("penalized")

#' 
#' 
#' 
#'                       
#' Least Square Estimation of $\tilde{w}$
#' 
#' Idea:
#' Step 1, Estimate Covariance of Y if Y is more than one replicates
#' Step 2, Transform to get All.Basis.Cor_matrix and Vectorize 
#' Step 3, do regression with or without positive constrains
#' Step 4, get w_tilde by scaling w by sum of w the
#' 
#' @param Y, a matrix; dim (r,n); assume r>1 and need to check
#' @param Log.All.Basis.Cor_matrix, a matrix, dim (n,n,m)
#' @param pos.coef, TRUE or FALSE, whether you put positive constraints on lse coef
#' @export
#' @return w_tilde
## ------------------------------------------------------------------------

LSE.w_tilde <- function(Y, Log.All.Basis.Cor_matrix, pos.coef){
  if(nrow(Y) == 1){
    stop("Y has no replicates and don't have LSE !")
  }
  cov_matrix = cov(Y, metho="pearson")
  cov_vector = Vectorize(cov_matrix)
  
  
  All.Basis.Cor_matrix = exp(Log.All.Basis.Cor_matrix)
  basis_vector = apply(All.Basis.Cor_matrix, 3, Vectorize) #each col is for one basis
  #no constrain
  lm(cov_vector~basis_vector-1)
  
  
  
  model.lse <- penalized(cov_vector, penalized = ~ basis_vector, unpenalized = ~ 0, 
                         lambda1=0, lambda2=0, positive = pos.coef, model = "linear", trace = FALSE)
  coef = as.vector(coefficients(model.lse, "penalized"))
  return(coef/sum(coef))
  
}

#' 
#' 
#' 
#' 
#' 
#' Vectorize a matrix
#' 
#' Keep upper triangle of a square matrix without diag and the first diag element 
#' 
#' @param a square matrix, (n,n)
#' @return a vector of length (n-1)*n/2+1
#' @export
## ------------------------------------------------------------------------
Vectorize <- function(A){
  keep.idx = upper.tri(A, diag = FALSE)
  keep.idx[1,1] = TRUE
  return(as.vector(A[keep.idx]))
}


