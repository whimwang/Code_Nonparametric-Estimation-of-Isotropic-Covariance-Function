#' ---
#' title: "Sparse.w_tilde.Rmd"
#' author: "Yiming Wang"
#' date: "10/27/2019"
#' output: html_document
#' ---
#' 
#' 
#' #' Function
#' #' This function change small weights  in w_tilde to be zero.
#' @param: w_tilde, a vector of length m
#' @param: threshold.tol = 10^(-5)
#' #'Step 1: convert weights to zero, if smaller than threshold.tol
#' #'Step 2: convert to new weights so that sum is 1
#' 
#' 
## ------------------------------------------------------------------------
Sparse.weights <- function(w_tilde,threshold.tol=10^(-5)){
  w_tilde.new = w_tilde
  w_tilde.new[w_tilde < threshold.tol] = 0
  w_tilde.new = w_tilde.new/sum(w_tilde.new)
  return(w_tilde.new)
}

