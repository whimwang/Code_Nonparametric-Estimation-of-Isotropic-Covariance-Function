#' ---
#' title: "BSpline"
#' author: "Whim Wang"
#' date: "7/23/2020"
#' output: html_document
#' ---
#' 
#' 
#' 
#' Boor Algorithm for calculating derivative of Bsplines and Bsplines
#' 
#' https://stackoverflow.com/questions/57507696/b-spline-derivative-using-de-boors-algorithm
#' 
#' 
#' Calculate tau value given index
#' @param index from -p to m+p+1
#' @return value $\tau(index)=index/(m+1)$
#' 
## -------------------------------------------------------------------------------------------------
tau_value <- function(index,p,m){
  if(index > m+p+1 || index < -p){
    stop('Index should be between -p and m+p+1 !')
  }else{
    return((index)/(m+1))
  }
}

#' 
#' 
#' This function calculate $f_j^l(x)$
#' $f_j^l(x)=\int_0^1(m+1)t^xB_{j+1}^l(t)dt$
#' 
#' @param x
#' @param j the jth basis,from 1 to m+p+1
#' @param l 
#' @param p order
#' @param m equally spaced knots with space $1/(m+1)$
#' @return value 
## -------------------------------------------------------------------------------------------------
Bspline_f <- function(x,j,l,p,m){
  if(l==0){
    if(0 <= tau_value(j-p,p,m) && tau_value(j-p+1,p,m) <= 1){ #tau_(-p)=t_1
      return((m+1)/(x+1)*(tau_value(j-p+1,p,m)^(x+1) 
                          - tau_value(j-p,p,m)^(x+1)))
    }else{
      return(0)
    }
  }else{
    return((m+1)/l * (Bspline_f(x+1,j,l-1,p,m) 
                      - tau_value(j-p,p,m) * Bspline_f(x,j,l-1,p,m)
                     + tau_value(j-p+l+1,p,m) * Bspline_f(x,j+1,l-1,p,m) 
                     - Bspline_f(x+1,j+1,l-1,p,m)))
  }
}


Cov_BsplineChoi <- function(rho,beta,p,m){
  if(length(beta)!=m+p){
    stop('vector beta should be of length m+p')
  }
  basis_value_matrix = matrix(0,nrow=length(rho),ncol=m+p)
  for(j in 1:(m+p)){
    basis_value_matrix[,j]=sapply(rho^2,Bspline_f,j=j,l=p-1,p,m)
  }
  return(as.vector(basis_value_matrix %*% beta))
}



Variog_BsplineChoi <- function(rho,beta,p,m){
  return(2*Cov_BsplineChoi(0,beta,p,m) - 2*Cov_BsplineChoi(rho,beta,p,m))
}


#' 
#' 
#' Plot in Choi's Paper
## ----eval=FALSE-----------------------------------------------------------------------------------
## h_seq = seq(from=0.1,to=2,length.out = 100)
## m=2;p=3
## for(j in seq(from=1,to=m+p,by=1)){
##   if(j==1){
##     plot(h_seq,sapply(h_seq^2, Bspline_f,j=j,l=2,p,m),
##          ylim=c(0,1),col=j,type='l',lty='dashed')
##   }else{
##     lines(h_seq,sapply(h_seq^2, Bspline_f,j=j,l=2,p,m),col=j,lty='dashed')
##   }
## }
## 
## #covariance 1
## #lines(h_seq,sapply(h_seq^2,function(x){5*Bspline_f(x,j=1,l=2,p,m)+0.5*Bspline_f(x,j=5,l=2,p,m)}))
## 
## #covariance 2
## #lines(h_seq,sapply(h_seq^2,function(x){0.7*Bspline_f(x,j=3,l=2,p,m)+0.1*Bspline_f(x,j=4,l=2,p,m)}))
## 
## #covariance 1 and covariance 2
## lines(h_seq,Cov_BsplineChoi(h_seq,beta=c(5,0,0,0,0.5),p,m))
## lines(h_seq,Cov_BsplineChoi(h_seq,beta=c(0,0,0.7,0.1,0),p,m))
## 
## 
## legend('topright',col=1:(m+p),legend=1:(m+p),lty=rep(1,m+p))
## 

#' 
## -------------------------------------------------------------------------------------------------
##knitr::purl('Bspline_Choi.Rmd',documentation = 2L)

