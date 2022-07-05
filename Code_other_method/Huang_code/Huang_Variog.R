#' ---
#' title: "Variog_Huang."
#' author: "Whim Wang"
#' date: "7/30/2020"
#' output: html_document
#' ---
#' 
#' 
#' 
#' 
#' 
#' For Huang's code, 
#' Huang's paper: Sph,[0,24] with range 1
#' [0,20] with range 1 and 1.5 also works
#' 
#' Huang's paper: Mat,[0,24] with range 0.8
#' [0,20] with range 1/sqrt(s)
#' 
#' 
#' 
## ------------------------------------------------------------------------------------------------
Specify_CovFun <- function(option,par.cov){
  all_options=c('Sph','Gaussian','Cauchy','ChoiNP','Matern','LinearMatern','Wave','GenCauchy')
  if(!option %in% all_options){
    stop('Option of covariance function is not defined !')
  }else if(option == 'Cauchy'){ #convert 1 to 0.83
    arguments = list(theta_1 = par.cov[1],theta_2 = par.cov[2])
    Cov_Generation <- function(rho){Cov_Cauchy(rho,theta_1 = arguments$theta_1, 
                                                 theta_2 = arguments$theta_2)}
    
    Variog_Generation <- function(rho){Variog_Cauchy(rho,theta_1 = arguments$theta_1, 
                                                 theta_2 = arguments$theta_2)}
    
    Spectral_Generation <- function(w){Spectral_Cauchy(w,theta_1 = arguments$theta_1,
                                                        theta_2 = arguments$theta_2)}
    
  }else if(option == 'Matern'){ 
    matern_arguments = list(theta_1 = par.cov[1],theta_2 = par.cov[2],theta_3 = par.cov[3])
    arguments = list(theta_1=1,
                     theta_2 = matern_arguments$theta_2,
                     theta_3 = matern_arguments$theta_3)
    Cov_Generation <- function(rho){Cov_Matern(rho,theta_1 = arguments$theta_1, 
                                               theta_2 = arguments$theta_2,
                                               theta_3 = arguments$theta_3)}
    Variog_Generation <- function(rho){Variog_TransMatern(rho,theta_1 = arguments$theta_1, 
                                               theta_2 = arguments$theta_2,
                                               theta_3 = arguments$theta_3)}
    
    Spectral_Generation <- function(rho){Spectral_TransMatern(rho,theta_1 = arguments$theta_1, 
                                               theta_2 = arguments$theta_2,
                                               theta_3 = arguments$theta_3)}
  }else if(option == 'LinearMatern'){
    arguments = list(theta_1 = 1,theta_2 = 2.5,
                                theta_3 = 1,
                                theta_4 = 3, 
                                 theta_5 = 2)
    Cov_Generation <- function(rho){Cov_LinearTransMatern(rho,theta_1 = arguments$theta_1,
                                                     theta_2 = arguments$theta_2,
                                                     theta_3 = arguments$theta_3,
                                                     theta_4 = arguments$theta_4, 
                                                     theta_5 = arguments$theta_5)}
    
    Variog_Generation <- function(rho){Variog_LinearTransMatern(rho,theta_1 = arguments$theta_1,
                                                     theta_2 = arguments$theta_2,
                                                     theta_3 = arguments$theta_3,
                                                     theta_4 = arguments$theta_4, 
                                                     theta_5 = arguments$theta_5)}
    
    Spectral_Generation <- function(rho){Spectral_LinearTransMatern(rho,theta_1 = arguments$theta_1,
                                                     theta_2 = arguments$theta_2,
                                                     theta_3 = arguments$theta_3,
                                                     theta_4 = arguments$theta_4, 
                                                     theta_5 = arguments$theta_5)}
    
  }else if(option == 'Wave'){ #convert 0.5 to 5/12
    #arguments = list(theta_1 = 1,theta_2 = 0.5)
    arguments = list(theta_1 = 1,theta_2 = 0.75) #updated 0811
    Cov_Generation <- function(rho){Cov_Wave(rho,theta_1 = arguments$theta_1, 
                                               theta_2 = arguments$theta_2)}
    
    Variog_Generation <- function(rho){Variog_Wave(rho,theta_1 = arguments$theta_1, 
                                               theta_2 = arguments$theta_2)}
    
    Spectral_Generation <- function(rho){Spectral_Wave(rho,theta_1 = arguments$theta_1, 
                                               theta_2 = arguments$theta_2)}
  }else if(option == 'Sph'){
    #arguments = list(theta_1 = 1,theta_2 = 1)
    arguments = list(theta_1 = 1,theta_2 = 5) #updated 0811
    Cov_Generation <- function(rho){Cov_Spherical(rho,theta_1 = arguments$theta_1, theta_2 = arguments$theta_2)}
    Variog_Generation <- function(rho){Variog_Spherical(rho,theta_1 = arguments$theta_1, theta_2 = arguments$theta_2)}
    Spectral_Generation = NULL
    
  }else if(option == 'Gaussian'){
    #arguments = list(theta_1 = 1,theta_2 = 1)
    arguments = list(theta_1 = 1,theta_2 = 3) #updated 0813
    Cov_Generation <- function(rho){Cov_Gaussian(rho,theta_1 = arguments$theta_1, 
                                                 theta_2 = arguments$theta_2)}
    Variog_Generation <- function(rho){Variog_Gaussian(rho,theta_1 = arguments$theta_1, 
                                                 theta_2 = arguments$theta_2)}
    Spectral_Generation = NULL
    
  }else if(option == 'ChoiNP'){
    arguments = list(beta = rep(1/6,6),p=3,m=3)
    scale_value = Cov_BsplineChoi(0,beta = arguments$beta,p=arguments$p,m=arguments$m)
    arguments$beta = arguments$beta/scale_value
    
    Cov_Generation <- function(rho){Cov_BsplineChoi(rho,beta = arguments$beta,p=arguments$p,m=arguments$m)}
    Variog_Generation <- function(rho){Variog_BsplineChoi(rho,beta = arguments$beta,p=arguments$p,m=arguments$m)}
    Spectral_Generation <- NULL
    
  }else if(option == 'GenCauchy'){
    arguments = list(theta_1 = 1,theta_2 = 0.3,theta_3=0.5,theta_4=2)
    Cov_Generation <- function(rho){Cov_GenCauchy(rho,theta_1 = arguments$theta_1, 
                                                  theta_2 = arguments$theta_2, 
                                                  theta_3 = arguments$theta_3,
                                                  theta_4 = arguments$theta_4)}
    
    Variog_Generation <- function(rho){Variog_GenCauchy(rho,theta_1 = arguments$theta_1, 
                                                  theta_2 = arguments$theta_2, 
                                                  theta_3 = arguments$theta_3,
                                                  theta_4 = arguments$theta_4)}
    
    Spectral_Generation <- NULL
    
  }else{
     stop('error !')
  }
  arguments$option = option
  return(list(cov_function = Cov_Generation,
              variog_function = Variog_Generation,
              spectral_function = Spectral_Generation,
              arguments = arguments))
}


#' 
#' 
#' Calculate empirical variogram for observations for r independent replications of n locations
#' 
#' @param coord matrix,n by d
#' @param observation Y,r by n
#' @return list,sum_varioghat $2\hat{r}(h_i)\cdot N(h_i)$;n_h $N(h_i)$;dist $h_i$
#' 
## ------------------------------------------------------------------------------------------------
Huang_Empirical_Variogram <- function(Y,coord_matrix,dist_type='norm2'){
  r = nrow(Y); n=ncol(Y)
  if(nrow(coord_matrix)!=n){
    stop('col of Y is not equal to row of coord matrix')
  }
  
  
  #Y_mean = mean(Y)
  #Y = Y - Y_mean
  
  #nonzero dist=n*(n-1)/2 + zero dist=n
  result = matrix(0,nrow=n*(n-1)/2+1,ncol=3)
  
  ## nonzero dist
  start_idx = 1
  for(i in 1:(n-1)){
    #ith sample with (i+1):n sample
    end_idx = start_idx + (n-i)-1 #save result for (i,i+1) to (i,n),n-i in total
    #dist
    if(dist_type=='norm2'){
      dist = sapply(seq(from=i+1,to=n,by=1), function(idx){norm(coord_matrix[i,]-coord_matrix[idx,],type='2')}) 
    }else if(dist_type == 'distCosine'){
      dist = sapply(seq(from=i+1,to=n,by=1), function(idx){distm(coord_matrix[i,],coord_matrix[idx,],fun=distCosine)})
    }
    result[start_idx:end_idx,1] = dist
    #sum of variog_hat
    cov_hat = sapply(seq(from=i+1,to=n,by=1), function(idx){sum((Y[,i] - Y[,idx])^2)})
    result[start_idx:end_idx,2] = cov_hat
    #N_h
    num = rep(r,end_idx-start_idx+1)
    result[start_idx:end_idx,3] = num
    #
    #print(paste0('start: ',start_idx,' end: ',end_idx))
    start_idx = end_idx+1
  }
  
  ## zero dist
  result[start_idx,1] = 0
  result[start_idx,2] = 0
  result[start_idx,3] = n * r
  

   
  ## check unique dist
  unique_dist = unique(result[,1])
  summary_result = matrix(0,nrow=length(unique_dist),ncol=3)
  for(i in 1:length(unique_dist)){
    idx = which(result[,1]==unique_dist[i])
    if(length(idx)!=1){
      summary_result[i,] = c(unique_dist[i],apply(result[idx,2:3],2,sum))
    }else{
      summary_result[i,] = result[idx,]
    }
  }
  
  order_index = order(summary_result[,1])
  
  return(list(dist = summary_result[order_index,1],
              sum_varioghat = summary_result[order_index,2],
              n_h = summary_result[order_index,3]))
}


#' 
#' 
#' 
#' 
#' 
#' @param dist, vector of length n_fit,unique dist for fitting
#' @param v, scaler, range of spectral
#' @param L, num of grid
#' @param nugget, False when no nugget effect, True when nugget exists
#' @return matrix, n_fit by L
#' 
#' X matrix, n_dist by L, with (i,j) element 
#' $$X_{i,j}=\frac{v}{L}(1-J_0(w_i\cdot h_j))\frac{1+w_i^2}{w_i^2}$$
## ------------------------------------------------------------------------------------------------
Get_Bmatrix <- function(dist,v,L,nugget=FALSE){
  n_fit = length(dist)
  grid = seq(1,L)/L * v
  
  #v/L *(1-J_0(w*h_i))*r(w) 
  X <- matrix(0,n_fit,L) 
  for(i in 1:n_fit){ #
    X[i,] <- (1-besselJ(grid*dist[i],0))*(1+grid^2)/(grid^2)
  }
  X <- X * v/L 
  
  if(nugget){ #if nugget exists
    X = cbind(X,rep(1,n_fit))
  }
  return(X)
}

#' 
#' 
#' 
#' Calculate empirical variog at rho, $2\hat{\lambda}(h)$
#' @param rho
#' @param f_est
#' @param v
#' @param L
#' $$2\hat{\gamma}(h)\approx\frac{v}{L}\sum_i(1-J_0(w_ih))r(w_i)\hat{f}(w_i)$$
## ------------------------------------------------------------------------------------------------

Empirical_Variog_est <- function(dist,f_est,v,L,nugget=FALSE){
  sapply(dist, function(x){sum(Get_Bmatrix(dist=x,v=v,L=L,nugget=nugget) %*% f_est)})
}


#' 
#' @param v range
#' @param L num grid
#' @return K matrix in paper
#' 
## ------------------------------------------------------------------------------------------------
Get_PenalizeKmatrix <- function(v,L,nugget=FALSE){
#L_grid is the number of knots of w; 
# Qc is a matrix of J by J-2; 
# with each column, the diag is 1, the element under diag is -2 and the element next is 1
Qc <- matrix(0,L,L-2)
Qc[1,1] <- 1
Qc[2,1] <- -2
Qc[3,1] <- 1
for(j in 2:(L-2)){
  Qc[,j] <- c(0,Qc[1:(L-1),j-1])
}

Qc <- L*Qc/v
Rc <- matrix(0,L-2,L)
Rc[1,1] <- 1/6
Rc[1,2] <- 2/3
Rc[1,3] <- 1/6
for(i in 2:(L-2)){
  Rc[i,] <- c(0,Rc[i-1,1:(L-1)])
}
Rc <- Rc[,2:(L-1)]*v/L
Dc <- Qc%*%solve(Rc)%*%t(Qc) # This is (2.3) in Green and Silverman
#Dc is K matrix in paper for calculating gKg

if(nugget){ #if nugget exists
  Dc = cbind(Dc,rep(0,L))
  Dc = rbind(Dc,rep(0,L+1))
}
return(Dc)
}

#' 
#' 
#' 
#' Given svd decomposition, return $$\sum_i d_i^{a}u_i u_i^T$$
## ------------------------------------------------------------------------------------------------
gencov <- function(eva,eve,a){
  J <- dim(eve)[1]
  gpower <- matrix(0,J,J)
  q <- sum(eva > 1e-4)
  #	q <- sum(eva > 1e-2)
  for(i in 1:q){
    gpower <- gpower+eva[i]^a*eve[,i]%*%t(eve[,i])
  }
  gpower
}

#' 
#' 
#' 
#' Select best penalize lambda value for
#' @param dist, vector 
#' @param empirical_variog, vector
#' @param fit_ratio,value between 0 and 1
#' @param v, spectral range
#' @param L, num of grid
#' @param ngcv,number of lambda
#' 
#' @param X, matrix of n_fit by L
#' @param K, penalized matrix K of n_fit by n_fit
#' @param weight, vector of length n_fit
#' 
## ----eval=FALSE----------------------------------------------------------------------------------
## ## QP parameters
## function(dist,empirical_varioghat,GenCov_list,fit_ratio=0.75,
##          v=20,L=500,ngcv=20,do_plot=FALSE){
## 
## #dist = Huang_result_orig$dist
## #fit_ratio = 0.75
## 
## n_dist = length(dist)
## n_fit = as.integer(n_dist* fit_ratio)
## fitting_variog = empirical_variog[1:n_fit]
## 
## 
## n_fit = length(fitting_dist)
## X <- Xnew <- Get_Bmatrix(dist = fitting_dist,v=v,L=L,nugget = FALSE) #B in paper
## Dc <- Dnew <- Get_PenalizeKmatrix(v=v,L=L,nugget = FALSE)
## XnewWWXnew <- t(Xnew)%*%WW%*%Xnew
## weight = Huang_result_orig$n_h[1:n_fit]
## WW = diag(weight)
## dvec = t(Xnew) %*% WW %*% fitting_variog
## Amat = diag(L)
## bvec = rep(0,L)
## 
## QP_list = list(X = X, Xnew = Xnew, Dc = Dc, Dnew = Dnew, weight = weight, WW = WW,
##                dvec = dvec, Amat = Amat, bvec = bvec)
## 
## ngcv = 20
## lossgamma <- gcvnew <- tr1 <- tr2 <- numeric(ngcv)
## lambda_seq = sapply(1:ngcv, function(x){10^(6*(x-1)/(ngcv-1))})
## 
## 
## for(igcv in 1:ngcv){
## 
##   print(paste0('cross-validation: ',igcv))
##   print(paste0('lambda value: ',lambda_seq[igcv]))
##   lambda = lambda_seq[igcv]
## 
##   ### solve QP without nugget
##   Dmat <- XnewWWXnew + lambda*Dnew
##   QPout <- solve.QP(Dmat,dvec,Amat,bvec=bvec)
##   fest <- QPout$solution
##   #twosigma22 <- QPout$solution[J+1]#;print(twosigma22)
##   variog_hat <- Xnew%*%QPout$solution #in paper: Bv or say l_ij g_w
## 
##   if(do_plot){
##     if(igcv == 1){
##       plot(dist, empirical_varioghat,
##            main=paste0('Variog of ',GenCov_list$arguments$option))
##       points(dist, sapply(Huang_result_orig$dist,GenCov_list$variog_function),col='red',pch=20)
##       lines(dist, sapply(Huang_result_orig$dist,GenCov_list$variog_function),col='red')
##     }
##     lines(fitting_dist,variog_hat,col='red')
##   }
## 
##   ## Find constraints which are non-active constraints
##   index <- (fest < 1e-6)
##   remove <- which(index) #index non-active
##   keep <- which(!index) #index active
##   Xprime <- X[,keep] #B_tilde
##   Dcprime <- Dc[keep,keep]
## 
##   XprimeW = t(Xprime) %*% WW
##   ## A(lambda)
##   BprimeWBprime = t(Xprime)%*%WW%*%Xprime
##   A_lambda = Xprime %*% solve(BprimeWBprime + lambda * Dcprime) %*% XprimeW
## 
##   ## A(0)
##   eig <- eigen(t(Xprime)%*%WW%*%Xprime)
##   eva <- eig$values
##   eve <- eig$vectors
##   inv <- gencov(eva,eve,-1) ## Inv of (BprimeWBprime)
##   A_0 = Xprime %*% inv %*% XprimeW
##   p_value <- sum(diag(WW%*%A_0))  #p_value in paper
## 
##   gcvnew[igcv] <- sum(weight*(variog_hat-fitting_variog)^2)/((1-sum(diag(WW%*%A_lambda))/p_value)^2)
##   #lossf[igcv] <- sum((fest-ftrue)^2) #est spectral vs true spectral
##   lossgamma[igcv] <- sum((variog_hat-sapply(fitting_dist,GenCov_list$variog_function))^2)
## }
## 
## 
## best_igcv = which(gcvnew == min(gcvnew))
## best_lambda = lambda_seq[best_igcv]
## best_gamma = lossgamma[best_igcv]
## print(paste0('cross-validation: ',igcv,'  best lambda value: ', lambda_seq[best_igcv]))
## 
## ### solve QP first time
## Dmat <- t(Xnew)%*%WW%*%Xnew + best_lambda*Dnew
## QPout <- solve.QP(Dmat,dvec,Amat,bvec=bvec)
## fest <- QPout$solution
## gammahat <- Xnew%*%QPout$solution #in paper: Bv or say l_ij g_w
## 
## layout(matrix(1))
## #empirical variog
## d_seq = seq(from=0,to=max(Huang_result_orig$dist),length.out = 100)
## plot(Huang_result_orig$dist, Huang_result_orig$empirical_varioghat,
##      xlab='dist',ylab='variog value',
##      main=paste0('Huang method of ',GenCov_list$arguments$option))
## #true
## lines(Huang_result_orig$dist, sapply(Huang_result_orig$dist,GenCov_list$variog_function),col='black',lty='solid',lwd=2)
## #Huang's estimation
## lines(fitting_dist,gammahat,col='red',lwd=2)
## legend('topright',legend=c('empirical cov','true cov','Huang method'),
##        col = c('black','black','red'),
##        pch = c(20,NA,NA),
##        lty = c(NA,'solid','solid'))
## 
## 
##  lambda_metric = list(lambda = lambda_seq,
##              gcvnew = gcvnew,
##              lossgamma = gcvnew)
## 
##  solution_list = list(Huang_coef = QPout$solution,
##                       best_gamma = best_gamma,
##                       best_
##                       best_igcv = best_igcv,
##                       best_lambda = best_lambda)
##  return(list(
##             optout_Huang = QPout$solution,
##             QP_list = QP_list,
##             lambda_metric = lambda_metric,
##             ))
## 
## }
## 
## 

#' 
#' 
#' 
#' 
#' 
## ------------------------------------------------------------------------------------------------
##knitr::purl('Huang_Variog.Rmd',documentation = 2L)

