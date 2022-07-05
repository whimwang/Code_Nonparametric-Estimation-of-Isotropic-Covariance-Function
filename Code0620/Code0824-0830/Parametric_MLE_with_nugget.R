#' ---
#' title: "Parametric_MLE_nugget"
#' author: "Whim Wang"
#' date: "9/4/2020"
#' output: html_document
#' ---
#' 
#' 
#' Calculate -2loglik value of Gaussian
#' @param par,length 3, [1] nugget effect [2-3] theta_1, theta_2
#' @param Y
#' @param rho_matrix
#' @param approx.option
#' @param approx.par
## ---------------------------------------------------------------------------------------------------------------------------------------
Two.Neg.Log.Gaussian.given_pars.with.nugget <- function(par,Y,rho_matrix,approx.option=2,approx.par = 10^(-8)){
  nugget = par[1];
  par.cov = par[2:length(par)]
  Sigma = matrix(sapply(rho_matrix, Cov_Gaussian, theta_1 = par.cov[1], theta_2 = par.cov[2]),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Sigma = Sigma + diag(1,nrow=nrow(rho_matrix))*nugget
  Two.Neg.Log.Likelihood.Use.Sigma(Y,Sigma,approx.option = approx.option,approx.par = approx.par)
} 

Two.Neg.Log.Cauchy.given_pars.with.nugget <- function(par,Y,rho_matrix,approx.option=2,approx.par = 10^(-8)){
  nugget = par[1];
  par.cov = par[2:length(par)]
  Sigma = matrix(sapply(rho_matrix, Cov_Cauchy, theta_1 = par.cov[1], theta_2 = par.cov[2]),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Sigma = Sigma + diag(1,nrow=nrow(rho_matrix))*nugget
  Two.Neg.Log.Likelihood.Use.Sigma(Y,Sigma,approx.option = approx.option,approx.par = approx.par)
} 


Two.Neg.Log.Exp.given_pars.with.nugget <- function(par,Y,rho_matrix,approx.option=2,approx.par = 10^(-8)){
  nugget = par[1];
  par.cov = par[2:length(par)]
  Sigma = matrix(sapply(rho_matrix, Cov_Exp, theta_1 = par.cov[1], theta_2 = par.cov[2]),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Sigma = Sigma + diag(1,nrow=nrow(rho_matrix))*nugget
  Two.Neg.Log.Likelihood.Use.Sigma(Y,Sigma,approx.option = approx.option,approx.par = approx.par)
} 


Two.Neg.Log.Matern.given_pars.with.nugget <- function(par,Y,rho_matrix,theta_3 = NULL,approx.option=2,approx.par = 10^(-8)){
  nugget = par[1];
  par.cov = par[2:length(par)]
  Sigma = matrix(sapply(rho_matrix, Cov_Matern, 
                        theta_1 = par.cov[1], theta_2 = par.cov[2],theta_3 = theta_3),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Sigma = Sigma + diag(1,nrow=nrow(rho_matrix))*nugget
  Two.Neg.Log.Likelihood.Use.Sigma(Y,Sigma,approx.option = approx.option,approx.par = approx.par)
} 



Two.Neg.Log.GenCauchy.given_pars.with.nugget <- function(par,Y,rho_matrix,approx.option=2,approx.par = 10^(-8)){
  nugget = par[1];
  par.cov = par[2:length(par)]
  Sigma = matrix(sapply(rho_matrix, Cov_GenCauchy, 
                        theta_1 = par.cov[1], theta_2 = par.cov[2],theta_3 = 0.5,theta_4 = 2),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Sigma = Sigma + diag(1,nrow=nrow(rho_matrix))*nugget
  Two.Neg.Log.Likelihood.Use.Sigma(Y,Sigma,approx.option = approx.option,approx.par = approx.par)
} 



Two.Neg.Log.RQ.given_pars.with.nugget <- function(par,Y,rho_matrix,approx.option=2,approx.par = 10^(-8)){
  nugget = par[1];
  par.cov = par[2:length(par)]
  Sigma = matrix(sapply(rho_matrix, Cov_RQ, 
                        theta_1 = par.cov[1], theta_2 = par.cov[2]),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Sigma = Sigma + diag(1,nrow=nrow(rho_matrix))*nugget
  Two.Neg.Log.Likelihood.Use.Sigma(Y,Sigma,approx.option = approx.option,approx.par = approx.par)
} 


#' 
#' 
#' Parametric MLE with nugget
## ---------------------------------------------------------------------------------------------------------------------------------------

MLE_est_Exp_with_nugget <- function(par.init,Y,rho_matrix){
  if(length(par.init)!=3){
    stop('len of par init of Exp should be 3, nugget,theta_1,theta_2')
  }
  
  tmp = optim(par=par.init,fn = Two.Neg.Log.Exp.given_pars.with.nugget, method = "L-BFGS-B",
              lower = rep(10^(-6),3),
              Y = Y, rho_matrix = rho_matrix)
  Sigma = matrix(sapply(rho_matrix, Cov_Exp, 
                        theta_1 = tmp$par[2], theta_2 = tmp$par[3]),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Sigma = Sigma + diag(1,nrow=nrow(rho_matrix))*tmp$par[1]
  
  return(list(par = tmp$par,
              fvalue = tmp$value,
              Sigma = Sigma))
  
}


MLE_est_Cauchy_with_nugget <- function(par.init,Y,rho_matrix){
  if(length(par.init)!=3){
    stop('len of par init of Exp should be 3, nugget,theta_1,theta_2')
  }
  
  tmp = optim(par=par.init,fn = Two.Neg.Log.Cauchy.given_pars.with.nugget, method = "L-BFGS-B",
              lower = rep(10^(-6),3),
              Y = Y, rho_matrix = rho_matrix)
  Sigma = matrix(sapply(rho_matrix, Cov_Cauchy, 
                        theta_1 = tmp$par[2], theta_2 = tmp$par[3]),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Sigma = Sigma + diag(1,nrow=nrow(rho_matrix))*tmp$par[1]
  
  return(list(par = tmp$par,
              fvalue = tmp$value,
              Sigma = Sigma))
  
}



MLE_est_Gaussian_with_nugget <- function(par.init,Y,rho_matrix){
  if(length(par.init)!=3){
    stop('len of par init of Gaussian should be 3, nugget,theta_1,theta_2')
  }
  
  tmp = optim(par=par.init,fn = Two.Neg.Log.Gaussian.given_pars.with.nugget, method = "L-BFGS-B",
              lower = rep(10^(-6),3),
              Y = Y, rho_matrix = rho_matrix)
  
   Sigma = matrix(sapply(rho_matrix, Cov_Gaussian, 
                        theta_1 = tmp$par[2], theta_2 = tmp$par[3]),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Sigma = Sigma + diag(1,nrow=nrow(rho_matrix))*tmp$par[1]
  
  return(list(par = tmp$par,
              fvalue = tmp$value,
              Sigma = Sigma))
  
}



MLE_est_Matern_with_nugget <- function(par.init,Y,rho_matrix,theta_3 = NULL){
  if(is.null(theta_3)){#theta_3 is not given
    theta_3_seq = 1:10
    fvalue_seq = numeric(length(theta_3_seq))
    par_seq = matrix(0,nrow=length(theta_3_seq),ncol=3)
    
    for(i in 1:length(theta_3_seq)){
      tmp = optim(par=par.init,fn = Two.Neg.Log.Matern.given_pars.with.nugget, method = "L-BFGS-B",
              lower = rep(10^(-6),3),
              Y = Y, rho_matrix = rho_matrix,theta_3 = theta_3_seq[i])
      fvalue_seq[i] = tmp$value
      par_seq[i,] = tmp$par
    }
    index=which.min(fvalue_seq)
    
    Sigma = matrix(sapply(rho_matrix, Cov_Matern, 
                        theta_1 = par_seq[index,2], theta_2 = par_seq[index,3],theta_3 = theta_3_seq[index]),
                   nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
    Sigma = Sigma + diag(1,nrow=nrow(rho_matrix))*par_seq[index,1]
    
    return(list(par = c(par_seq[index,],theta_3_seq[index]),
                fvalue = fvalue_seq[index],
                Sigma = Sigma))
    
  }else{#theta_3 is given
     tmp = optim(par=par.init,fn = Two.Neg.Log.Matern.given_pars.with.nugget, method = "L-BFGS-B",
              lower = rep(10^(-6),3),
              Y = Y, rho_matrix = rho_matrix,theta_3 = theta_3)
     
     Sigma = matrix(sapply(rho_matrix, Cov_Matern, 
                        theta_1 = tmp$par[2], theta_2 = tmp$par[3],theta_3 = theta_3),
                   nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
     Sigma = Sigma + diag(1,nrow=nrow(rho_matrix))*tmp$par[1]
     return(list(par=c(tmp$par,theta_3),
                 fvalue = tmp$value,
                 Sigma = Sigma))
  }
}





MLE_est_GenCauchy_with_nugget <- function(par.init,Y,rho_matrix){
  if(length(par.init)!=3){
    stop('len of par init of GenCauchy should be 3, nugget,theta_1,theta_2')
  }
  
  tmp = optim(par=par.init,fn = Two.Neg.Log.GenCauchy.given_pars.with.nugget, method = "L-BFGS-B",
              lower = rep(10^(-6),3),
              Y = Y, rho_matrix = rho_matrix)
  
  Sigma = matrix(sapply(rho_matrix, Cov_GenCauchy, 
                        theta_1 = tmp$par[2], theta_2 = tmp$par[3],theta_3 = 0.5,theta_4=2),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Sigma = Sigma + diag(1,nrow=nrow(rho_matrix))*tmp$par[1]
  
  return(list(par = c(tmp$par,0.5,2),
              fvalue = tmp$value,
              Sigma = Sigma))
  
}




MLE_est_RQ_with_nugget <- function(par.init,Y,rho_matrix){
  if(length(par.init)!=3){
    stop('len of par init of RQ should be 3, nugget,theta_1,theta_2')
  }
  
  tmp = optim(par=par.init,fn = Two.Neg.Log.RQ.given_pars.with.nugget, method = "L-BFGS-B",
              lower = rep(10^(-6),3),
              Y = Y, rho_matrix = rho_matrix)
  
  Sigma = matrix(sapply(rho_matrix, Cov_RQ, 
                        theta_1 = tmp$par[2], theta_2 = tmp$par[3]),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Sigma = Sigma + diag(1,nrow=nrow(rho_matrix))*tmp$par[1]
  
  return(list(par = tmp$par,
              fvalue = tmp$value,
              Sigma = Sigma))
  
}
  
 



#' 
#' 
#' 
#' 
#' test
## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------
## option = 'Exp';par.cov = c(1,2);nugget = 0.2;
## 
## l_max=5
## Data.Preparation.result = Data_Preparation_Fixed_range(option,par.cov,d=2,N=2,r=200,n=40,l_max=l_max,seed = 0)
## 
## Y = Data.Preparation.result$Output.list$Y_matrix[1,,]
## rho_matrix = Data.Preparation.result$Output.list$rho_matrix.orig
## 
## Y = Y + rnorm(length(Y),mean=0,sd=sqrt(nugget))
## true_Sigma = Data.Preparation.result$Hidden.list$True.Sigma_matrix +diag(1,nrow=nrow(rho_matrix)) * nugget
## true_fvalue = Two.Neg.Log.Likelihood.Use.Sigma(Y,Sigma = true_Sigma,approx.option=2,approx.par = 10^(-8))
## 
## layout(matrix(1:2,1,2))
## plot(rho_matrix,true_Sigma_R)
## plot(rho_matrix,cov(Y))
## 
## 
## 
## MLE_est_Exp = MLE_est_Exp_with_nugget(par.init = c(0.8,0.8,0.8),Y,rho_matrix)
## MLE_est_Exp$par
## MLE_est_Exp$fvalue
## 
## 
## MLE_est_Gaussian = MLE_est_Gaussian_with_nugget(par.init = c(0.8,0.8,0.8),Y,rho_matrix)
## MLE_est_Gaussian$par
## MLE_est_Gaussian$fvalue
## 
## 
## MLE_est_Matern = MLE_est_Matern_with_nugget(par.init = c(0.8,0.8,0.8),Y,rho_matrix,theta_3 = NULL)
## MLE_est_Matern$par
## MLE_est_Matern$fvalue
## 
## 
## MLE_est_GenCauchy = MLE_est_GenCauchy_with_nugget(par.init = c(0.8,0.8,0.8),Y,rho_matrix)
## MLE_est_GenCauchy$par
## MLE_est_GenCauchy$fvalue
## 
## 
## MLE_est_RQ = MLE_est_RQ_with_nugget(par.init = c(0.8,0.8,0.8),Y,rho_matrix)
## MLE_est_RQ$par
## MLE_est_RQ$fvalue
## 
## 
## 
## plot(rho_matrix,true_Sigma,type='p',pch=20,ylim=c(0,2),
##      main = paste0(option,' Cov'))
## lines(rho_seq,sapply(rho_seq,Cov_function,option=option,par.cov = par.cov),lty='dashed',lwd=2)
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='Exp',par.cov=MLE_est_Exp$par[2:length(MLE_est_Exp$par)]),col='red',lwd=2)
## points(x=0,y=sum(MLE_est_Exp$par[1:2]),col='red',pch=16)
## #text(1,MLE_est_Exp$par[1],labels = round(MLE_est_Exp$par[1:2],digits=4))
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='Gaussian',par.cov=MLE_est_Gaussian$par[2:length(MLE_est_Gaussian$par)]),col='green',lwd=2)
## points(x=0,y=sum(MLE_est_Gaussian$par[1:2]),col='green',pch=16)
## #text(1,Gaussian_est$par[1],labels = round(Gaussian_est$par[1],digits=4))
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='Matern',par.cov=MLE_est_Matern$par[2:length(MLE_est_Matern$par)]),col='blue',lwd=2)
## points(x=0,y=sum(MLE_est_Matern$par[1:2]),col='blue',pch=16)
## #text(1,Matern_est$par[1],labels = round(Matern_est$par[1],digits=4))
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='GenCauchy',par.cov=MLE_est_GenCauchy$par[2:length(MLE_est_GenCauchy$par)]),
##       col='orange',lwd=2)
## points(x=0,y=sum(MLE_est_GenCauchy$par[1:2]),col='orange',pch=16)
## #text(1,GenCauchy_est$par[1],labels = round(GenCauchy_est$par[1],digits=4))
## 
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='RQ',par.cov=MLE_est_RQ$par[2:length(MLE_est_RQ$par)]),col='purple',lwd=2)
## points(x=0,y=sum(MLE_est_RQ$par[1:2]),col='purple',pch=16)
## #text(1,RQ_est$par[1],labels = round(RQ_est$par[1],digits=4))
## if(FALSE){
## legend('topright',legend=c('true','Exp_est','Gaussian_est','Matern_est','GenCauchy_est','RQ_est'),
##                   col = c('black','red','green','blue','orange','purple'),
##                   pch = c(20,rep(NA,5)),
##                   lty = c(NA,rep(1,5)))
## }
## 
## legend('topright',legend=c(paste0("true ",option, " (",paste0(round(c(nugget,par.cov),3),collapse = ','),") fvalue:",round(true_fvalue,3)),
## paste0("Exp_est (",paste0(round(MLE_est_Exp$par,3),collapse = ','),") fvalue:",round(MLE_est_Exp$fvalue,3)),
## paste0("Gaussian_est (",paste0(round(MLE_est_Gaussian$par,3),collapse = ','),") fvalue:",round(MLE_est_Gaussian$fvalue,3)),
## paste0("Matern_est (",paste0(round(MLE_est_Matern$par,3),collapse = ','),") fvalue:",round(MLE_est_Matern$fvalue,3)),
## paste0("GenCauchy_est (",paste0(round(MLE_est_GenCauchy$par,3),collapse = ','),") fvalue:",round(MLE_est_GenCauchy$fvalue,3)),
## paste0("RQ_est (",paste0(round(MLE_est_RQ$par,3),collapse = ','),") fvalue:",round(MLE_est_RQ$fvalue,3))
## ),
##                   col = c('black','red','green','blue','orange','purple'),
##                   pch = c(20,rep(NA,5)),
##                   lty = c(NA,rep(1,5)))
## 
## 
## 

#' 
#' 
#' 
#' 
#' 
## ----eval=FALSE-------------------------------------------------------------------------------------------------------------------------
## knitr::purl('Parametric_MLE_with_nugget.Rmd',documentation = 2L)

