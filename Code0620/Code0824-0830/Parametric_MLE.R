#' ---
#' title: "Parametric_MLE"
#' author: "Whim Wang"
#' date: "8/27/2020"
#' output: html_document
#' ---
#' 
#' Calculate -2loglik for Gaussian 
#' Given range parameter, first calculate Correlation function; then calculate $\hat{\sigma^2}$ and corresponding -2loglik
#' 
#' @param:theta:range parameter and
#' @param:Y,r by n matrix
#' @param:rho_matrix, n by n matrix
#' @param:theta_1:$\sigma^2$, default is 1,
#' @return
#' 
#' 
#' Calculate -2loglik for Matern
#' @param theta,length 1 vector; range value
#' @param:Y
#' @param:rho_matrix
#' @param:theta_1 $\sigma^2$
#' @param:theta_3 $kappa$ need to be given
#' 
#' 
#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------

Two.Neg.Log.Gaussian.given_pars <- function(theta,Y,rho_matrix,theta_1=1,approx.option = 2,approx.par = 10^(-8)){
  Sigma_R = matrix(sapply(rho_matrix, Cov_Gaussian, theta_1 = theta_1, theta_2 = theta),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Two.Neg.Log.Likelihood.Use.Sigma_R(Y,Sigma_R,approx.option=approx.option,approx.par=approx.par)
}



Two.Neg.Log.Exp.given_pars <- function(theta,Y,rho_matrix,theta_1=1,approx.option = 2,approx.par = 10^(-8)){
  Sigma_R = matrix(sapply(rho_matrix, Cov_Exp, theta_1 = theta_1, theta_2 = theta),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Two.Neg.Log.Likelihood.Use.Sigma_R(Y,Sigma_R,approx.option=approx.option,approx.par=approx.par)
}


Two.Neg.Log.Matern.given_pars <- function(theta,Y,rho_matrix,theta_1=1,theta_3=NULL,approx.option = 2,approx.par = 10^(-8)){
  Sigma_R = matrix(sapply(rho_matrix, Cov_Matern, theta_1 = theta_1, theta_2 = theta[1],theta_3=theta_3),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Two.Neg.Log.Likelihood.Use.Sigma_R(Y,Sigma_R,approx.option=approx.option,approx.par=approx.par)
}


Two.Neg.Log.GenCauchy.given_pars <- function(theta,Y,rho_matrix,theta_1 = 1,approx.option = 2,approx.par = 10^(-8)){
  Sigma_R = matrix(sapply(rho_matrix, Cov_GenCauchy, theta_1 = theta_1, theta_2 = theta[1],theta_3=0.5,theta_4=2),
                   nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Two.Neg.Log.Likelihood.Use.Sigma_R(Y,Sigma_R,approx.option=approx.option,approx.par=approx.par)
}

Two.Neg.Log.RQ.given_pars <- function(theta,Y,rho_matrix,theta_1=1,approx.option = 2,approx.par = 10^(-8)){
  Sigma_R = matrix(sapply(rho_matrix, Cov_RQ, theta_1 = theta_1, theta_2 = theta),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Two.Neg.Log.Likelihood.Use.Sigma_R(Y,Sigma_R,approx.option=approx.option,approx.par=approx.par)
}

Two.Neg.Log.Cauchy.given_pars <- function(theta,Y,rho_matrix,theta_1=1,approx.option = 2,approx.par = 10^(-8)){
  Sigma_R = matrix(sapply(rho_matrix, Cov_Cauchy, theta_1 = theta_1, theta_2 = theta),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  Two.Neg.Log.Likelihood.Use.Sigma_R(Y,Sigma_R,approx.option=approx.option,approx.par=approx.par)
}


  




#' 
#' 
#' 
#' Calculate MLE estimation of Gaussian
#' 
#' @param par.init,length of 1, initial value of range value
#' @param Y
#' @param rho_matrix
#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------

MLE_est_Gaussian <- function(par.init,Y,rho_matrix){
  if(length(par.init)!=1){
    stop('len of par init of Gaussian should be 1')
  }
  tmp = optim(par=par.init,fn = Two.Neg.Log.Gaussian.given_pars,method='L-BFGS-B',
      lower = 10^(-6),
      Y = Y, rho_matrix = rho_matrix)
  
  Sigma_R = matrix(sapply(rho_matrix, Cov_Gaussian, theta_1 = 1, theta_2 = tmp$par),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  
  sigma2 = Cov.Zero.Approx(Y,Sigma_R,tol=10^(-8))
  return(list(par=c(sigma2,tmp$par),
              Sigma_R_est = Sigma_R,
              fvalue = tmp$value))
  
}



MLE_est_Exp <- function(par.init,Y,rho_matrix){
  if(length(par.init)!=1){
    stop('len of par init of Gaussian should be 1')
  }
  tmp = optim(par=par.init,fn = Two.Neg.Log.Exp.given_pars,method='L-BFGS-B',
      lower = 10^(-6),
      Y = Y, rho_matrix = rho_matrix)
  
  Sigma_R = matrix(sapply(rho_matrix, Cov_Exp, theta_1 = 1, theta_2 = tmp$par),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  
  sigma2 = Cov.Zero.Approx(Y,Sigma_R,tol=10^(-8))
  return(list(par=c(sigma2,tmp$par),
              Sigma_R_est = Sigma_R,
              fvalue = tmp$value))
  
}


MLE_est_Matern <- function(par.init,Y,rho_matrix,theta_3 = NULL){
  if(is.null(theta_3)){ #theta_3 is not given
  theta_3_seq = 1:10
  fvalue_seq = numeric(length(theta_3_seq))
  range_seq = numeric(length(theta_3_seq))
  
  for(i in 1:length(theta_3_seq)){
    tmp = optim(par=par.init,fn = Two.Neg.Log.Matern.given_pars,method='L-BFGS-B',
        lower = 10^(-6),theta_3=theta_3_seq[i],
        Y = Y, rho_matrix = rho_matrix)
    fvalue_seq[i] = tmp$value
    range_seq[i] = tmp$par
  }
  
  index=which.min(fvalue_seq)
  
  Sigma_R = matrix(sapply(rho_matrix, Cov_Matern, theta_1 = 1, theta_2 = range_seq[index],theta_3=theta_3_seq[index]),
                   nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  sigma2 = Cov.Zero.Approx(Y,Sigma_R,tol=10^(-8))
  
  
  return(list(par=c(sigma2,range_seq[index],theta_3_seq[index]),
       Sigma_R_est = Sigma_R,
       fvalue=fvalue_seq[index]))
    
  }else{#theta_3 is given
  tmp = optim(par=par.init,fn = Two.Neg.Log.Matern.given_pars,method='L-BFGS-B',
        lower = 10^(-6),theta_3=theta_3,
        Y = Y, rho_matrix = rho_matrix)
  
  Sigma_R = matrix(sapply(rho_matrix, Cov_Matern, theta_1 = 1, theta_2 = tmp$par,theta_3=theta_3),
                   nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  sigma2 = Cov.Zero.Approx(Y,Sigma_R,tol=10^(-8))
  
  return(list(par=c(sigma2,tmp$par,theta_3),
       Sigma_R_est = Sigma_R,
       fvalue=tmp$value))
  }
  
}



MLE_est_GenCauchy <- function(par.init,Y,rho_matrix){
  if(length(par.init)!=1){
    stop('len of par init of Gaussian should be 1')
  }
  tmp = optim(par=par.init,fn = Two.Neg.Log.GenCauchy.given_pars,method='L-BFGS-B',
      lower = 10^(-6),
      Y = Y, rho_matrix = rho_matrix)
  
  Sigma_R = matrix(sapply(rho_matrix, Cov_GenCauchy, theta_1 = 1, theta_2 = tmp$par,theta_3=0.5,theta_4=2),
                   nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  
  sigma2 = Cov.Zero.Approx(Y,Sigma_R,tol=10^(-8))
  return(list(par=c(sigma2,tmp$par,0.5,2),
              Sigma_R_est = Sigma_R,
              fvalue = tmp$value))
  
}



MLE_est_RQ <- function(par.init,Y,rho_matrix){
  if(length(par.init)!=1){
    stop('len of par init of Gaussian should be 1')
  }
  tmp = optim(par=par.init,fn = Two.Neg.Log.RQ.given_pars,method='L-BFGS-B',
      lower = 10^(-6),
      Y = Y, rho_matrix = rho_matrix)
  
  Sigma_R = matrix(sapply(rho_matrix, Cov_RQ, theta_1 = 1, theta_2 = tmp$par),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  
  sigma2 = Cov.Zero.Approx(Y,Sigma_R,tol=10^(-8))
  return(list(par=c(sigma2,tmp$par),
              Sigma_R_est = Sigma_R,
              fvalue = tmp$value))
  
}




MLE_est_Cauchy <- function(par.init,Y,rho_matrix){
  if(length(par.init)!=1){
    stop('len of par init of Gaussian should be 1')
  }
  tmp = optim(par=par.init,fn = Two.Neg.Log.Cauchy.given_pars,method='L-BFGS-B',
      lower = 10^(-6),
      Y = Y, rho_matrix = rho_matrix)
  
  Sigma_R = matrix(sapply(rho_matrix, Cov_Cauchy, theta_1 = 1, theta_2 = tmp$par),nrow=nrow(rho_matrix),ncol=ncol(rho_matrix))
  
  sigma2 = Cov.Zero.Approx(Y,Sigma_R,tol=10^(-8))
  return(list(par=c(sigma2,tmp$par),
              Sigma_R_est = Sigma_R,
              fvalue = tmp$value))
  
}
  

#' 
#' Calculate $\sqrt(\int_0^\Inf (f1-f2)^2dx)$
#' @param f1, function 1
#' @param f2, function 2
#' @param hmin
#' @param hmax
#' @param ..., parameters for f2
#' @return
#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------
Norm2_between_f <- function(f1,f2,hmin,hmax,...){
  diff_f <- function(x){
    (f1(x)-f2(x,...))^2
  }
  return(sqrt(integrate(diff_f,lower = hmin,upper = hmax)$value))
}

#' 
#' 
#' Calculate sup norm
## ----------------------------------------------------------------------------------------------------------------------------------------------------
SupNorm_between_f <- function(f1,f2,hmin,hmax,...){
   diff_f <- function(x){
    f1(x)-f2(x,...)
  }
  x_seq = seq(from=hmin,to=hmax,length.out = min(150,(hmax-hmin)*10))
  return(max(abs(diff_f(x_seq))))
}

#' 
#' 
## ----------------------------------------------------------------------------------------------------------------------------------------------------
Parametric_Estimation_Metric <- function(true_option,true_par.cov,est_option,Y,rho_matrix.orig,
                                         metric_upper_bound,metric_grid_space,metric_grid){
 
    if(est_option == 'Gaussian'){
      if(true_option == 'Gaussian'){
        Param_est = MLE_est_Gaussian(true_par.cov[2],Y,rho_matrix.orig)
      }else{
        Param_est = MLE_est_Gaussian(0.8,Y,rho_matrix.orig)
      }
    }else if(est_option == 'Matern'){
        Param_est = MLE_est_Matern(0.8,Y,rho_matrix.orig,theta_3=NULL)
    }else if(est_option == 'GenCauchy'){
        Param_est = MLE_est_GenCauchy(0.8,Y,rho_matrix.orig)
    }else if(est_option == 'Cauchy'){
      Param_est = MLE_est_Cauchy(0.8,Y,rho_matrix.orig)
    }else{
      stop('est_option is not known !')
    }
    
    #Define true cov and cor function
    f_cov_true <- function(x){
      return(Cov_function(x,option=true_option,par.cov=true_par.cov))
    }

    f_cor_true <- function(x){
      return(Cov_function(x,option=true_option,par.cov=c(1,true_par.cov[-1])))
    }
    
    
      
    #cor metric
    Param_est$norm2_cor_integration =  try(Norm2_between_f(f_cor_true,Cov_function,
                                                           hmin=0,hmax=metric_upper_bound,
                                                           option=est_option,par.cov = c(1,Param_est$par[-1]))/sqrt(metric_upper_bound),
                                           silent = TRUE)
    Param_est$norm2_cor_grid =  sqrt(sum(sapply(metric_grid,function(x){
      (f_cor_true(x)-Cov_function(x,option = est_option,par.cov=c(1,Param_est$par[-1])))^2
      }
      )) * metric_grid_space/metric_upper_bound)
    
    Param_est$supnorm_cor_grid =  max(sapply(metric_grid,function(x){
      abs(f_cor_true(x)-Cov_function(x,option = est_option,par.cov=c(1,Param_est$par[-1])))
      }
      ))
      
    #cov metric
    Param_est$norm2_cov_integration =  try(Norm2_between_f(f_cov_true,Cov_function,
                                                           hmin=0,hmax=metric_upper_bound,
                                                           option=est_option,par.cov = Param_est$par)/sqrt(metric_upper_bound),
                                              silent = TRUE)
    Param_est$norm2_cov_grid =  sqrt(sum(sapply(metric_grid,function(x){
      (f_cov_true(x)-Cov_function(x,option = est_option,par.cov=Param_est$par))^2})) * metric_grid_space/metric_upper_bound)
    Param_est$supnorm_cov_grid =  max(sapply(metric_grid,function(x){
      abs(f_cov_true(x)-Cov_function(x,option = est_option,par.cov=Param_est$par))
      }
      ))
    
    if('try-error' %in% class(Param_est$norm2_cor_integration)){
      Param_est$norm2_cor_integration = NA
    }
    if('try-error' %in% class(Param_est$norm2_cov_integration)){
      Param_est$norm2_cov_integration = NA
    }
    
    #variog metric
    f_variog_true <- function(x){
      return(Cov_function(0,option=true_option,par.cov = true_par.cov) - Cov_function(x,option=true_option,par.cov=true_par.cov))
    }
    f_variog_est <- function(x){
      return(Cov_function(0, option = est_option,par.cov=Param_est$par) - Cov_function(x,option = est_option, par.cov = Param_est$par))
    }
    
    Param_est$norm2_variog_integration =  try(Norm2_between_f(f_variog_true,f_variog_est,
                                                           hmin=0,hmax=metric_upper_bound)/sqrt(metric_upper_bound),
                                           silent = TRUE)
    Param_est$norm2_variog_grid =  sqrt(sum(sapply(metric_grid,function(x){
      (f_variog_true(x)-f_variog_est(x))^2
      }
      )) * metric_grid_space/metric_upper_bound)
    
    Param_est$supnorm_variog_grid =  max(sapply(metric_grid,function(x){
      abs(f_variog_true(x)-f_variog_est(x))
    }
    ))
    
    if('try-error' %in% class(Param_est$norm2_variog_integration)){
      Param_est$norm2_variog_integration = NA
    }
    
    
    return(list(option = est_option,par = Param_est$par,fvalue = Param_est$fvalue,
                                                     norm2_cov_grid = Param_est$norm2_cov_grid,
                                                     norm2_cov_integration = Param_est$norm2_cov_integration,
                                                     supnorm_cov_grid = Param_est$supnorm_cov_grid,
                                                     norm2_cor_grid = Param_est$norm2_cor_grid,
                                                     norm2_cor_integration = Param_est$norm2_cor_integration,
                                                     supnorm_cor_grid = Param_est$supnorm_cor_grid,
                                                     norm2_variog_grid = Param_est$norm2_variog_grid,
                                                     norm2_variog_integration = Param_est$norm2_variog_integration,
                                                     supnorm_variog_grid = Param_est$supnorm_variog_grid))
    
}

  

#' 
#' 
#' test
## ----eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------
## option = 'Exp'; par.cov=c(1,2)
## #option = 'Matern'; par.cov = c(1,1.5,1)
## #option = 'Gaussian'; par.cov = c(1,4)
## #option = 'GenCauchy'; par.cov=c(1,0.5,0.5,2)
## #option = 'RQ'; par.cov = c(1,2)
## l_max = 5
## #range [0,2*l_max]^d
## Data.Preparation.result = Data_Preparation_Fixed_range(option,par.cov,d=2,N=2,r=200,n=60,l_max=l_max,seed = 0)
## 
## Y = Data.Preparation.result$Output.list$Y_matrix[2,,]
## rho_matrix = Data.Preparation.result$Output.list$rho_matrix.orig
## true_Sigma_R = Data.Preparation.result$Hidden.list$True.Sigma_matrix
## 
## Gaussian_est = MLE_est_Gaussian(0.8,Y,rho_matrix)
## Gaussian_est$par
## Gaussian_est$fvalue
## 
## 
## Matern_est = MLE_est_Matern(0.8,Y,rho_matrix,theta_3=NULL)
## Matern_est$par
## Matern_est$fvalue
## 
## 
## 
## Exp_est = MLE_est_Exp(0.8,Y,rho_matrix)
## Exp_est$par
## Exp_est$fvalue
## 
## 
## 
## GenCauchy_est = MLE_est_GenCauchy(0.8,Y,rho_matrix)
## GenCauchy_est$par
## GenCauchy_est$fvalue
## 
## 
## 
## RQ_est = MLE_est_RQ(0.8,Y,rho_matrix)
## RQ_est$par
## RQ_est$fvalue
## 
## 
## Two.Neg.Log.Likelihood.Use.Sigma_R(Y,Data.Preparation.result$Hidden.list$True.Sigma_matrix,approx.option = 2,approx.par = 10^(-8))
## Cov.Zero.Approx(Y,Sigma_R = Data.Preparation.result$Hidden.list$True.Sigma_matrix,tol=approx.par)
## 
## 
## layout(matrix(1:2,1,2))
## rho_seq = seq(from=0,to=max(rho_matrix),length.out = 100)
## 
## plot(rho_matrix,Data.Preparation.result$Hidden.list$True.Sigma_matrix,type='p',pch=20,ylim=c(0,1),
##      main = paste0(option,' Cor'))
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='Exp',par.cov=c(1,Exp_est$par[2])),col='red',lwd=2)
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='Gaussian',par.cov=c(1,Gaussian_est$par[2])),col='green',lwd=2)
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='Matern',par.cov=c(1,Matern_est$par[2],Matern_est$par[3])),col='blue',lwd=2)
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='GenCauchy',par.cov=c(1,GenCauchy_est$par[2],GenCauchy_est$par[3],GenCauchy_est$par[4])),col='orange',lwd=2)
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='RQ',par.cov=c(1,RQ_est$par[2])),col='purple',lwd=2)
## 
## 
## legend('topright',legend=c('true','Exp_est','Gaussian_est','Matern_est','GenCauchy_est','RQ_est'),
##                   col = c('black','red','green','blue','orange','purple'),
##                   pch = c(20,rep(NA,5)),
##                   lty = c(NA,rep(1,5)))
## 
## 
## plot(rho_matrix,Data.Preparation.result$Hidden.list$True.Sigma_matrix,type='p',pch=20,ylim=c(0,1),
##      main = paste0(option,' Cov'))
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='Exp',par.cov=Exp_est$par),col='red',lwd=2)
## text(1,Exp_est$par[1],labels = round(Exp_est$par[1],digits=4))
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='Gaussian',par.cov=Gaussian_est$par),col='green',lwd=2)
## text(1,Gaussian_est$par[1],labels = round(Gaussian_est$par[1],digits=4))
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='Matern',par.cov=Matern_est$par),col='blue',lwd=2)
## text(1,Matern_est$par[1],labels = round(Matern_est$par[1],digits=4))
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='GenCauchy',par.cov=GenCauchy_est$par),col='orange',lwd=2)
## text(1,GenCauchy_est$par[1],labels = round(GenCauchy_est$par[1],digits=4))
## 
## 
## lines(rho_seq,sapply(rho_seq, Cov_function,option='RQ',par.cov=RQ_est$par),col='purple',lwd=2)
## text(1,RQ_est$par[1],labels = round(RQ_est$par[1],digits=4))
## 
## legend('topright',legend=c('true','Exp_est','Gaussian_est','Matern_est','GenCauchy_est','RQ_est'),
##                   col = c('black','red','green','blue','orange','purple'),
##                   pch = c(20,rep(NA,5)),
##                   lty = c(NA,rep(1,5)))
## 
## 
## 

#' 
#' 
## ----eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------
## knitr::purl('Parametric_MLE.Rmd',documentation = 2L)

