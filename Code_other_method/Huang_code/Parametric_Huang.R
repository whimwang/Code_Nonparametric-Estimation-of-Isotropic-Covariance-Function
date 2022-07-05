
## Matern ####
Matern_Parametric_est_Huang_model <- function(h,Variog_hat,weight){
  theta_3.vector <- seq(from=1,to=10,by=1)
  value.vector <- numeric(length(theta_3.vector))
  for(i in 1:length(theta_3.vector)){
    tmp <- try(optim(c(1,0.8),Matern_Variog_WeightedMSE,
                     h = h,
                     Variog_hat = Variog_hat,
                     weight = weight,
                     theta_3 = theta_3.vector[i],
                     lower=c(1e-6,1e-6),method="L-BFGS-B"),
               silent = TRUE)
    if('try-error' %in% class(tmp)){
      value.vector[i] = NA
    }else{
      value.vector[i] = tmp$value
    }
  }
  index = which.min(value.vector)
  optout_Matern <- try(optim(c(1,0.8),Matern_Variog_WeightedMSE,
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
Gaussian_Parametric_est_Huang_model <- function(h,Variog_hat,weight){
  optout_Gaussian <- try(optim(c(1,3),Gaussian_Variog_WeightedMSE,
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
Cauchy_Parametric_est_Huang_model <- function(h,Variog_hat,weight){
  optout_Cauchy <- try(optim(c(1,0.8),Cauchy_Variog_WeightedMSE,
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
GenCauchy_Parametric_est_Huang_model <- function(h,Variog_hat,weight){
  optout_GenCauchy <- try(optim(c(1,0.3),GenCauchy_Variog_WeightedMSE,
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




Parametric_est_Huang_model_metric <- function(true_option,true_par.cov,est_option,
                                             h,Variog_hat,weight,
                                             metric_upper_bound,metric_grid_space,metric_grid,
                                             v,L){
  
  if(est_option == 'Gaussian'){
    Param_est = Gaussian_Parametric_est_Huang_model(h=h,Variog_hat=Variog_hat,weight=weight)
  }else if(est_option == 'Matern'){
    Param_est = Matern_Parametric_est_Huang_model(h=h,Variog_hat=Variog_hat,weight=weight)
  }else if(est_option == 'GenCauchy'){
    Param_est = GenCauchy_Parametric_est_Huang_model(h=h,Variog_hat=Variog_hat,weight=weight)
  }else if(est_option == 'Cauchy'){
    Param_est = Cauchy_Parametric_est_Huang_model(h=h,Variog_hat=Variog_hat,weight=weight)
  }else{
    stop('est_option is not known !')
  }
  
  if(!is.null(Param_est$par)){ #success
    
    #Define true cov and cor function
    f_cov_true <- function(x){
      return(Cov_function(x,option=true_option,par.cov=true_par.cov))
    }
    
    f_cor_true <- function(x){
      return(Cov_function(x,option=true_option,par.cov=c(1,true_par.cov[-1])))
    }
    
    f_variog_true <- function(x){
      return(par.cov[1] - Cov_function(x,option = true_option, par.cov=true_par.cov))
    }
    
    #est
    f_cov_est <- function(x){
      return(Cov_function(x,option=est_option,par.cov=Param_est$par))
    }
    f_cor_est <- function(x){
      return(Cov_function(x,option=est_option,par.cov=c(1,Param_est$par[-1])))
    }
    f_variog_est <- function(x){
      return(Param_est$par[1] - Cov_function(x,option = est_option, par.cov=Param_est$par))
    }
    
    
    if(!is.null(Spectral_function(1,option=true_option,par.cov=true_par.cov)) &
      !is.null(Spectral_function(1,option=est_option,par.cov=Param_est$par))){ #spectral of est valid #true spectral exists
         w_seq = seq(from=1,to=L,by=1) * v/L
         spectral_true = sapply(w_seq, Spectral_function,option = true_option,par.cov=true_par.cov)
         spectral_est = sapply(w_seq, Spectral_function,option = est_option, par.cov = Param_est$par)
         Param_est$norm2_spectral_grid = sqrt(sum((spectral_est-spectral_true)^2) *v/L)
     
    }else{#
      Param_est$norm2_spectral_grid = NA
    }
    
    
    
    #cor metric
    Param_est$norm2_cor_integration =  try(Norm2_between_f(f_cor_true,f_cor_est,
                                                           hmin=0,hmax=metric_upper_bound)/sqrt(metric_upper_bound),
                                           silent = TRUE)
    Param_est$norm2_cor_grid =  sqrt(sum(sapply(metric_grid,function(x){
      (f_cor_true(x)-f_cor_est(x))^2
    }
    )) * metric_grid_space/metric_upper_bound)
    
    Param_est$supnorm_cor_grid =  max(sapply(metric_grid,function(x){
      abs(f_cor_true(x)-f_cor_est(x))
    }
    ))
    
    #cov metric
    Param_est$norm2_cov_integration =  try(Norm2_between_f(f_cov_true,f_cov_est,
                                                           hmin=0,hmax=metric_upper_bound)/sqrt(metric_upper_bound),
                                           silent = TRUE)
    Param_est$norm2_cov_grid =  sqrt(sum(sapply(metric_grid,function(x){
      (f_cov_true(x)-f_cov_est(x))^2
    }
    )) * metric_grid_space/metric_upper_bound)
    
    Param_est$supnorm_cov_grid =  max(sapply(metric_grid,function(x){
      abs(f_cov_true(x)-f_cov_est(x))
    }
    ))
    
    if('try-error' %in% class(Param_est$norm2_cor_integration)){
      Param_est$norm2_cor_integration = NA
    }
    if('try-error' %in% class(Param_est$norm2_cov_integration)){
      Param_est$norm2_cov_integration = NA
    }
    
    #variog metric
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
    
  }
  
  Param_est$option = est_option
  return(Param_est)
  
}


