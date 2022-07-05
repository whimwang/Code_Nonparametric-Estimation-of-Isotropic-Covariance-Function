#
code.path = '/Users/wangyiming/Documents/2019 Research Daily/research_code/Code_2020/Updated0620/'
#Server
#code.path = "/mnt/home/ywang225/Cov_research/Updated0919/"

source(paste0(code.path,'Code0620/',"Code0912-0919/Part-I.Spectral.Cov.R"))
source(paste0(code.path,'Code0620/',"Code0912-0919/Part-II.Data_Generation.R"))
source(paste0(code.path,'Code0620/',"Code1004-1010/Support.R"))
source(paste0(code.path,'Code0620/',"Code1018-1025/Data_Preparation.R"))

source(paste0(code.path,'Code_other_method/','Huang_code/Huang_Variog.R'))
source(paste0(code.path,'Code_other_method/','Bspline_Choi/Bspline_Choi.R'))
source(paste0(code.path,'Code_other_method/','Bspline_Choi/NP_Choi_Estimation.R'))
source(paste0(code.path,'Code_other_method/','Cov_Variog_Weighted_MSE.R'))

source(paste0(code.path,'Code_other_method/','Bspline_Choi/Parametric_Choi.R')) #norm2 function
source(paste0(code.path,'Code_other_method/','Huang_code/Parametric_Huang.R'))


dependency_packages = c('MASS','geoR','quadprog')

library('quadprog')
library('MASS')
library('geoR')
library('foreach')
library('doParallel')

## convert a list to a vector: c(name,value,name,value,...)
list_to_arg_vector <- function(arg_list){
  as.vector(sapply(1:(length(arg_list)), 
                   function(i){c(names(arg_list[i]),as.character(arg_list[i]))}))
}


#getwd()

cores <- detectCores(logical=T)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)


######################################################################################
#####################    Specify Arguments   #################################################

do_plot = TRUE; 


#0916:Consider Gau,Cauchy_small,LinearTransMat,TransMat,GenCauchy
N = 2
option = 'Matern'; par.cov=c(1,1.25,1);metric_epsilon = 0.001; ## 10.28812
#option = 'Cauchy'; par.cov = c(1,0.8); metric_epsilon = 0.05; ## 0.01:79.996 0.03:26.65466  0.04:19.98 0.05:15:97999 when dist=23.3026,cov value is 0.034
#option = 'Gaussian'; par.cov = c(1,3); metric_epsilon = 0.001; ## 7.88
#option = 'LinearMatern'; par.cov=c(1,1.25,1,3/2/sqrt(2),2); metric_epsilon = 0.001;m_seq = m_seq[m_seq < 150] ## 10.515
#option = 'GenCauchy'; par.cov=c(1,0.3,0.5,2); metric_epsilon=0.15 # is Cov_function(max(observed_dist))
##GenCauchy: 0.15: 13.32996 0.125:19.19766 0.1:29.9985 0.1134602:23.3026




## specify coordinates
whether_irregular = TRUE
### regular
if(!whether_irregular){
  coord_arg = list(whether_irregular=whether_irregular,r=1,d=2,space=0.5,l_min=0,l_max=20)
}else{
  ### irregular
  coord_arg = list(whether_irregular=whether_irregular,r=200,d=2,n=60,l_min=0,l_max=20)
}


## specify fitting args
HuangFitted_arg = list(v=10,L=200)
if(whether_irregular){ #when irregular,use all
  HuangFitted_arg$fit_ratio=1
  #1:23.3026; 0.9:16.23687  0.85:15.07886 0.8:14.18993
}else{ #when regular
  HuangFitted_arg$fit_ratio=1
}
#Since irregular, r=200 replicates should ensure a relative stable estimation of variogram.


#
#result_dir = paste0(getwd(),'/')
#server
result_dir = '/mnt/home/ywang225/Cov_research/Result_0919/Huang/'
filename = paste(c(option,list_to_arg_vector(HuangFitted_arg),list_to_arg_vector(coord_arg),format(Sys.time(), "%b%d_%X")),collapse = '_')
result_path = paste0(result_dir,filename)
dir.create(result_path)




######################################################################################
#####################    Coord Generation    #################################################

Data.Preparation.result = Data_Preparation_Fixed_range(option,par.cov,d=coord_arg$d,N,
                                                       r=coord_arg$r,n=coord_arg$n,
                                                       l_max=coord_arg$l_max,seed = 0)

save(Data.Preparation.result, file=paste0(result_path,"/Data.Preparation.result.RData"))


##### for metric grid setting #####
metric_upper_bound = Support.function(par.cov = c(1,par.cov[-1]),option = option, epsilon = metric_epsilon)
metric_grid_count = 2000
metric_grid = seq(from=0,to=metric_upper_bound,length.out = metric_grid_count)
metric_grid_space = metric_upper_bound/metric_grid_count


start_time = Sys.time()


all.result <- foreach(k = seq(from=1,to=N,by=1),.packages=dependency_packages,.combine = 'rbind') %dopar%
  {
    tryCatch({
      print(k)
    
    
    print(paste('Simulation',k,'starts!'))
    start_iteration = Sys.time()
    print(Sys.time())
    
    if(do_plot){
      pdf(file=paste0(result_path,'/',k,'_',format(Sys.time(), "%b%d_%X"),'_plot.pdf'))
    }
    
    ####################################################################################
    ################## Cov Generation and Observation Generation #####################
    
    #generate random param and corresponding cov function
    coord_matrix = Data.Preparation.result$Output.list$coord_matrix
    Y = Data.Preparation.result$Output.list$Y_matrix[k,,]
    
    ## empirical stat
    Huang_result_orig = Huang_Empirical_Variogram(Y = Y,coord_matrix = coord_matrix)
    Huang_result_orig$empirical_varioghat = Huang_result_orig$sum_varioghat/Huang_result_orig$n_h
    Huang_result_orig$weight = Huang_result_orig$n_h
    
    ## specify fitting args
    v_range <- HuangFitted_arg$v
    L_grid <- HuangFitted_arg$L
    if(!is.null(Spectral_function(1,option=option,par.cov = par.cov))){ #spectral exist
      ftrue <- sapply(seq(from=1,to=L_grid,by=1) * v_range/L_grid, Spectral_function,option=option,par.cov=par.cov)
    }else{ 
      ftrue <- NULL
    }
   
    
    n_dist = length(Huang_result_orig$dist)
    n_fit = as.integer(n_dist* HuangFitted_arg$fit_ratio)
    
    fitting_dist = Huang_result_orig$dist[1:n_fit]
    fitting_variog = Huang_result_orig$empirical_varioghat[1:n_fit]
    
    
    ### QP paramters
    X <- Xnew <- Get_Bmatrix(dist = fitting_dist,v=v_range,L=L_grid,nugget = FALSE) #B in paper
    Dc <- Dnew <- Get_PenalizeKmatrix(v=v_range,L=L_grid,nugget = FALSE)
    weight = Huang_result_orig$n_h[1:n_fit]
    WW = diag(weight)
    XnewWWXnew <- t(Xnew)%*%WW%*%Xnew
    dvec = t(Xnew) %*% WW %*% fitting_variog
    Amat = diag(L_grid)
    bvec = rep(0,L_grid)
    
    ## choose lamdba
    ngcv = 20
    lossf <- lossgamma <- gcvnew <- tr1 <- tr2 <- numeric(ngcv)
    lambda_seq = sapply((1:ngcv)-1, function(x){10^(6*(x-1)/(ngcv-1))})
    
    #for(igcv in 1:ngcv){
    for(igcv in 1:ngcv){
      
      print(paste0('cross-validation: ',igcv))
      print(paste0('lambda value: ',lambda_seq[igcv]))
      lambda = lambda_seq[igcv]
      
      ### solve QP without nugget
      Dmat <- XnewWWXnew + lambda*Dnew
      QPout <- solve.QP(Dmat,dvec,Amat,bvec=bvec)
      fest <- QPout$solution
      variog_hat <- Xnew%*%QPout$solution #in paper: Bv or say l_ij g_w
      
      if(do_plot){
        if(igcv == 1){
          plot(Huang_result_orig$dist, Huang_result_orig$empirical_varioghat,
               xlab='dist',ylab='variog value',
               main=paste0('Variog of ',option))
          #true
          #points(Huang_result_orig$dist, sapply(Huang_result_orig$dist,GenCov_list$variog_function),col='red',pch=20)
          lines(Huang_result_orig$dist, sapply(Huang_result_orig$dist,function(x){
            2*Cov_function(0,option=option,par.cov=par.cov) - 2*Cov_function(x,option=option,par.cov=par.cov)}),col='black',lwd=5)
        
          }
        lines(fitting_dist,variog_hat,col='red',lwd=5)
        legend('topright',legend = c('empirical variog','true variog','est variog'),
               col = c('black','black','red'),
               pch = c(1,NA,NA),
               lty=c(NA,1,1),bty='n')
      }
      
      ## Find constraints which are non-active constraints
      index <- (fest < 1e-6) 
      remove <- which(index) #index non-active
      keep <- which(!index) #index active
      Xprime <- X[,keep] #B_tilde
      Dcprime <- Dc[keep,keep] 
      
      XprimeW = t(Xprime) %*% WW
      ## A(lambda)
      BprimeWBprime = t(Xprime)%*%WW%*%Xprime
      A_lambda = Xprime %*% solve(BprimeWBprime + lambda * Dcprime) %*% XprimeW
      
      ## A(0)
      eig <- eigen(t(Xprime)%*%WW%*%Xprime)
      eva <- eig$values
      eve <- eig$vectors
      inv <- gencov(eva,eve,-1) ## Inv of (BprimeWBprime)
      A_0 = Xprime %*% inv %*% XprimeW
      p_value <- sum(diag(WW%*%A_0))  #p_value in paper
      
      gcvnew[igcv] <- sum(weight*(variog_hat-fitting_variog)^2)/((1-sum(diag(WW%*%A_lambda))/p_value)^2)
      if(!is.null(ftrue)){ #ISE of spectral
        lossf[igcv] <- sqrt(sum((fest-ftrue)^2) * HuangFitted_arg$v/HuangFitted_arg$L) #est spectral vs true spectral
      }
      #ISE of variog
      f_variog_est <- function(dist){
        Empirical_Variog_est(dist,f_est=QPout$solution,v=HuangFitted_arg$v,L=HuangFitted_arg$L)}
      
      lossgamma[igcv] <- Norm2_between_f(f1=f_variog_est,f2=Variog_function,
                      hmin=0,hmax=metric_upper_bound,
                      option=option,par.cov=par.cov)
    }
    
    ## save lambda
    lambda_list = list(lambda = lambda_seq,ISE_f = lossf,ISE_gamma = lossgamma,gcvnew = gcvnew)
    
    ## save QP parameters
    QP_list = list(X=X,Xnew=Xnew,Dc=Dc,Dnew=Dnew,weight=weight,WW=WW,dvec=dvec,Amat=Amat,bvec=bvec)
    
    ########################################
    ####### Metric of Best solution ##########
    if(option == 'Wave'){
      best_igcv = 1
    }else{
      best_igcv = which(gcvnew == min(gcvnew))
    }
    best_lambda = lambda_seq[best_igcv]
    
    ### solve QP first time
    Dmat <- t(Xnew)%*%WW%*%Xnew + best_lambda*Dnew
    QPout <- solve.QP(Dmat,dvec,Amat,bvec=bvec)
    Huang_fest <- QPout$solution
    gammahat <- Xnew%*%QPout$solution #in paper: Bv or say l_ij g_w
    
    ################################ 
    ######## Define est function  ##### 
   
    f_variog_est <- function(dist){
      0.5*Empirical_Variog_est(dist,f_est=Huang_fest,v=HuangFitted_arg$v,L=HuangFitted_arg$L)
    }
    
    #Version 1
    #best_C_zero = f_variog_est(metric_upper_bound)/2
    #Version 2
    best_C_zero = f_variog_est(max(Data.Preparation.result$Output.list$rho_matrix.orig))

    f_cov_est <- function(dist){
      best_C_zero - f_variog_est(dist)
    }
    
    f_cor_est <- function(dist){
      f_cov_est(dist)/best_C_zero
    }
   
    f_variog_true <- function(dist){
      0.5*Variog_function(dist,option=option,par.cov=par.cov)
    }
    
    f_cov_true <- function(dist){
      Cov_function(dist,option=option,par.cov=par.cov)
    }
    
    f_cor_true <- function(dist){
      Cov_function(dist,option=option,par.cov=c(1,par.cov[-1]))
    }
    
    ### spectral metric: discretized norm-2
    best_ISE_spectral = lossf[best_igcv]
    
    
    ##### variog metric ######
    best_norm2_variog_integration =  try(Norm2_between_f(f1=f_variog_est,f2=f_variog_true,
                                                                       hmin=0,hmax=metric_upper_bound)/sqrt(metric_upper_bound),
                                                       silent = FALSE)
    if('try-error' %in% class(best_norm2_variog_integration)){
      best_norm2_variog_integration = NA
    }
   
    best_norm2_variog_grid = sqrt(sum((f_variog_true(metric_grid) - f_variog_est(metric_grid))^2 * metric_grid_space)/metric_upper_bound) 
    best_supnorm_variog_grid = max(abs(f_variog_true(metric_grid) - f_variog_est(metric_grid))) 
    
    
    ##### cov metric
    best_norm2_cov_integration =  try(Norm2_between_f(f1=f_cov_est,f2=f_cov_true,
                                                         hmin=0,hmax=metric_upper_bound)/sqrt(metric_upper_bound),
                                         silent = FALSE)
    if('try-error' %in% class(best_norm2_cov_integration)){
      best_norm2_cov_integration = NA
    }
    
    best_norm2_cov_grid = sqrt(sum((f_cov_true(metric_grid) - f_cov_est(metric_grid))^2 * metric_grid_space)/metric_upper_bound) 
    best_supnorm_cov_grid = max(abs(f_cov_true(metric_grid) - f_cov_est(metric_grid))) 
    
    ##### cor metric
    best_norm2_cor_integration =  try(Norm2_between_f(f1=f_cor_est,f2=f_cor_true,
                                                      hmin=0,hmax=metric_upper_bound)/sqrt(metric_upper_bound),
                                      silent = FALSE)
    if('try-error' %in% class(best_norm2_cor_integration)){
      best_norm2_cor_integration = NA
    }
    
    best_norm2_cor_grid = sqrt(sum((f_cor_true(metric_grid) - f_cor_est(metric_grid))^2 * metric_grid_space)/metric_upper_bound) 
    best_supnorm_cor_grid = max(abs(f_cor_true(metric_grid) - f_cor_est(metric_grid))) 
    
    print(paste0('cross-validation: ',best_igcv,'  best lambda value: ', best_lambda))
    
    
    ######################################################################################
    #####################    Parametric Method    #################################################
    
    Matern_est = Parametric_est_Huang_model_metric(true_option=option,true_par.cov=par.cov,est_option='Matern',
                                                  h=Huang_result_orig$dist,
                                                  Variog_hat=Huang_result_orig$empirical_varioghat,
                                                  weight=Huang_result_orig$weight,
                                                  metric_upper_bound=metric_upper_bound,metric_grid_space=metric_grid_space,metric_grid = metric_grid,
                                                  v = HuangFitted_arg$v,
                                                  L = HuangFitted_arg$L)
   
    Cauchy_est = Parametric_est_Huang_model_metric(true_option=option,true_par.cov=par.cov,est_option='Cauchy',
                                                   h=Huang_result_orig$dist,
                                                   Variog_hat=Huang_result_orig$empirical_varioghat,
                                                   weight=Huang_result_orig$weight,
                                                   metric_upper_bound=metric_upper_bound,metric_grid_space=metric_grid_space,metric_grid = metric_grid,
                                                   v = HuangFitted_arg$v,
                                                   L = HuangFitted_arg$L)
    
    Gaussian_est = Parametric_est_Huang_model_metric(true_option=option,true_par.cov=par.cov,est_option='Gaussian',
                                                   h=Huang_result_orig$dist,
                                                   Variog_hat=Huang_result_orig$empirical_varioghat,
                                                   weight=Huang_result_orig$weight,
                                                   metric_upper_bound=metric_upper_bound,metric_grid_space=metric_grid_space,metric_grid = metric_grid,
                                                   v = HuangFitted_arg$v,
                                                   L = HuangFitted_arg$L)
    
    GenCauchy_est = Parametric_est_Huang_model_metric(true_option=option,true_par.cov=par.cov,est_option='GenCauchy',
                                                     h=Huang_result_orig$dist,
                                                     Variog_hat=Huang_result_orig$empirical_varioghat,
                                                     weight=Huang_result_orig$weight,
                                                     metric_upper_bound=metric_upper_bound,metric_grid_space=metric_grid_space,metric_grid = metric_grid,
                                                     v = HuangFitted_arg$v,
                                                     L = HuangFitted_arg$L)
    
    Parametric_result_list <- list(Matern_est = Matern_est,
                                   Cauchy_est = Cauchy_est,
                                   Gaussian_est = Gaussian_est,
                                   GenCauchy_est = GenCauchy_est)
    

    ######################################################################################
    #####################    Record Results     #################################################
    
    if(do_plot){
      
      ## First plot: variog
      d_seq = seq(from=0,to=max(Huang_result_orig$dist)*1.2,length.out = 100)
      ## empirical variog_hat
      plot(Huang_result_orig$dist,Huang_result_orig$empirical_varioghat,type='p',pch=1,
           xlab='dist',ylab='variog value',
           main=paste('Huang method of ',option),
           ylim=range(Huang_result_orig$empirical_varioghat))
      ## true variogram
      lines(d_seq,sapply(d_seq, Variog_function,option=option,par.cov=par.cov),col='black',lty='solid',lwd=5)
      ## Nonparametric
      lines(d_seq,2*sapply(d_seq, f_variog_est),col='red',lwd=5,lty='dashed')
      all_color = c('red','blue','green','purple','orange','pink','yellow')
      for(i in 1:length(Parametric_result_list)){
        Current_list = Parametric_result_list[[i]]
        if(!is.null(Current_list$par)){
          lines(d_seq, sapply(d_seq,Variog_function,option = Current_list$option,par.cov=Current_list$par),col=all_color[i+1],lty='dashed',lwd=2)
        }
      }
      
      legend('topright',legend=c('empirical variog','true variog','Huang NP','Matern','Cauchy','Gaussian','GenCauchy'),
             col = c('black','black',all_color),
             pch = c(20,NA,rep(NA,length(all_color))),
             lty = c(NA,'solid',rep('dashed',length(all_color))),bty='n')
      
      ## Second plot: cov
      range1 = range(sapply(d_seq, Cov_function,option=option,par.cov=par.cov))
      range2 = range(Huang_result_orig$empirical_varioghat)/2
      plot(d_seq, sapply(d_seq, Cov_function,option=option,par.cov=par.cov),type='l',
           xlab='dist',ylab='cov value',lwd=5,
           main=paste('Huang method of ',option),
           ylim=range(c(range1,range2)))
      ## Nonparametric
      lines(d_seq,sapply(d_seq, f_cov_est),col='red',lwd=5,lty='dashed')
      for(i in 1:length(Parametric_result_list)){
        Current_list = Parametric_result_list[[i]]
        if(!is.null(Current_list$par)){
          lines(d_seq, sapply(d_seq,Cov_function,option = Current_list$option,par.cov=Current_list$par),col=all_color[i+1],lty='dashed',lwd=2)
        }
      }
      legend('topright',legend=c('empirical variog','true variog','Huang NP','Matern','Cauchy','Gaussian','GenCauchy'),
             col = c('black','black',all_color),
             pch = c(20,NA,rep(NA,length(all_color))),
             lty = c(NA,'solid',rep('dashed',length(all_color))),bty='n')
      
      ## Third plot: cor
      plot(d_seq, sapply(d_seq, Cov_function,option=option,par.cov=c(1,par.cov[-1])),type='l',
           xlab='dist',ylab='cor value',lwd=5,
           main=paste('Huang method of ',option),
           ylim=c(0,1.2))
      ## Nonparametric
      lines(d_seq,sapply(d_seq, f_cor_est),col='red',lwd=5,lty='dashed')
      for(i in 1:length(Parametric_result_list)){
        Current_list = Parametric_result_list[[i]]
        if(!is.null(Current_list$par)){
          lines(d_seq, sapply(d_seq,Cov_function,option = Current_list$option,par.cov=c(1,Current_list$par[-1])),
                col=all_color[i+1],lty='dashed',lwd=2)
        }
      }
      legend('topright',legend=c('empirical variog','true variog','Huang NP','Matern','Cauchy','Gaussian','GenCauchy'),
             col = c('black','black',all_color),
             pch = c(20,NA,rep(NA,length(all_color))),
             lty = c(NA,'solid',rep('dashed',length(all_color))),bty='n')
      
      ## Fourth plot: spectral
      w_seq = (1:L_grid) * v_range/L_grid
      #Huang
      plot(w_seq,Huang_fest,col='red',main=paste0('Huang method of ',option),
           xlab='w',ylab='spectral',pch=20,
           ylim=c(min(c(Huang_fest,ftrue))-0.1,max(c(Huang_fest,ftrue))+0.1))
      lines(w_seq,Huang_fest,col='red',lty='dashed',lwd=5)
      lines(w_seq,ftrue,col='black',lwd=5)
      for(i in 1:length(Parametric_result_list)){
        Current_list = Parametric_result_list[[i]]
        if(!is.null(Current_list$par) & (!is.null(Spectral_function(1,option=Current_list$option,par.cov = Current_list$par)))){
          lines(w_seq, sapply(w_seq,Spectral_function,option=Current_list$option,par.cov = Current_list$par),
                col=all_color[i+1],lty='dashed',lwd=2)
        }
      }
      legend('topright',legend=c('empirical variog','true variog','Huang NP','Matern','Cauchy','Gaussian','GenCauchy'),
             col = c('black','black',all_color),
             pch = c(20,NA,rep(NA,length(all_color))),
             lty = c(NA,'solid',rep('dashed',length(all_color))),bty='n')
      
      
      dev.off()
    }
    
    #####################################################################################
    #####################    Save Results    ############################################
    optout_Huang = list(par=QPout$solution,
                        C_zero = best_C_zero,
                        norm2_cov_grid = best_norm2_cov_grid,
                        norm2_cov_integration = best_norm2_cov_integration,
                        supnorm_cov_grid = best_supnorm_cov_grid,
                        norm2_cor_grid = best_norm2_cor_grid,
                        norm2_cor_integration = best_norm2_cor_integration,
                        supnorm_cor_grid = best_supnorm_cor_grid,
                        norm2_variog_grid = best_norm2_variog_grid,
                        norm2_variog_integration = best_norm2_variog_integration,
                        supnorm_variog_grid = best_supnorm_variog_grid,
                        norm2_spectral_grid = best_ISE_spectral,
                        lambda = best_lambda)
    
    Huang_result_list <-  list(optout_Huang = optout_Huang,
                               Huang_result_orig = Huang_result_orig,
                               lambda_list = lambda_list,
                               QP_list = QP_list,
                               args = list(HuangFitted_arg = HuangFitted_arg, coord_arg = coord_arg))
    
    save(Huang_result_list,Parametric_result_list,
         file=paste0(result_path,'/',k,'_',format(Sys.time(), "%b%d_%X"),'.RData'))
    
    one.sample.list = list()
    one.sample.list[[1]] = list(Huang_result_list = Huang_result_list,
                                Parametric_result_list = Parametric_result_list)
    
    return(one.sample.list)
    },
    error=function(e){
      content = conditionMessage(e)
      print(paste0(k,' Error: ',content))
      write(content,file=paste0(result_path,'/',k,'_error_',format(Sys.time(), "%b%d_%X"),'.txt'))
    }
    )
}



end_time <- Sys.time()
print(end_time - start_time)




metric_colnames = c('C_zero','norm2_cov_grid','norm2_cov_integration','supnorm_cov_grid',
                    'norm2_cor_grid','norm2_cor_integration','supnorm_cor_grid',
                    'norm2_variog_grid','norm2_variog_integration','supnorm_variog_grid','norm2_spectral_grid')

Nonparametric_metric_matrix = t(sapply(all.result, function(x){
  optout_Huang = x$Huang_result_list$optout_Huang
  if(is.null(optout_Huang$par)){
    return(rep(NA,11))
  }else{
    return(c(optout_Huang$C_zero,
             optout_Huang$norm2_cov_grid,optout_Huang$norm2_cov_integration,optout_Huang$supnorm_cov_grid,
             optout_Huang$norm2_cor_grid,optout_Huang$norm2_cor_integration,optout_Huang$supnorm_cor_grid,
             optout_Huang$norm2_variog_grid,optout_Huang$norm2_variog_integration,optout_Huang$supnorm_variog_grid,
             optout_Huang$norm2_spectral_grid))
  }
}))
colnames(Nonparametric_metric_matrix) = metric_colnames
nonparametric_metric.file = paste0(result_path,'/Huang_metric_all','.txt')
write.table(Nonparametric_metric_matrix, file = nonparametric_metric.file,
            row.names = paste0('result.',1:length(all.result)),col.names = TRUE)





## save Parametric estimate value
parametric_order = c('Matern_est','Cauchy_est','Gaussian_est','GenCauchy_est')

param_metric_matrix = t(sapply(all.result, function(x){
  if(!is.null(x$Huang_result_list$optout_Huang$par)){ #sucess
    c(x$Parametric_result_list$Matern_est$par,
      x$Parametric_result_list$Cauchy_est$par,
      x$Parametric_result_list$Gaussian_est$par,
      x$Parametric_result_list$GenCauchy_est$par,
      sapply(x$Parametric_result_list, function(opt_est){opt_est$norm2_cov_grid}),
      sapply(x$Parametric_result_list, function(opt_est){opt_est$norm2_cov_integration}),
      sapply(x$Parametric_result_list, function(opt_est){opt_est$supnorm_cov_grid}),
      sapply(x$Parametric_result_list, function(opt_est){opt_est$norm2_cor_grid}),
      sapply(x$Parametric_result_list, function(opt_est){opt_est$norm2_cor_integration}),
      sapply(x$Parametric_result_list, function(opt_est){opt_est$supnorm_cor_grid}),
      sapply(x$Parametric_result_list, function(opt_est){opt_est$norm2_variog_grid}),
      sapply(x$Parametric_result_list, function(opt_est){opt_est$norm2_variog_integration}),
      sapply(x$Parametric_result_list, function(opt_est){opt_est$supnorm_variog_grid}),
      sapply(x$Parametric_result_list, function(opt_est){opt_est$norm2_spectral_grid})
      )
  }else{
    return(rep(NA,47))
  }
}))

colnames(param_metric_matrix) = c(paste0('Matern.',1:3),paste0('Cauchy.',1:2),paste0('Gaussian.',1:2),paste0('GenCauchy.',1:4),
                                  paste0(parametric_order,'.norm2_cov_grid'),
                                  paste0(parametric_order,'.norm2_cov_integration'),
                                  paste0(parametric_order,'.supnorm_cov_grid'),
                                  paste0(parametric_order,'.norm2_cor_grid'),
                                  paste0(parametric_order,'.norm2_cor_integration'),
                                  paste0(parametric_order,'.supnorm_cor_grid'),
                                  paste0(parametric_order,'.norm2_variog_grid'),
                                  paste0(parametric_order,'.norm2_variog_integration'),
                                  paste0(parametric_order,'.supnorm_variog_grid'),
                                  paste0(parametric_order,'.norm2_spectral_grid'))
                                  
mean_metric= apply(param_metric_matrix,2,function(x){mean(x,na.rm=TRUE)})
sd_metric= apply(param_metric_matrix,2,function(x){sd(x,na.rm=TRUE)})

param_metric_matrix = rbind(param_metric_matrix,mean_metric)
param_metric_matrix = rbind(param_metric_matrix,sd_metric)

print(option)
print('mean')
print(mean_metric)
print('sd')
print(sd_metric)


write.table(param_metric_matrix,file=paste0(result_path,'/Parametric_metric_result.txt'),row.names = FALSE)


save(all.result,file=paste0(result_path,'/all.result.Total_',N,'.RData'))
save.image(file=paste0(result_path,'/Data_',N,'.RData'))


stopCluster(cl)

