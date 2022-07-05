## Compare BsplineNP vs Parametric 


code.path = '/Users/wangyiming/Documents/2019 Research Daily/research_code/Code_2020/Updated0620/'
#Server
#code.path = "/mnt/home/ywang225/Cov_research/Updated0923/"


source(paste0(code.path,'Code0620/',"Code0912-0919/Part-I.Spectral.Cov.R"))
source(paste0(code.path,'Code0620/',"Code0912-0919/Part-II.Data_Generation.R"))
source(paste0(code.path,'Code0620/',"Code1004-1010/Support.R"))
source(paste0(code.path,'Code0620/',"Code1018-1025/Data_Preparation.R"))

source(paste0(code.path,'Code_other_method/','Huang_code/Huang_Variog.R'))
source(paste0(code.path,'Code_other_method/','Bspline_Choi/Parametric_Choi.R'))
source(paste0(code.path,'Code_other_method/','Bspline_Choi/Bspline_Choi.R'))
source(paste0(code.path,'Code_other_method/','Bspline_Choi/NP_Choi_Estimation.R'))

source(paste0(code.path,'Code_other_method/','Cov_Variog_Weighted_MSE.R'))


dependency_packages = c('MASS','geoR')

library('MASS')
library('geoR')
library('foreach')
library('doParallel')

## convert a list to a vector: c(name,value,name,value,...)
list_to_arg_vector <- function(arg_list){
  as.vector(sapply(1:(length(arg_list)), function(i){c(names(arg_list[i]),as.character(arg_list[i]))}))
}



#getwd()
cores <- detectCores(logical=T)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)



######################################################################################
#####################    Specify Arguments   #################################################
do_plot = TRUE

## specify coordinates


#0916:Consider Gau,Cauchy_small,LinearTransMat,TransMat,GenCauchy
N = 100
#option = 'Matern'; par.cov=c(1,1.25,1);metric_epsilon = 0.001;
#option = 'Cauchy'; par.cov = c(1,0.8); metric_epsilon = 0.05;
#option = 'Gaussian'; par.cov = c(1,3); metric_epsilon = 0.001;
#option = 'LinearMatern'; par.cov=c(1,1.25,1,3/2/sqrt(2),2); metric_epsilon = 0.001;
option = 'GenCauchy'; par.cov=c(1,0.3,0.5,2); metric_epsilon = 0.15



##### for metric grid setting #####
metric_upper_bound = Support.function(par.cov = c(1,par.cov[-1]),option = option, epsilon = metric_epsilon)
metric_grid_count = 2000
metric_grid = seq(from=0,to=metric_upper_bound,length.out = metric_grid_count)
metric_grid_space = metric_upper_bound/metric_grid_count



## arguments for Choi Nonpametric Fitting
#m:3,5,7 
#p:2,3,5
ChoiFitted_arg = list(m =5,p=3)



whether_irregular = TRUE
### regular
if(!whether_irregular){
  coord_arg = list(whether_irregular=whether_irregular,r=1,d=2,space=0.5,l_min=0,l_max=20)
}else{
  ### irregular
  coord_arg = list(whether_irregular=whether_irregular,r=200,d=2,n=60,l_min=0,l_max=20)
}



#result_dir = paste0(getwd(),"/")
#server
result_dir = '/mnt/home/ywang225/Cov_research/Result_0919/Choi/'
filename = paste(c(option,list_to_arg_vector(ChoiFitted_arg),list_to_arg_vector(coord_arg),format(Sys.time(), "%b%d_%X")),collapse = '_')
result_path = paste0(result_dir,filename)
dir.create(result_path)


######################################################################################
#####################    Coord Generation    #################################################

Data.Preparation.result = Data_Preparation_Fixed_range(option,par.cov,d=coord_arg$d,N,
                                                       r=coord_arg$r,n=coord_arg$n,
                                                       l_max=coord_arg$l_max,seed = 0)

save(Data.Preparation.result, file=paste0(result_path,"/Data.Preparation.result.RData"))

start_time = Sys.time()


all.result <- foreach(k = seq(from=1,to=N,by=1),.packages=dependency_packages,.combine='rbind') %dopar%
  #for(k in seq(from=1,to=N,by=1))
  {
    
    tryCatch({
      print(k)
      
      ######################################################################################
      #####################    Choi's Method    #################################################
      
      coord_matrix = Data.Preparation.result$Output.list$coord_matrix
      Y = Data.Preparation.result$Output.list$Y_matrix[k,,]
      
      ## empirical stat
      Choi_result_orig = Choi_Empirical_Cov(Y,coord_matrix)
      Choi_result_orig$empirical_covhat = Choi_result_orig$sum_covhat/Choi_result_orig$n_h
      Choi_result_orig$weight = Choi_result_orig$n_h/(1-Choi_result_orig$empirical_covhat)^2
      
      ## basis value
      Bspline_basis <- matrix(0,nrow=length(Choi_result_orig$dist),ncol=ChoiFitted_arg$m+ChoiFitted_arg$p)
      for(j in 1:(ChoiFitted_arg$m+ChoiFitted_arg$p)){
        Bspline_basis[,j] = sapply(Choi_result_orig$dist^2,Bspline_f,j,
                                   l=ChoiFitted_arg$p-1,p=ChoiFitted_arg$p,m=ChoiFitted_arg$m)
      }
      
      
      optout_BsplineChoi <- try(optim(par=rep(1/(ncol(Bspline_basis)),ncol(Bspline_basis)),Weighted_MSE,method='L-BFGS-B',
                                      y = Choi_result_orig$empirical_covhat, X=Bspline_basis,weight = Choi_result_orig$weight,
                                      lower=rep(0,ncol(Bspline_basis))),silent = FALSE)
      
      if('try-error' %in% class(optout_BsplineChoi)){
        optout_BsplineChoi = list(par=NULL)
        #success_indicator[1] = FALSE
      }else{
        
        ## to calculate metric
        f_cov_true <- function(x){
          Cov_function(x,option = option,par.cov=par.cov)
        }
        f_cor_true <- function(x){
          Cov_function(x,option = option,par.cov=c(1,par.cov[-1]))
        }
        f_variog_true <- function(x){
          par.cov[1] - Cov_function(x,option = option,par.cov=par.cov)
        }
        
        
        f_cov_est <- function(rho){
          Cov_BsplineChoi(rho,beta=optout_BsplineChoi$par,p=ChoiFitted_arg$p,m=ChoiFitted_arg$m)
        }
        f_cor_est <- function(rho){
          Cov_BsplineChoi(rho,beta=optout_BsplineChoi$par,p=ChoiFitted_arg$p,m=ChoiFitted_arg$m)/Cov_BsplineChoi(0,beta=optout_BsplineChoi$par,p=ChoiFitted_arg$p,m=ChoiFitted_arg$m)
        }
        
        f_variog_est <- function(rho){
          Cov_BsplineChoi(0,beta=optout_BsplineChoi$par,p=ChoiFitted_arg$p,m=ChoiFitted_arg$m) - 
            Cov_BsplineChoi(rho,beta=optout_BsplineChoi$par,p=ChoiFitted_arg$p,m=ChoiFitted_arg$m)
        }
        
        optout_BsplineChoi$C_zero = f_cov_est(0)
        
        #cov metric
        optout_BsplineChoi$norm2_cov_integration =  try(Norm2_between_f(f_cov_true,f_cov_est,0,metric_upper_bound)/sqrt(metric_upper_bound),silent = FALSE)
        if('try-error' %in% class(optout_BsplineChoi$norm2_cov_integration)){
          optout_BsplineChoi$norm2_cov_integration = NA
        }
        optout_BsplineChoi$norm2_cov_grid = sqrt(sum((f_cov_true(metric_grid) - f_cov_est(metric_grid))^2 * metric_grid_space)/metric_upper_bound) 
        optout_BsplineChoi$supnorm_cov_grid = max(abs(f_cov_true(metric_grid) - f_cov_est(metric_grid)))
        
        #cor metric
        optout_BsplineChoi$norm2_cor_integration =  try(Norm2_between_f(f_cor_true,f_cor_est,0,metric_upper_bound)/sqrt(metric_upper_bound),silent = FALSE)
        if('try-error' %in% class(optout_BsplineChoi$norm2_cor_integration)){
          optout_BsplineChoi$norm2_cor_integration = NA
        }
        optout_BsplineChoi$norm2_cor_grid = sqrt(sum((f_cor_true(metric_grid) - f_cor_est(metric_grid))^2 * metric_grid_space)/metric_upper_bound) 
        optout_BsplineChoi$supnorm_cor_grid = max(abs(f_cor_true(metric_grid) - f_cor_est(metric_grid)))
        
        #variog metric
        optout_BsplineChoi$norm2_variog_integration =  try(Norm2_between_f(f_variog_true,f_variog_est,0,metric_upper_bound)/sqrt(metric_upper_bound),
                                                           silent = FALSE)
        if('try-error' %in% class(optout_BsplineChoi$norm2_variog_integration)){
          optout_BsplineChoi$norm2_variog_integration = NA
        }
        optout_BsplineChoi$norm2_variog_grid = sqrt(sum((f_variog_true(metric_grid) - f_variog_est(metric_grid))^2 * metric_grid_space)/metric_upper_bound) 
        optout_BsplineChoi$supnorm_variog_grid = max(abs(f_variog_true(metric_grid) - f_variog_est(metric_grid))) 
        
      }
      
      
      ######################################################################################
      #####################    Parametric Method    #################################################
      
      Matern_est = Parametric_est_Choi_model_metric(true_option=option,true_par.cov = par.cov,
                                                    est_option = 'Matern', 
                                                    h = Choi_result_orig$dist,Cov_hat = Choi_result_orig$empirical_covhat,weight = Choi_result_orig$weight,
                                                    metric_upper_bound = metric_upper_bound,metric_grid_space = metric_grid_space, metric_grid = metric_grid)
      
      Cauchy_est = Parametric_est_Choi_model_metric(true_option=option,true_par.cov = par.cov,
                                                    est_option = 'Cauchy', 
                                                    h = Choi_result_orig$dist,Cov_hat = Choi_result_orig$empirical_covhat,weight = Choi_result_orig$weight,
                                                    metric_upper_bound = metric_upper_bound,metric_grid_space = metric_grid_space, metric_grid = metric_grid)
      
      Gaussian_est = Parametric_est_Choi_model_metric(true_option=option,true_par.cov = par.cov,
                                                      est_option = 'Gaussian', 
                                                      h = Choi_result_orig$dist,Cov_hat = Choi_result_orig$empirical_covhat,weight = Choi_result_orig$weight,
                                                      metric_upper_bound = metric_upper_bound,metric_grid_space = metric_grid_space, metric_grid = metric_grid)
      
      GenCauchy_est = Parametric_est_Choi_model_metric(true_option=option,true_par.cov = par.cov,
                                                       est_option = 'GenCauchy', 
                                                       h = Choi_result_orig$dist,Cov_hat = Choi_result_orig$empirical_covhat,weight = Choi_result_orig$weight,
                                                       metric_upper_bound = metric_upper_bound,metric_grid_space = metric_grid_space, metric_grid = metric_grid)
      
      Parametric_result_list <- list(Matern_est = Matern_est,
                                     Cauchy_est = Cauchy_est,
                                     Gaussian_est = Gaussian_est,
                                     GenCauchy_est = GenCauchy_est)
      ######################################################################################
      #####################    Record Results    #################################################
      
      if(do_plot & !is.null(optout_BsplineChoi$par)){
        pdf(file=paste0(result_path,'/',k,'_',format(Sys.time(), "%b%d_%X"),'_plot.pdf'))
        ## plot fitted result
        d_seq = seq(from=0,to=max(Choi_result_orig$dist)*1.2,length.out = 100)
        ## empirical cov_hat
        plot(Choi_result_orig$dist,Choi_result_orig$empirical_covhat,type='p',pch=1,
             xlab='dist',ylab='cov value',
             main=paste('ChoiNP method of ',option),
             ylim=c(min(Choi_result_orig$empirical_covhat),max(1,Choi_result_orig$empirical_covhat)))
        ## true
        lines(d_seq,sapply(d_seq, Cov_function,option = option, par.cov = par.cov),col='black',lty='solid',lwd=5)
        ## Nonparametric
        lines(d_seq,sapply(d_seq, f_cov_est),col='red',lty='dashed',lwd=5)
        ## Parametric
        all_color = c('red','blue','green','purple','orange','pink','yellow')
        for(i in 1:length(Parametric_result_list)){
          Current_list = Parametric_result_list[[i]]
          if(!is.null(Current_list$par)){
            lines(d_seq, sapply(d_seq,Cov_function,option = Current_list$option,par.cov=Current_list$par),col=all_color[i+1],lty='dashed',lwd=2)
          }
        }
        
        legend('topright',legend=c('empirical cov','true cov','Bspline Choi','Matern','Cauchy','Gaussian','GenCauchy'),
               col = c('black','black',all_color),
               pch = c(20,NA,rep(NA,length(all_color))),
               lty = c(NA,'solid',rep('dashed',length(all_color))),bty='n')
        
        dev.off()
      }
      
      #####################################################################################
      #####################    Save Results    ############################################
      
      
      BsplineChoi_result_list <- list(optout_BsplineChoi = optout_BsplineChoi,
                                      Choi_result_orig = Choi_result_orig,
                                      args = list(ChoiFitted_arg = ChoiFitted_arg,
                                                  coord_arg = coord_arg))
      
      
      save(BsplineChoi_result_list,Parametric_result_list,
           file=paste0(result_path,'/',k,'_',format(Sys.time(), "%b%d_%X"),'.RData'))
      
      one.sample.list = list()
      one.sample.list[[1]] = list(BsplineChoi_result_list = BsplineChoi_result_list,
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


## save coef of nonparametric estimation
Nonparametric_coef_matrix = t(sapply(all.result, function(x){
  if(is.null(x$BsplineChoi_result_list$optout_BsplineChoi$par)){
    return(rep(NA,ChoiFitted_arg$m+ChoiFitted_arg$p))
  }else{
    return(x$BsplineChoi_result_list$optout_BsplineChoi$par)
  }
}))

coef.file = paste0(result_path,'/Bspline_coef_est_all','.txt')
write.table(Nonparametric_coef_matrix, file = coef.file,
            row.names = paste0('result.',1:N))


## save metric of nonparametric estimation
metric_colnames = c('C_zero',
                    'norm2_cov_grid','norm2_cov_integration','supnorm_cov_grid',
                    'norm2_cor_grid','norm2_cor_integration','supnorm_cor_grid',
                    'norm2_variog_grid','norm2_variog_integration','supnorm_variog_grid')

Nonparametric_metric_matrix = t(sapply(all.result, function(x){
  optout_BsplineChoi = x$BsplineChoi_result_list$optout_BsplineChoi
  if(is.null(optout_BsplineChoi$par)){
    return(rep(NA,10))
  }else{
    return(c(optout_BsplineChoi$C_zero,
             optout_BsplineChoi$norm2_cov_grid,optout_BsplineChoi$norm2_cov_integration,optout_BsplineChoi$supnorm_cov_grid,
             optout_BsplineChoi$norm2_cor_grid,optout_BsplineChoi$norm2_cor_integration,optout_BsplineChoi$supnorm_cor_grid,
             optout_BsplineChoi$norm2_variog_grid,optout_BsplineChoi$norm2_variog_integration,optout_BsplineChoi$supnorm_variog_grid))
  }
}))
colnames(Nonparametric_metric_matrix) = metric_colnames
nonparametric_metric.file = paste0(result_path,'/Bspline_metric_all','.txt')
write.table(Nonparametric_metric_matrix, file = nonparametric_metric.file,
            row.names = paste0('result.',1:length(all.result)),col.names = TRUE)


## save Parametric estimate value
parametric_order = c('Matern_est','Cauchy_est','Gaussian_est','GenCauchy_est')
param_metric_matrix = t(sapply(all.result, function(x){
  if(!is.null(x$BsplineChoi_result_list$optout_BsplineChoi$par)){ #sucess
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
      sapply(x$Parametric_result_list, function(opt_est){opt_est$supnorm_variog_grid}))
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
                                  paste0(parametric_order,'.supnorm_variog_grid'))

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


#}

stopCluster(cl)


