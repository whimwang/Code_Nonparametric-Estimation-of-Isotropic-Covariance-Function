
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
source(paste0(code.path,'Code_other_method/','Huang_code/Parametric_Huang_est_with_nugget.R')) #param with nugget


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


## specify fitting args
HuangFitted_arg = list(v=10,L=200,fit_ratio = 1,dist.range_parameter = 3*(10^5),cov.distort = 100)
do_plot = TRUE
#Since irregular, r=200 replicates should ensure a relative stable estimation of variogram.


result_dir = '/mnt/home/ywang225/Cov_research/Result_1010/Huang_real_data/'
filename = paste(c(list_to_arg_vector(ChoiFitted_arg),format(Sys.time(), "%b%d_%X")),collapse = '_')
result_path = paste0(result_dir,filename)
dir.create(result_path)


input_data.path = "/Users/wangyiming/Documents/2019 Research Daily/research_code/Code_2020/Updated0620/Code0620/input_real_data/Year117_Station241_49_27_60_100/"
#Server
#input_data.path = '/mnt/home/ywang225/Cov_research/Updated1010/Code0620/input_real_data/Year117_Station241_49_27_60_100/'
load(paste0(input_data.path,"Year117_Station241_49_27_60_100_input_data.RData"))

lon_lat_data = input_data$lon_lat_data
x_y_scale = input_data$x_y_scale
observation_residual = input_data$observation_residual
dist_select = input_data$dist_select

n = ncol(observation_residual)
r = nrow(observation_residual)
idx = 1:20

#cov_with_nugget = mean(apply(observation_residual[,idx],2,var))
cov_with_nugget = mean(apply(observation_residual[,idx],2,function(x){mean(x^2)}))
dist.range_parameter = HuangFitted_arg$dist.range_parameter
cov.distort = HuangFitted_arg$cov.distort

      
      ## empirical stat
      Huang_result_orig_converted = Huang_Empirical_Variogram(Y=observation_residual[,idx],coord_matrix=lon_lat_data[idx,],dist_type='distCosine')
      Huang_result_orig_converted$dist = Huang_result_orig_converted$dist / dist.range_parameter
      Huang_result_orig_converted$sum_varioghat = Huang_result_orig_converted$sum_varioghat/cov.distort
      Huang_result_orig_converted$empirical_varioghat = Huang_result_orig_converted$sum_varioghat/Huang_result_orig_converted$n_h
      Huang_result_orig_converted$weight = Huang_result_orig_converted$n_h
      
      dist_non_zero_index = which(Huang_result_orig_converted$dist!=0)
      
      Huang_result_orig_converted = list(dist = Huang_result_orig_converted$dist[dist_non_zero_index],
                               sum_varioghat = Huang_result_orig_converted$sum_varioghat[dist_non_zero_index],
                               n_h = Huang_result_orig_converted$n_h[dist_non_zero_index],
                               empirical_varioghat = Huang_result_orig_converted$empirical_varioghat[dist_non_zero_index],
                               weight = Huang_result_orig_converted$weight[dist_non_zero_index])
      
    
      ## specify fitting args
      v_range <- HuangFitted_arg$v
      L_grid <- HuangFitted_arg$L
    
      
      
      n_dist = length(Huang_result_orig_converted$dist)
      n_fit = as.integer(n_dist* HuangFitted_arg$fit_ratio)
      
      fitting_dist = Huang_result_orig_converted$dist[1:n_fit]
      fitting_variog = Huang_result_orig_converted$empirical_varioghat[1:n_fit]
      
      X <- Xnew <- Get_Bmatrix(dist = fitting_dist,v=v_range,L=L_grid,nugget = TRUE) #B in paper
      Dc <- Dnew <- Get_PenalizeKmatrix(v=v_range,L=L_grid,nugget = TRUE)
      weight = Huang_result_orig_converted$n_h[1:n_fit]
      WW = diag(weight)
      XnewWWXnew <- t(Xnew)%*%WW%*%Xnew
      dvec = t(Xnew) %*% WW %*% fitting_variog
      Amat = diag(L_grid+1)
      bvec = rep(0,L_grid+1)
      
      ## choose lamdba
      ngcv = 20
      f_semivariog_with_nugget_est_zero <- nugget_est <- lossf <- lossgamma <- gcvnew <- tr1 <- tr2 <- numeric(ngcv)
      lambda_seq = sapply((1:ngcv)-1, function(x){10^(6*(x-1)/(ngcv-1))})
      
      
      pdf(file=paste0(result_path,'/Huang_tuning_parameters_',format(Sys.time(), "%b%d_%X"),'_plot.pdf'))
      
      for(igcv in 1:ngcv){
        
        print(paste0('cross-validation: ',igcv))
        print(paste0('lambda value: ',lambda_seq[igcv]))
        lambda = lambda_seq[igcv]
        
        ### solve QP without nugget
        Dmat <- XnewWWXnew + lambda*Dnew
        QPout <- solve.QP(Dmat,dvec,Amat,bvec=bvec)
        nugget_est[igcv] = QPout$solution[L_grid+1]/2 #2nugget
        fest <- QPout$solution[1:L_grid]
        
        variog_hat <- Xnew%*%QPout$solution #in paper: Bv or say l_ij g_w
        
        
        if(do_plot){
          if(igcv == 1){
            plot(Huang_result_orig_converted$dist*dist.range_parameter/10^4, 
                 Huang_result_orig_converted$empirical_varioghat/2 * cov.distort,
                 xlab='dist/(10^4)',ylab='semi-variog value',ylim=c(0,max(Huang_result_orig_converted$empirical_varioghat * cov.distort/2)*1.2),
                 main=paste0('Semivariog'))
            legend('topright',legend = c('empirical variog','true variog','est variog','min.observed.dist'),
                   col = c('black','black','red','blue'),
                   pch = c(1,NA,NA,NA),
                   lty=c(NA,1,1,1),bty='n')
            abline(v=min(Huang_result_orig_converted$dist)*dist.range_parameter/10^4,col='blue',lwd=3)
            #true
            #points(Huang_result_orig_converted$dist, sapply(Huang_result_orig_converted$dist,GenCov_list$variog_function),col='red',pch=20)
          }
          
          d_seq = seq(from=0,to=max(Huang_result_orig_converted$dist)*dist.range_parameter,length.out = 100)
          f_semivariog_est <- function(dist){ #semivariog including nugget value 
            Empirical_Variog_est(dist/dist.range_parameter,f_est=QPout$solution,v=HuangFitted_arg$v,L=HuangFitted_arg$L,nugget = TRUE) * cov.distort/2 
          }
          lines(d_seq/10^4,sapply(d_seq, f_semivariog_est),col='red',lwd=2)
          
          #lines(fitting_dist,variog_hat,col='red',lwd=5)
        }
        
        f_semivariog_with_nugget_est_zero[igcv] = f_semivariog_est(0)
        
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
        
      }
      
      dev.off()
      
      ### QP paramters
      
      ## save lambda
      lambda_list = list(lambda = lambda_seq)
      
      ## save QP parameters
      QP_list = list(X=X,Xnew=Xnew,Dc=Dc,Dnew=Dnew,weight=weight,WW=WW,dvec=dvec,Amat=Amat,bvec=bvec)
      
      nugget_est = nugget_est * cov.distort
      best_igcv = which(gcvnew == min(gcvnew))
      best_lambda = lambda_seq[best_igcv]
      best_nugget_est = nugget_est[best_igcv]
      best_semivariog_est_zero = f_semivariog_with_nugget_est_zero[best_igcv] #with nugget
      
      
      ########################################
      ####### Metric of Best solution ##########
      
      ### solve QP first time
      Dmat <- t(Xnew)%*%WW%*%Xnew + best_lambda*Dnew
      QPout <- solve.QP(Dmat,dvec,Amat,bvec=bvec)
      Huang_fest <- QPout$solution[1:L_grid]
      gammahat <- Xnew%*%QPout$solution #in paper: Bv or say l_ij g_w
      
      ################################ 
      ######## Define est function  ##### 
      f_semivariog_est_Huang <- function(dist){ #semivariog with nugget value 
        Empirical_Variog_est(dist/dist.range_parameter,f_est=QPout$solution,v=HuangFitted_arg$v,L=HuangFitted_arg$L,nugget = TRUE) * cov.distort/2 
      }
      
      #Version 1: empirical at max obs dist
      best_C_zero_with_nugget_version1 = f_semivariog_est_Huang(max(Huang_result_orig_converted$dist * dist.range_parameter))
      best_C_zero_version1 = best_C_zero_with_nugget_version1 - best_nugget_est
      f_cov_est_Huang_version1 <- function(dist){
        f_semivariog_est_Huang(0) - f_semivariog_est_Huang(dist) + best_C_zero_version1
      }
      #Version 2: cov using 
      best_C_zero_with_nugget_version2 = cov_with_nugget
      best_C_zero_version2 = cov_with_nugget - best_nugget_est
      f_cov_est_Huang_version2 <- function(dist){
        f_semivariog_est_Huang(0) - f_semivariog_est_Huang(dist) + best_C_zero_version2
      }
      
      
      
      
   
     
      
      ######################################################################################
      #####################    Parametric Method    #################################################
      if(TRUE){
      
      optout_Matern_Huang = Matern_Parametric_est_with_nugget_Huang_model(h = Huang_result_orig_converted$dist,
                                                                 Variog_hat = Huang_result_orig_converted$empirical_varioghat,
                                                                 weight = Huang_result_orig_converted$weight)
      
      optout_Gaussian_Huang = Gaussian_Parametric_est_with_nugget_Huang_model(h = Huang_result_orig_converted$dist,
                                                                 Variog_hat = Huang_result_orig_converted$empirical_varioghat,
                                                                 weight = Huang_result_orig_converted$weight)
      
      
      optout_Cauchy_Huang = Cauchy_Parametric_est_with_nugget_Huang_model(h = Huang_result_orig_converted$dist,
                                                                 Variog_hat = Huang_result_orig_converted$empirical_varioghat,
                                                                 weight = Huang_result_orig_converted$weight)
      
      
      optout_GenCauchy_Huang = GenCauchy_Parametric_est_with_nugget_Huang_model(h = Huang_result_orig_converted$dist,
                                                                 Variog_hat = Huang_result_orig_converted$empirical_varioghat,
                                                                 weight = Huang_result_orig_converted$weight)
      
      
      param_est_option = c('Matern','Cauchy','Gaussian','GenCauchy')
      param_est_nugget_Huang = rep(NA,4)
      param_est_vector_Huang = rep(NA,17)
      param_est_C_zero_with_nugget_Huang = rep(NA,4)
  
      
      param_est_vector_Huang[c(1,2)] = c(dist.range_parameter,cov.distort)
      if(!is.na(optout_Matern_Huang$par)){
        param_est_vector_Huang[3:6] = optout_Matern_Huang$par
        param_est_nugget_Huang[1] = optout_Matern_Huang$par[1]
        param_est_C_zero_with_nugget_Huang[1] = sum(optout_Matern_Huang$par[1:2])
      }
      if(!is.na(optout_Cauchy_Huang$par)){
        param_est_vector_Huang[7:9] = optout_Cauchy_Huang$par
        param_est_nugget_Huang[2] = optout_Cauchy_Huang$par[1]
        param_est_C_zero_with_nugget_Huang[2] = sum(optout_Cauchy_Huang$par[1:2])
      }
      if(!is.na(optout_Gaussian_Huang$par)){
        param_est_vector_Huang[10:12] = optout_Gaussian_Huang$par
        param_est_nugget_Huang[3] = optout_Gaussian_Huang$par[1]
        param_est_C_zero_with_nugget_Huang[3] = sum(optout_Gaussian_Huang$par[1:2])
      }
      if(!is.na(optout_GenCauchy_Huang$par)){
        param_est_vector_Huang[13:17] = optout_GenCauchy_Huang$par
        param_est_nugget_Huang[4] = optout_GenCauchy_Huang$par[1]
        param_est_C_zero_with_nugget_Huang[4] = sum(optout_GenCauchy_Huang$par[1:2])
      }
      
      names(param_est_vector_Huang) = c('dist.range_parameter','cov.distort',
                                                paste0('Matern.',1:4),paste0('Cauchy.',1:3),paste0('Gaussian.',1:3),paste0('GenCauchy.',1:5))
      names(param_est_nugget_Huang) = param_est_option
      names(param_est_C_zero_with_nugget_Huang) = param_est_option
      param_est_C_zero_Huang = param_est_C_zero_with_nugget_Huang - param_est_nugget_Huang
      }
      
      
      ######## Plot Huang vs Parametric Variogram #########
      if(TRUE){
        
        pdf(file=paste0(result_path,'/Huang_vs_Parametric_',format(Sys.time(), "%b%d_%X"),'_plot.pdf'))
        
        d_seq = seq(from=0,to=5*(10^6),length.out = 100)
      
        #empirical vs HuangFit:
        plot(Huang_result_orig_converted$dist * dist.range_parameter/10^4,
             Huang_result_orig_converted$empirical_varioghat * cov.distort,
             ylim=c(0,max(Huang_result_orig_converted$empirical_varioghat) * cov.distort*1.2),
             xlim=c(0,max(d_seq)/10^4),
             pch=20,xlab='dist/10^4',ylab='Variog',main='Variog with nugget:Huang')
        
        all_color = 2:6
        
        lines(d_seq/10^4,
             (sapply(d_seq/dist.range_parameter, Variog_function, option = 'Matern',par.cov = optout_Matern_Huang$par[-1]) + 2*optout_Matern_Huang$par[1]) * cov.distort,
             lty='dashed',col=all_color[1],lwd=2)

        lines(d_seq/(10^4),(sapply(d_seq/dist.range_parameter, Variog_function, option = 'Cauchy', par.cov = optout_Cauchy_Huang$par[-1]) + 2*optout_Cauchy_Huang$par[1]) * cov.distort,
              lwd=2,lty='dashed',col=all_color[2])
        
        lines(d_seq/(10^4),(sapply(d_seq/dist.range_parameter, Variog_function, option = 'Gaussian', par.cov = optout_Gaussian_Huang$par[-1])  + 2*optout_Gaussian_Huang$par[1]) * cov.distort,
              lwd=2,lty='dashed',col=all_color[3])
        
        lines(d_seq/(10^4),(sapply(d_seq/dist.range_parameter, Variog_function, option = 'GenCauchy', par.cov = optout_GenCauchy_Huang$par[-1]) + 2*optout_GenCauchy_Huang$par[1]) * cov.distort,
              lwd=2,lty='dashed',col=all_color[4])
        
        lines(d_seq/(10^4),sapply(d_seq, f_semivariog_est_Huang) * 2,lwd=3,col=all_color[5])
        
        legend('topleft',
               legend = c('empirical variog',paste0(c(param_est_option,'Huang'),' nugget:',c(round(c(param_est_nugget_Huang * cov.distort,best_nugget_est),digits = 3)))),
               pch = c(20,rep(NA,5)),
               lty = c(NA,rep('dashed',4),'solid'),
               col = c('black',all_color),bty='n')
        
        
        #Cov without nugget
        all_color = 2:7
        plot(d_seq/10^4,
             sapply(d_seq/dist.range_parameter, Cov_function, option = 'Matern',par.cov = optout_Matern_Huang$par[-1]) * cov.distort,
             type='l',lty='dashed',lwd=2,col=all_color[1],
             xlab='dist/10^4',ylab='Cov',main="Cov without nugget",
             ylim=c(0,max(param_est_C_zero_Huang)*cov.distort * 1.2))
        
        lines(d_seq/(10^4),sapply(d_seq/dist.range_parameter, Cov_function, option = 'Cauchy', par.cov = optout_Cauchy_Huang$par[-1]) * cov.distort,
              lwd=2,lty='dashed',col=all_color[2])
        
        lines(d_seq/(10^4),sapply(d_seq/dist.range_parameter, Cov_function, option = 'Gaussian', par.cov = optout_Gaussian_Huang$par[-1]) * cov.distort,
              lwd=2,lty='dashed',col=all_color[3])
        
        lines(d_seq/(10^4),sapply(d_seq/dist.range_parameter, Cov_function, option = 'GenCauchy', par.cov = optout_GenCauchy_Huang$par[-1]) * cov.distort,
              lwd=2,lty='dashed',col=all_color[4])
        
        #nonparametric
        lines(d_seq/(10^4),sapply(d_seq,f_cov_est_Huang_version1), 
              lwd=2,lty='dashed',col=all_color[5])
        
        lines(d_seq/(10^4),sapply(d_seq,f_cov_est_Huang_version2), 
              lwd=2,lty='dashed',col=all_color[6])
        
        legend('topleft',
               legend = c(paste0(c(param_est_option,'Huang max dist','Huang cov'),' nugget:',c(round(c(param_est_nugget_Huang * cov.distort,best_nugget_est,best_nugget_est),digits = 3)))),
               lty = c(rep('dashed',4),'solid','solid'),
               col = c('black',all_color),bty='n')
        
        dev.off()
        
      }
      
      
      
      ###### save results #######
      
      Huang_result = list(dist.range_parameter = dist.range_parameter,cov.distort = cov.distort,
                                          par = QPout$solution,cov_with_nugget = cov_with_nugget,
                                          HuangFitted_arg = HuangFitted_arg)
      
      Huang_tuning_result = list(lambda_seq = lambda_seq, nugget_est = nugget_est)
      
      
      write_Huang_result_func <- function(Huang_result,filename){
        cat(paste("dist.range_paramter",Huang_result$dist.range_parameter),file=filename, sep="\n")
        cat(paste("cov.distort",Huang_result$cov.distort),file=filename, sep="\n",append = TRUE)
        cat(paste("cov_with_nugget",Huang_result$cov_with_nugget),file=filename, sep="\n",append = TRUE)
        cat(paste("v",Huang_result$HuangFitted_arg$v,"L",Huang_result$HuangFitted_arg$L,'fit_ratio',Huang_result$HuangFitted_arg$fit_ratio),file=filename, sep="\n", append = TRUE)
      }
      
      write_Huang_result_func(Huang_result,filename=paste0(result_path,'/','Huang_result.txt'))
      save(Huang_result,file=paste0(result_path,'/','Huang_result.RData'))
      save(Huang_result_orig_converted,file=paste0(result_path,'/','Huang_result_orig_converted.RData'))
      save(Huang_tuning_result,file=paste0(result_path,'/','Huang_tuning_result.RData'))
     
      write.table(param_est_C_zero_with_nugget_Huang,file=paste0(result_path,'/','Huang_Parametric_C_zero_with_nugget_result.txt'),col.names = FALSE,row.names = TRUE)
      write.table(param_est_nugget_Huang,file=paste0(result_path,'/','Huang_Parametric_nugget_result.txt'),col.names = FALSE,row.names = TRUE)
      write.table(param_est_vector_Huang,file=paste0(result_path,'/','Huang_Parametric_vector_result.txt'),col.names = FALSE,row.names = TRUE)
      
      
      