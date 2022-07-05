library('Rcpp')
library('inline')
library('gstat')
library('geoR')
library('sp')
#
code.path = '/Users/wangyiming/Documents/2019 Research Daily/research_code/Code_2020/Updated0620/'
#Server
#code.path = "/mnt/home/ywang225/Cov_research/Updated1011/"
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



code.path = '/Users/wangyiming/Documents/2019 Research Daily/research_code/Code_2020/Updated0620/Code0620/'
#server
#code.path = "/mnt/home/ywang225/Cov_research/Updated1022/Code0620/"
source(paste0(code.path,"Code0912-0919/Part-III.Objective.functions.R"))


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
HuangFitted_arg = list(v=10,L=200,fit_ratio = 1,cov.distort = 4^2,dist.range_parameter = 1*(10^5))
do_plot = TRUE
approx.option = 2
approx.par = 10^(-9)
#Since irregular, r=200 replicates should ensure a relative stable estimation of variogram.


#result_dir = '/mnt/home/ywang225/Cov_research/Result_1022/Huang_real_data/'
result_dir = '/Users/wangyiming/Documents/2019 Research Daily/research_code/Code_2020/Result_1022/Huang_real_data/'
filename = paste(c(list_to_arg_vector(HuangFitted_arg),format(Sys.time(), "%b%d_%X")),collapse = '_')
result_path = paste0(result_dir,filename)
dir.create(result_path)

print(result_path)


input_data.path = "/Users/wangyiming/Documents/2019 Research Daily/research_code/Code_2020/Updated0620/Code0620/input_real_data/Year40_Station189_46_30_80_100/"
#Server
#input_data.path = '/mnt/home/ywang225/Cov_research/Updated1016/Code0620/input_real_data/Year30_Station276_46_30_80_100/'
load(paste0(input_data.path,"Year40_Station189_46_30_80_100.RData"))

lon_lat_data = input_data$lon_lat_data
x_y_scale = input_data$x_y_scale
observation_residual = input_data$observation_residual
dist_select = input_data$dist_select

n = ncol(observation_residual)
r = nrow(observation_residual)
print(paste0('n:  ',n,'  r:',r))
idx = 1:n




start.time = Sys.time()

cov_with_nugget = mean(apply(observation_residual[,idx],2,var))
#cov_with_nugget = mean(apply(observation_residual[,idx],2,function(x){mean(x^2)}))
dist.range_parameter = HuangFitted_arg$dist.range_parameter
cov.distort = HuangFitted_arg$cov.distort

pdf(file=paste0(result_path,'/semivariogram_',format(Sys.time(), "%b%d_%X"),'_plot.pdf'))

  
cutoff = 210*(10^4) ; width = cutoff/100
i=1
one_replication = data.frame(obs = observation_residual[i,idx],x=x_y_scale$LONGITUDE.1[idx],y=x_y_scale$LATITUDE.1[idx])
coordinates(one_replication) = ~x+y
one_replication_variogram = variogram(obs~1,one_replication,cutoff = cutoff, width = width)
#plot(one_replication_variogram$dist/10^4,one_replication_variogram$gamma)

all_variogram_value = matrix(0,ncol=r+1,nrow=length(one_replication_variogram$dist))
colnames(all_variogram_value) = c("dist",rownames(observation_residual))
all_variogram_value[,1] = one_replication_variogram[,2]
all_variogram_value[,2] = one_replication_variogram[,3]

#250
plot(one_replication_variogram$dist/10^(4),one_replication_variogram$gamma,type='p',pch=20,ylim=c(0,400),xlab="dist (10^4)",ylab="semivariogram (100)")
lines(one_replication_variogram$dist/10^(4),one_replication_variogram$gamma,lwd=3,lty=)

for(i in 2:nrow(observation_residual)){
  one_replication = data.frame(obs = observation_residual[i,idx],x=x_y_scale$LONGITUDE.1[idx],y=x_y_scale$LATITUDE.1[idx])
  coordinates(one_replication) = ~x+y
  one_replication_variogram = variogram(obs~1,one_replication,cutoff = cutoff, width = width)
  if(all.equal(one_replication_variogram$dist,all_variogram_value[,1])){
    all_variogram_value[,i+1] = one_replication_variogram$gamma
  }else{
    print(i)
    stop('dist is different')
  }
  points(one_replication_variogram$dist/10^(4),one_replication_variogram$gamma,pch=20)
  lines(one_replication_variogram$dist/10^(4),one_replication_variogram$gamma,col=i,lwd=3)
}
lines(all_variogram_value[,1]/10^4,apply(all_variogram_value[,2:ncol(all_variogram_value)],1,mean),col='black',lwd=10,lty='dashed')
points(all_variogram_value[,1]/10^4,apply(all_variogram_value[,2:ncol(all_variogram_value)],1,mean),col='red',pch=20)

dev.off()
      
      Huang_result_orig_converted = list(dist=all_variogram_value[,1],
                    empirical_varioghat = apply(all_variogram_value[,2:ncol(all_variogram_value)],1,mean))
      Huang_result_orig_converted$weight = rep(1,length(Huang_result_orig_converted$dist))
      Huang_result_orig_converted$n_h = rep(1,length(Huang_result_orig_converted$dist))
      Huang_result_orig_converted$dist = Huang_result_orig_converted$dist/dist.range_parameter
      Huang_result_orig_converted$empirical_varioghat = Huang_result_orig_converted$empirical_varioghat/cov.distort
        
      ## empirical stat
      #Huang_result_orig_converted = Huang_Empirical_Variogram(Y=observation_residual[,idx],coord_matrix=lon_lat_data[idx,],dist_type='distCosine')
      #Huang_result_orig_converted$dist = Huang_result_orig_converted$dist / dist.range_parameter
      #Huang_result_orig_converted$sum_varioghat = Huang_result_orig_converted$sum_varioghat/cov.distort
      #Huang_result_orig_converted$empirical_varioghat = Huang_result_orig_converted$sum_varioghat/Huang_result_orig_converted$n_h
      #Huang_result_orig_converted$weight = Huang_result_orig_converted$n_h
      
      dist_non_zero_index = which(Huang_result_orig_converted$dist!=0)
      
      Huang_result_orig_converted = list(dist = Huang_result_orig_converted$dist[dist_non_zero_index],
                               #sum_varioghat = Huang_result_orig_converted$sum_varioghat[dist_non_zero_index],
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
        nugget_est[igcv] = QPout$solution[L_grid+1] #nugget
        fest <- QPout$solution[1:L_grid]
        
        variog_hat <- Xnew%*%QPout$solution #in paper: Bv or say l_ij g_w
        
        
        if(do_plot){
          if(igcv == 1){
            plot(Huang_result_orig_converted$dist*dist.range_parameter/10^4, 
                 Huang_result_orig_converted$empirical_varioghat * cov.distort,
                 xlim=c(0,max(dist_select)/10^4),
                 xlab='dist/(10^4)',ylab='semi-variog value',ylim=c(0,max(Huang_result_orig_converted$empirical_varioghat * cov.distort)*1.2),
                 main=paste0('Semivariog'))
            legend('topright',legend = c('empirical variog','true variog','est variog','min.observed.dist'),
                   col = c('black','black','red','blue'),
                   pch = c(1,NA,NA,NA),
                   lty=c(NA,1,1,1),bty='n')
            abline(v=min(Huang_result_orig_converted$dist)*dist.range_parameter/10^4,col='blue',lwd=3)
            #true
            #points(Huang_result_orig_converted$dist, sapply(Huang_result_orig_converted$dist,GenCov_list$variog_function),col='red',pch=20)
          }
          
          d_seq = seq(from=0,to=max(dist_select),length.out = 100)
          f_semivariog_est <- function(dist){ #semivariog including nugget value 
            Empirical_Variog_est(dist/dist.range_parameter,f_est=QPout$solution,v=HuangFitted_arg$v,L=HuangFitted_arg$L,nugget = TRUE) * cov.distort
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
      
     
      
      print(paste0('best igcv: ', best_igcv, ' best lambda: ', best_lambda,' best nugget: ', best_nugget_est))
      
      
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
        Empirical_Variog_est(dist/dist.range_parameter,f_est=QPout$solution,v=HuangFitted_arg$v,L=HuangFitted_arg$L,nugget = TRUE) * cov.distort
      }
      
      Huang_optout = list(best_igcv = best_igcv,
                          best_lambda = best_lambda,
                          best_nugget_est = best_nugget_est,
                          f_est = QPout$solution)
      
      #Version 1: empirical at max obs dist
      best_C_zero_with_nugget_version1 = f_semivariog_est_Huang(max(dist_select))
      best_C_zero_version1 = best_C_zero_with_nugget_version1 - best_nugget_est
      f_cov_est_Huang_version1 <- function(dist){
        f_semivariog_est_Huang(0) - f_semivariog_est_Huang(dist) + best_C_zero_version1
      }
      
      Huang_optout$best_C_zero_version1 = best_C_zero_version1
      
      Sigma_cov_est_Huang_version1 <- matrix(sapply(dist_select, f_cov_est_Huang_version1),nrow=n,ncol=n)/cov.distort+ 
        diag(rep(best_nugget_est,n),nrow=n,ncol=n)/cov.distort
      Huang_optout$deviance1  = Two.Neg.Log.Likelihood.Use.Sigma(Y=observation_residual/sqrt(cov.distort),Sigma = Sigma_cov_est_Huang_version1,
                                          approx.option = approx.option,approx.par = approx.par)
      
      #Version 2: cov using 
      best_C_zero_with_nugget_version2 = cov_with_nugget
      best_C_zero_version2 = cov_with_nugget - best_nugget_est
      f_cov_est_Huang_version2 <- function(dist){
        f_semivariog_est_Huang(0) - f_semivariog_est_Huang(dist) + best_C_zero_version2
      }
      Huang_optout$best_C_zero_version2 = best_C_zero_version2
      Huang_optout$cov_with_nugget = cov_with_nugget
      
      Sigma_cov_est_Huang_version2 <- matrix(sapply(dist_select, f_cov_est_Huang_version2),nrow=n,ncol=n)/cov.distort+ 
        diag(rep(best_nugget_est,n),nrow=n,ncol=n)/cov.distort
      Huang_optout$deviance2 = Two.Neg.Log.Likelihood.Use.Sigma(Y=observation_residual/sqrt(cov.distort),Sigma = Sigma_cov_est_Huang_version2,
                                       approx.option = approx.option,approx.par = approx.par)
      
      end.time = Sys.time()
      print(end.time-start.time)
      Huang_optout$time = end.time - start.time
      
   
     
      
      ######################################################################################
      #####################    Parametric Method    #################################################
      if(TRUE){
      #Matern
      time1 = Sys.time()
      optout_Matern_Huang = Matern_Parametric_est_with_nugget_Huang_model(h = Huang_result_orig_converted$dist,
                                                                 Variog_hat = Huang_result_orig_converted$empirical_varioghat * 2,
                                                                 weight = Huang_result_orig_converted$weight)
      Sigma_Matern <- matrix(sapply(dist_select/dist.range_parameter, Cov_function,option='Matern',par=optout_Matern_Huang$par[-1]),nrow=n,ncol=n)+ 
        diag(rep(optout_Matern_Huang$par[1],n),nrow=n,ncol=n)
      optout_Matern_Huang$deviance = Two.Neg.Log.Likelihood.Use.Sigma(Y=observation_residual/sqrt(cov.distort),Sigma = Sigma_Matern,
                                       approx.option = approx.option,approx.par = approx.par)
      optout_Matern_Huang$time = Sys.time()-time1
      
      print(time1 - Sys.time())
      
      #Cauchy
      time1 = Sys.time()
      optout_Cauchy_Huang = Cauchy_Parametric_est_with_nugget_Huang_model(h = Huang_result_orig_converted$dist,
                                                                          Variog_hat = Huang_result_orig_converted$empirical_varioghat * 2,
                                                                          weight = Huang_result_orig_converted$weight)
      Sigma_Cauchy <- matrix(sapply(dist_select/dist.range_parameter, Cov_function,option='Cauchy',par=optout_Cauchy_Huang$par[-1]),nrow=n,ncol=n)+ 
        diag(rep(optout_Cauchy_Huang$par[1],n),nrow=n,ncol=n)
      optout_Cauchy_Huang$deviance = Two.Neg.Log.Likelihood.Use.Sigma(Y=observation_residual/sqrt(cov.distort),Sigma = Sigma_Cauchy,
                                                                      approx.option = approx.option,approx.par = approx.par)
      
      optout_Cauchy_Huang$time = Sys.time()-time1
      print(time1 - Sys.time())
      
      
      #Gaussian
      time1 = Sys.time()
      optout_Gaussian_Huang = Gaussian_Parametric_est_with_nugget_Huang_model(h = Huang_result_orig_converted$dist,
                                                                          Variog_hat = Huang_result_orig_converted$empirical_varioghat * 2,
                                                                          weight = Huang_result_orig_converted$weight)
      Sigma_Gaussian <- matrix(sapply(dist_select/dist.range_parameter, Cov_function,option='Gaussian',par=optout_Gaussian_Huang$par[-1]),nrow=n,ncol=n)+ 
        diag(rep(optout_Gaussian_Huang$par[1],n),nrow=n,ncol=n)
      optout_Gaussian_Huang$deviance = Two.Neg.Log.Likelihood.Use.Sigma(Y=observation_residual/sqrt(cov.distort),Sigma = Sigma_Gaussian,
                                                                      approx.option = approx.option,approx.par = approx.par)
      optout_Gaussian_Huang$time = Sys.time()-time1
      print(time1 - Sys.time())
      
      #GenCauchy
      time1 = Sys.time()
      optout_GenCauchy_Huang = GenCauchy_Parametric_est_with_nugget_Huang_model(h = Huang_result_orig_converted$dist,
                                                                          Variog_hat = Huang_result_orig_converted$empirical_varioghat * 2,
                                                                          weight = Huang_result_orig_converted$weight)
      Sigma_GenCauchy <- matrix(sapply(dist_select/dist.range_parameter, Cov_function,option='GenCauchy',par=optout_GenCauchy_Huang$par[-1]),nrow=n,ncol=n)+ 
        diag(rep(optout_GenCauchy_Huang$par[1],n),nrow=n,ncol=n)
      optout_GenCauchy_Huang$deviance = Two.Neg.Log.Likelihood.Use.Sigma(Y=observation_residual/sqrt(cov.distort),Sigma = Sigma_GenCauchy,
                                                                      approx.option = approx.option,approx.par = approx.par)
      optout_GenCauchy_Huang$time = Sys.time()-time1
      print(time1 - Sys.time())

      
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
        
        d_seq = seq(from=0,to=max(dist_select)*1.2,length.out = 100)
      
        #empirical vs HuangFit:
        plot(Huang_result_orig_converted$dist * dist.range_parameter/10^4,
             Huang_result_orig_converted$empirical_varioghat * cov.distort,
             ylim=c(0,max(Huang_result_orig_converted$empirical_varioghat) * cov.distort*1.2),
             xlim=c(0,max(d_seq)/10^4),
             pch=20,xlab='dist/10^4',ylab='SemiVariog',main='SemiVariog with nugget:Huang')
        
        all_color = 2:6
        
        lines(d_seq/10^4,
             (sapply(d_seq/dist.range_parameter, Variog_function, option = 'Matern',par.cov = optout_Matern_Huang$par[-1]) + 2*optout_Matern_Huang$par[1]) * cov.distort/2,
             lty='dashed',col=all_color[1],lwd=2)

        lines(d_seq/(10^4),(sapply(d_seq/dist.range_parameter, Variog_function, option = 'Cauchy', par.cov = optout_Cauchy_Huang$par[-1]) + 2*optout_Cauchy_Huang$par[1]) * cov.distort/2,
              lwd=2,lty='dashed',col=all_color[2])
        
        lines(d_seq/(10^4),(sapply(d_seq/dist.range_parameter, Variog_function, option = 'Gaussian', par.cov = optout_Gaussian_Huang$par[-1])  + 2*optout_Gaussian_Huang$par[1]) * cov.distort/2,
              lwd=2,lty='dashed',col=all_color[3])
        
        lines(d_seq/(10^4),(sapply(d_seq/dist.range_parameter, Variog_function, option = 'GenCauchy', par.cov = optout_GenCauchy_Huang$par[-1]) + 2*optout_GenCauchy_Huang$par[1]) * cov.distort/2,
              lwd=2,lty='dashed',col=all_color[4])
        
        lines(d_seq/(10^4),sapply(d_seq, f_semivariog_est_Huang),lwd=3,col=all_color[5])
        
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
               col = c(all_color),bty='n')
        
        dev.off()
        
      }
      
      
      
      ###### save results #######
      
      Huang_result = list(dist.range_parameter = dist.range_parameter,cov.distort = cov.distort,
                                          par = QPout$solution,cov_with_nugget = cov_with_nugget,
                                          HuangFitted_arg = HuangFitted_arg)
      
      Huang_tuning_result = list(lambda_seq = lambda_seq, nugget_est = nugget_est)
      
      Huang_parametric_result_list = list(optout_Matern_Huang = optout_Matern_Huang,
                                          optout_Cauchy_Huang = optout_Cauchy_Huang,
                                          optout_Gaussian_Huang = optout_Gaussian_Huang,
                                          optout_GenCauchy_Huang = optout_GenCauchy_Huang)
      
      
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
      save(Huang_parametric_result_list,file=paste0(result_path,'/','Huang_parametric_result_list.RData'))
      save(Huang_optout,file=paste0(result_path,'/','Huang_optout.RData'))
     
      write.table(param_est_C_zero_with_nugget_Huang,file=paste0(result_path,'/','Huang_Parametric_C_zero_with_nugget_result.txt'),col.names = FALSE,row.names = TRUE)
      write.table(param_est_nugget_Huang,file=paste0(result_path,'/','Huang_Parametric_nugget_result.txt'),col.names = FALSE,row.names = TRUE)
      write.table(param_est_vector_Huang,file=paste0(result_path,'/','Huang_Parametric_vector_result.txt'),col.names = FALSE,row.names = TRUE)
      
      
      
      
      
      
      