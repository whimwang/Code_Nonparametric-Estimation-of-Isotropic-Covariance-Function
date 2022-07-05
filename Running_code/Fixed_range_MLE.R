

library("nloptr")
library("ggplot2")
library("reshape2")#melt
library("Matrix")
library("penalized")
library(foreach)
library(doParallel)
library("base")
library("msos") #logdet function
library('Rcpp') #basis matrix
library('inline') #Rcpp basis matrix

required_packages = c('nloptr','psych','psych','RcppZiggurat','Rfast','ggplot2','reshape2','Matrix','penalized','base','msos','Rcpp','inline','geoR')

#c("nloptr", "psych","mvtnorm","RcppZiggurat","Rfast","Rsolnp","pracma","penalized",
#  "reshape2","ggplot2","base","Matrix","msos")


code.path = '/Users/wangyiming/Documents/2019 Research Daily/research_code/Code_2020/Updated0620/Code0620/'
#server
#code.path = "/mnt/home/ywang225/Cov_research/Updated0923/Code0620/"
source(paste0(code.path,"Code0912-0919/Part-I.Spectral.Cov.R"))
source(paste0(code.path,"Code0912-0919/Part-II.Data_Generation.R"))
source(paste0(code.path,"Code0912-0919/Part-III.Objective.functions.R"))
source(paste0(code.path,"Code0912-0919/Part-IIII.Check.Cov.Est_and_Plot.R"))
source(paste0(code.path,"Code1004-1010/Support.R"))

source(paste0(code.path,"Code1018-1025/LSE.w_tilde.R"))
source(paste0(code.path,"Code1018-1025/Additional_Covariance_Not_Inf_dim.R"))
#source(paste0(code.path,"Code1018-1025/Decide_Num_Basis.R"))
source(paste0(code.path,"Code1018-1025/Data_Preparation.R"))
source(paste0(code.path,"Code1018-1025/optimization.functions.R"))
source(paste0(code.path,"Code1018-1025/Plot.List.Result.R"))
source(paste0(code.path,"Code1018-1025/Write.List.Result.R"))
source(paste0(code.path,"Code1018-1025/Boxplot.all.seed.est.par.R"))
source(paste0(code.path,"Code1026-1031/optimize_with_range_w.R"))
source(paste0(code.path,"Code1026-1031/Sparse.weights.R"))
source(paste0(code.path,"Code1101-1107/Plot.List.Result.with.range.R"))
source(paste0(code.path,"Code1101-1107/Write.List.result.with.range.R"))

source(paste0(code.path,"Code0824-0830/Parametric_MLE.R"))
#source(paste0(code.path,'Code0824-0830/MLE_Berstein_with_nugget.R'))



##############################################
################# Input ######################
##############################################
d = 2; n=20;
r = 10; N = 20;
l_max = 20


seed = 0;#seed for coord generation

#### for m_seq  #####
a_index = seq(from=0.05,to=0.9,by=0.05)
m_seq = unique(1+floor((n*r)^a_index))
m_seq = m_seq[m_seq < 500] #

#option = 'Matern'; par.cov=c(1,1.25,1);metric_epsilon = 0.001; m_seq = m_seq[m_seq < 150] ##10.28
option = 'Cauchy'; par.cov = c(1,0.8); metric_epsilon = 0.05;  ## 15.9
#option = 'Gaussian'; par.cov = c(1,3); metric_epsilon = 0.001; ## 7.88
#option = 'LinearMatern'; par.cov=c(1,1.25,1,3/2/sqrt(2),2); metric_epsilon = 0.001;m_seq = m_seq[m_seq < 150] ##10.55
#option = 'GenCauchy'; par.cov=c(1,0.3,0.5,2); metric_epsilon = 0.15 ##13.32996



##### for metric grid setting #####
metric_upper_bound = Support.function(par.cov = c(1,par.cov[-1]),option = option, epsilon = metric_epsilon)
metric_grid_count = 2000
metric_grid = seq(from=0,to=metric_upper_bound,length.out = metric_grid_count)
metric_grid_space = metric_upper_bound/metric_grid_count




Fixed_l_max=TRUE; option.range = 'Fixed'
epsilon = Cov_function(l_max * 2 * sqrt(2),option,par.cov)


####################################################################
#################### Transformation and Basis ####################
####################################################################

m.m.tol = 0; m.1.tol=0
pd.adjust = TRUE;



#####################################################
###### generate seed and optimization ####################
####################################################################

n_seed = 1; doParallel = FALSE


#### setting parameters :DON NOT Change ######
if(option == 'Sph'|option == 'Wave'){
  pos.coef = FALSE
}else{
  pos.coef = TRUE; #whether estimation if positive or not
}
approx.option = 2; 



#if(option == 'Gaussian'){
#  approx.par = 10^(-6)
#}else{
#  approx.par = 10^(-9); 
#}
approx.par = 10^(-9)




################################### ################################### 
################################### Code START ##############################
################################### ################################### 



start.time = Sys.time()

################## Step 1 : Data Generation #############################

#Data.Preparation.result = Data_Preparation(option, par.cov, d, N, r=r, n=n, epsilon = epsilon, option.range = option.range, seed = 0)
#range [0,2*RQ]^d
Data.Preparation.result = Data_Preparation_Fixed_range(option,par.cov,d,N,r,n,l_max=l_max,seed = 0)

Output.list = Data.Preparation.result$Output.list
l_max = Data.Preparation.result$Hidden.list$l_max


################## Step 2 : set result path #############################

#result_dir = '/mnt/home/ywang225/Cov_research/Result_0923/Mymethod_100/'
filename = paste(c(option,round(par.cov,2),n,r,N,l_max,option.range,format(Sys.time(), "%b%d_%X")),collapse = '_')
result_path = paste0(result_dir,filename,"/")
dir.create(result_path)
save(Data.Preparation.result, file=paste0(result_path,"Data.Preparation.result.RData"))

task.name = paste0("Log_approx.option",approx.option,"_approx.par_",approx.par,
                   "_d",d,"_n",n,"_r",r,"_N",
                   N,"_n_seed",n_seed,"_option_",option,"_l_",l_max,
                   "_par.cov","(",paste0(par.cov,collapse = ","),")")


#all.best.seed.list = list()

time.total.start = Sys.time()



getwd()
cores <- detectCores(logical=T)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)




all.sample.list <- foreach(k = seq(from=1,to=N,by=1),.packages=required_packages,
                           .noexport = c('Cor.Basis.Log.rho_Rcpp_func','Log.All.Basis.Cor_Rcpp_func')) %dopar%
  {
    source(paste0(code.path,"Code0912-0919/Part-III.Objective.functions.R"))
    print(paste0("k=",k," starts!"))
    current.sample.start.time = Sys.time()
    
    one.sample.all.m.list <- list()
    
    print(k)
    Y = Output.list$Y_matrix[k,,]
    rho_matrix.orig = Output.list$rho_matrix.orig
    
    
    #Nonparamtric Method
    for(m in m_seq)                                           
    {
      w_tilde.init = rep(1/m, m)
      #w_tilde.init = runif(m, min=0,max=1)
      #w_tilde.init = w_tilde.init/sum(w_tilde.init)
      result.range.ratio = optimize_with_range_w(w_tilde.init = w_tilde.init, Y=Output.list$Y_matrix[k,,], 
                                                 rho_matrix = rho_matrix.orig, 
                                                 upper.iter = 10, upper.convergence = 10^(-3),
                                                 positive.check = TRUE, pd.check = TRUE,pd.adjust = TRUE,
                                                 approx.option, approx.par,pos.coef)
      
      
      ### Transformation #####
      #transformed.result = Transform.Distance(rho_matrix.orig, ratio = 1/range.new)# rho_matrix/ratio
      transformed.result = Transform.Distance(rho_matrix.orig, ratio = result.range.ratio$est.range)# rho_matrix/ratio
      rho_matrix = transformed.result$transformed_rho_matrix
      transform.ratio = transformed.result$ratio #max of rho_matrix.orig
      
      
      ### Calculate basis #### 
      #### Generate full rank and positive definite matrix
      Log.All.Basis.Cor_matrix = Log.Basis.Cor.Final(rho_matrix, m, positive.check = TRUE, pd.check = TRUE,pd.adjust = TRUE)
      
      V.R.result = Est.Sigma_and_Est.v_Log (est.par_matrix = matrix(result.range.ratio$est.w_tilde, nrow=1, ncol=m),
                                            Log.All.Basis.Cor_matrix,
                                            Y = Output.list$Y_matrix[k,,],
                                            approx.option=approx.option,approx.par= approx.par)
      
      result.range.ratio$est.V = V.R.result$est.V_vector
      result.range.ratio$est.R = V.R.result$est.Sigma_R_matrix
      
      ### record m
      result.range.ratio$m = m
      
      
      ### record MLE of true Sigma_R
      result.range.ratio$true_fvalue = Two.Neg.Log.Likelihood.Use.Sigma_R(Y,Sigma_R=Data.Preparation.result$Hidden.list$True.Sigma_matrix,approx.option,approx.par)
      
      
      ## to calculate metric
      f_cov_true <- function(x){
        Cov_function(x,option = option,par.cov=par.cov)
      }
      f_cor_true <- function(x){
        Cov_function(x,option = option,par.cov=c(1,par.cov[-1]))
      }
      f_cov_est <- function(x){
        tmp <- function(x){
          vector = sapply(1:m, function(k){
            return(-sum(log1p(rep(x^2/result.range.ratio$est.range^2,m-k+1)/seq(from=k,to=m,by=1))))}) 
          sum(exp(vector + log(result.range.ratio$est.w_tilde))) * result.range.ratio$est.V
        }
        sapply(x,tmp)
      }
      f_cor_est <- function(x){
        tmp <- function(x){
          vector = sapply(1:m, function(k){
            return(-sum(log1p(rep(x^2/result.range.ratio$est.range^2,m-k+1)/seq(from=k,to=m,by=1))))}) 
          sum(exp(vector + log(result.range.ratio$est.w_tilde)))
        }
        sapply(x,tmp)
      }
      
      
      
      #cov metric
      result.range.ratio$norm2_cov_integration =  try(Norm2_between_f(f_cov_true,f_cov_est,0,metric_upper_bound)/sqrt(metric_upper_bound),silent = FALSE)
      if('try-error' %in% class(result.range.ratio$norm2_cov_integration)){
        result.range.ratio$norm2_cov_integration = NA
      }
      result.range.ratio$norm2_cov_grid = sqrt(sum((f_cov_true(metric_grid) - f_cov_est(metric_grid))^2 * metric_grid_space)/metric_upper_bound) 
      result.range.ratio$supnorm_cov_grid = max(abs(f_cov_true(metric_grid) - f_cov_est(metric_grid)))
      
      #cor metric
      result.range.ratio$norm2_cor_integration =  try(Norm2_between_f(f_cor_true,f_cor_est,0,metric_upper_bound)/sqrt(metric_upper_bound),silent = FALSE)
      if('try-error' %in% class(result.range.ratio$norm2_cor_integration)){
        result.range.ratio$norm2_cor_integration = NA
      }
      result.range.ratio$norm2_cor_grid = sqrt(sum((f_cor_true(metric_grid) - f_cor_est(metric_grid))^2 * metric_grid_space)/metric_upper_bound) 
      result.range.ratio$supnorm_cor_grid = max(abs(f_cor_true(metric_grid) - f_cor_est(metric_grid)))
      
      
      ## semi-variog of norm2_variog_grid, supnorm_variog_grid, 
      f_variog_est <- function(x){
        tmp <- function(x){
          vector = sapply(1:m, function(k){
            return(-sum(log1p(rep(x^2/result.range.ratio$est.range^2,m-k+1)/seq(from=k,to=m,by=1))))}) 
          result.range.ratio$est.V - sum(exp(vector + log(result.range.ratio$est.w_tilde))) * result.range.ratio$est.V
        }
        sapply(x,tmp)
      }
      
      f_variog_true <- function(x){
        Cov_function(0,option=option,par.cov=par.cov) - Cov_function(x,option = option,par.cov=par.cov)
      }
      
      result.range.ratio$norm2_variog_integration =  try(Norm2_between_f(f_variog_true,f_variog_est,0,metric_upper_bound)/sqrt(metric_upper_bound),
                                                         silent = FALSE)
      if('try-error' %in% class(result.range.ratio$norm2_variog_integration)){
        result.range.ratio$norm2_variog_integration = NA
      }
      result.range.ratio$norm2_variog_grid = sqrt(sum((f_variog_true(metric_grid) - f_variog_est(metric_grid))^2 * metric_grid_space)/metric_upper_bound) 
      result.range.ratio$supnorm_variog_grid = max(abs(f_variog_true(metric_grid) - f_variog_est(metric_grid))) 
      
      
      one.sample.all.m.list[[length(one.sample.all.m.list)+1]] = result.range.ratio
      
    }
    
    
    
    
    ###########################################################################
    ##########################    save results      ##############################
    ###########################################################################
    #save RData of current sample's Nonparametric Estimation
    save(one.sample.all.m.list, file=paste0(result_path,task.name,"_sample_NonParam_",k,".RData"))
    
    
    current.sample.end.time = Sys.time()
    print(current.sample.end.time - current.sample.start.time)
    return(one.sample.all.m.list)
    
  }


time.total.end = Sys.time()
print(time.total.end - time.total.start)

## save Nonparametric results
rho_support = Data.Preparation.result$Hidden.list$rho_support
True.Sigma_matrix = Data.Preparation.result$Hidden.list$True.Sigma_matrix


pdf(paste0(result_path,task.name,".sample.pdf"),width = 10.5,height = 5)
layout(matrix(c(rep(1,3),2,3,4),nrow=2,ncol=3,byrow = TRUE), heights = c(1.2,3))

for(i in 1:length(m_seq)){
  ## save fig ### 
  m = all.sample.list[[1]][[i]]$m
  plot.new()
  text(0.5,0.5,paste0("m=",m),cex=2)
  #best.fig.file = paste0(task.name,"_m_",m,".jpg")
  #jpeg(best.fig.file,height = 600,width = 1800)
  
  Plot.List.Result.with.range(True.Sigma_matrix, rho_support, par.cov, option,
                              rho_matrix.orig = Data.Preparation.result$Output.list$rho_matrix.orig ,
                              lapply(1:N, function(idx){all.sample.list[[idx]][[i]]}),epsilon,
                              text.legend = "best ")
  
  
}

dev.off()


total.list = NULL
for(i in 1:N){
  total.list = c(total.list, all.sample.list[[i]])
}



### 
Write.List.result.with.range(total.list, task.name, whether.best = TRUE, text = "all",result_path=result_path)


########################################
######## Parametric Method ###############


all.sample.parametric.list <- foreach(k = seq(from=1,to=N,by=1),.packages=required_packages) %dopar%
  #for(k in seq(from=1,to=N,by=1))
  {   
    print(k)
    one.sample.all.para_est.list <- list()
    
    rho_matrix.orig = Output.list$rho_matrix.orig
    Y = Data.Preparation.result$Output.list$Y_matrix[k,,]
    
    one.sample.all.para_est.list$Matern_est = Parametric_Estimation_Metric(true_option = option,true_par.cov=par.cov,
                                                                           est_option = 'Matern',
                                                                           Y = Y, rho_matrix.orig = rho_matrix.orig,
                                                                           metric_upper_bound = metric_upper_bound, metric_grid_space = metric_grid_space, metric_grid = metric_grid)
    
    one.sample.all.para_est.list$Cauchy_est = Parametric_Estimation_Metric(true_option = option,true_par.cov=par.cov,
                                                                           est_option = 'Cauchy',
                                                                           Y = Y, rho_matrix.orig = rho_matrix.orig,
                                                                           metric_upper_bound = metric_upper_bound, metric_grid_space = metric_grid_space, metric_grid = metric_grid)
    
    one.sample.all.para_est.list$Gaussian_est = Parametric_Estimation_Metric(true_option = option,true_par.cov=par.cov,
                                                                             est_option = 'Gaussian',
                                                                             Y = Y, rho_matrix.orig = rho_matrix.orig,
                                                                             metric_upper_bound = metric_upper_bound, metric_grid_space = metric_grid_space, metric_grid = metric_grid)
    
    one.sample.all.para_est.list$GenCauchy_est = Parametric_Estimation_Metric(true_option = option,true_par.cov=par.cov,
                                                                              est_option = 'GenCauchy',
                                                                              Y = Y, rho_matrix.orig = rho_matrix.orig,
                                                                              metric_upper_bound = metric_upper_bound, metric_grid_space = metric_grid_space, metric_grid = metric_grid)
    
    
    #save Parametric Estimation
    return(one.sample.all.para_est.list)
  }
save(all.sample.parametric.list,file=paste0(result_path,task.name,"_NonParam.RData"))

## save Parametric estimate value
parametric_order = c('Matern_est','Cauchy_est','Gaussian_est','GenCauchy_est')
param_est_matrix = t(sapply(all.sample.parametric.list, function(x){
  c(x$Matern_est$par,x$Cauchy_est$par,x$Gaussian_est$par,x$GenCauchy_est$par,
    sapply(1:length(x), function(i){x[[i]]$norm2_cov_grid}),
    sapply(1:length(x), function(i){x[[i]]$norm2_cov_integration}),
    sapply(1:length(x), function(i){x[[i]]$supnorm_cov_grid}),
    sapply(1:length(x), function(i){x[[i]]$norm2_cor_grid}),
    sapply(1:length(x), function(i){x[[i]]$norm2_cor_integration}),
    sapply(1:length(x), function(i){x[[i]]$supnorm_cor_grid}),
    sapply(1:length(x), function(i){x[[i]]$norm2_variog_grid}),
    sapply(1:length(x), function(i){x[[i]]$norm2_variog_integration}),
    sapply(1:length(x), function(i){x[[i]]$supnorm_variog_grid}))
}
))
colnames(param_est_matrix) = c(paste0('Matern.',1:3),paste0('Cauchy.',1:2),paste0('Gaussian.',1:2),paste0('GenCauchy.',1:4),
                               paste0(parametric_order,'.norm2_cov_grid'),
                               paste0(parametric_order,'.norm2_cov_integration'),
                               paste0(parametric_order,'.supnorm_cov_grid'),
                               paste0(parametric_order,'.norm2_cor_grid'),
                               paste0(parametric_order,'.norm2_cor_integration'),
                               paste0(parametric_order,'.supnorm_cor_grid'),
                               paste0(parametric_order,'.norm2_variog_grid'),
                               paste0(parametric_order,'.norm2_variog_integration'),
                               paste0(parametric_order,'.supnorm_variog_grid'))

write.table(param_est_matrix,file=paste0(result_path,'Parametric_est_result.txt'),row.names = FALSE)



stopImplicitCluster()
stopCluster(cl)



end.time = Sys.time()
print(paste0("total time is ",end.time - start.time))




pdf(file = paste0(result_path,task.name,".summary.pdf"))

layout(matrix(1:4,2,2,byrow = TRUE))

# plot cor
for(i in 1:length(m_seq)){
  plot(NA,xlim=c(0,max(Data.Preparation.result$Hidden.list$True.Sigma_matrix)*1.2), ylim=c(0,1+0.2),
       main=paste0("m=", m_seq[i]),ylab="est cor",xlab="true cor")
  for(k in 1:N){
    #print(k)
    true.cor = as.matrix(Data.Preparation.result$Hidden.list$True.Sigma_matrix)/Data.Preparation.result$Hidden.list$True.Sigma_matrix[1,1]
    order.true.cor = order(true.cor,decreasing=FALSE)
    points(true.cor[order.true.cor], as.matrix(all.sample.list[[k]][[i]]$est.R)[order.true.cor],pch=16,cex=0.5)
    lines(true.cor[order.true.cor], as.matrix(all.sample.list[[k]][[i]]$est.R)[order.true.cor])
  }
  abline(a=0,b=1,col="red")
  legend("topright",
         legend=c("y=x",paste0("Sim:",1:N)),
         col=c("red",rep("black",N)),
         lty=rep(1,N+1),bty="n")
  
}

layout(matrix(1))

est.V.matrix = NULL
for(i in 1:length(m_seq)){
  est.V.vector = sapply(all.sample.list, function(x){x[[i]]$est.V})
  est.V.matrix = cbind(est.V.matrix,est.V.vector)
}
colnames(est.V.matrix)=paste0("m=", m_seq)
boxplot(est.V.matrix, main="boxplot of est V",use.cols = TRUE)
abline(h=1,col="red")

dev.off()


pdf(paste0(result_path,task.name,"_boxplot_par_range.pdf"))
layout(matrix(1:2,2,1),height=c(2,1.5))
final.par.matrix = NULL
for(i in 1:length(m_seq)){
  final.par.matrix = t(sapply(all.sample.list, function(x){x[[i]]$est.w_tilde}))
  colnames(final.par.matrix)=paste0("m=", 1:m_seq[i])
  boxplot(final.par.matrix, use.cols = TRUE, main=paste0("boxplot of final w_tilde when m=",m_seq[i]))
  
  
  range.vector = sapply(all.sample.list, function(x){x[[i]]$est.range})
  boxplot(range.vector, main="boxplot of est range", horizontal = TRUE)
}



dev.off()


save.image(paste0(result_path,task.name,".RData"))



