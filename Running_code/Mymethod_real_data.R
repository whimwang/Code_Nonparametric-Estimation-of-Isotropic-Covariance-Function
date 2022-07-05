#each time: submit a large m_seq to allow to run all the time: set do.parallel=TRUE and meanwhile run a small m to get parametric result : set do.parallel=FALSE

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
code.path = '/Users/wangyiming/Documents/2019 Research Daily/research_code/Code_2020/Updated0620/Code0620/'
#server
#code.path = "/mnt/home/ywang225/Cov_research/Updated1022/Code0620/"
source(paste0(code.path,"Code0912-0919/Part-I.Spectral.Cov.R"))
source(paste0(code.path,"Code0912-0919/Part-II.Data_Generation.R"))
source(paste0(code.path,"Code0912-0919/Part-III.Objective.functions.R"))
source(paste0(code.path,"Code1018-1025/optimization.functions.R"))
source(paste0(code.path,'Code0824-0830/MLE_Berstein_with_nugget.R'))
source(paste0(code.path,'Code0824-0830/Parametric_MLE_with_nugget.R'))


list_to_arg_vector <- function(arg_list){
  as.vector(sapply(1:(length(arg_list)), function(i){c(names(arg_list[i]),as.character(arg_list[i]))}))
}

Mymethod_arg = list(cov.distort = 4^2,dist.range_parameter = 1*(10^5))
do.parallel=FALSE

#Server
result_dir = '/mnt/home/ywang225/Cov_research/Result_1022/Mymethod_real_data/'
#result_dir = './'
filename = paste(c(list_to_arg_vector(Mymethod_arg),format(Sys.time(), "%b%d_%X")),collapse = '_')
result_path = paste0(result_dir,filename)
dir.create(result_path)

input_data.path = "/Users/wangyiming/Documents/2019 Research Daily/research_code/Code_2020/Updated0620/Code0620/input_real_data/Year40_Station189_46_30_80_100/"
#Server
#input_data.path = '/mnt/home/ywang225/Cov_research/Updated1016/Code0620/input_real_data/Year40_Station189_46_30_80_100/'
load(paste0(input_data.path,"Year40_Station189_46_30_80_100.RData"))

lon_lat_data = input_data$lon_lat_data
x_y_scale = input_data$x_y_scale
observation_residual = input_data$observation_residual
dist_select = input_data$dist_select

n = ncol(observation_residual)
r = nrow(observation_residual)
print(paste0('n:  ',n,'  r:',r))
idx = 1:n #241

cov_with_nugget = mean(apply(observation_residual[,idx],2,var))
#cov_with_nugget = mean(apply(observation_residual[,idx],2,function(x){mean(x^2)}))

cov.distort = Mymethod_arg$cov.distort
dist.range_parameter = Mymethod_arg$dist.range_parameter

a_index = seq(from=0.05,to=0.9,by=0.05)
m_seq = unique(1+floor((n*r)^a_index))
m_seq = m_seq[m_seq < 100] 
time_seq = numeric(length(m_seq))


if(do.parallel==TRUE){
  getwd()
  cores <- detectCores(logical=T)
  cl <- makeCluster(cores)
  registerDoParallel(cl, cores=cores)

  ### NP myMethod #####
  
  #all.m_seq.result = list()
  #for(i in 1:length(m_seq)){
  all.m_seq.result <- foreach(i = seq(from=1,to=length(m_seq),by=1),.packages=required_packages,
                             .noexport = c('Cor.Basis.Log.rho_Rcpp_func','Log.All.Basis.Cor_Rcpp_func')) %dopar%
    {
    source(paste0(code.path,"Code0912-0919/Part-III.Objective.functions.R"))
    print(paste0("i=",i," starts!"))
    current.m.start.time = Sys.time()
    
    m = m_seq[i]
    par.init = c(0.5,rep(1/m,times=m))
    
    rho_matrix.change = dist_select[idx,][,idx]/dist.range_parameter
    opt = optimize_with_range_w_nugget(par.init, Y=observation_residual[,idx]/sqrt(cov.distort), 
                                       rho_matrix= rho_matrix.change, 
                                       upper.iter = 10, upper.convergence = 10^(-3),
                                       positive.check = TRUE, pd.check = TRUE, pd.adjust = TRUE,
                                       approx.option=2, approx.par=10^(-9), pos.coef=TRUE)
    opt$m = m
    opt$cov.distort = cov.distort
    opt$dist.range_parameter = dist.range_parameter
    
    current.m.end.time = Sys.time()
    opt$time = difftime(current.m.end.time,current.m.start.time,units = 'hours')
   
    
    
    #save results of current m
    opt_vector = c(opt$m,opt$time,opt$dist.range_parameter, opt$cov.distort,opt$est.range,opt$est.fvalue,opt$est.par)
    names(opt_vector) = c('m','time','dist.range_parameter','cov.distort','est.range','est.fvalue','nugget',paste0('coef',1:opt$m))
    write.table(opt_vector,file = paste0(result_path,'/Mymethod_m_',m,'.txt'),col.names = FALSE)
    
    time_seq[i] = opt$time
    print(paste0("##############  m=",m,"     #################"))
    print(current.m.end.time - current.m.start.time)
    
    #all.m_seq.result[[i]] = opt
    return(opt)
    }
    stopImplicitCluster()
    stopCluster(cl)
  
}else{
  ### NP myMethod #####
  all.m_seq.result=list()
  for(i in 1:length(m_seq)){
      print(paste0("i=",i," starts!"))
      current.m.start.time = Sys.time()
      
      m = m_seq[i]
      par.init = c(0.5,rep(1/m,times=m))
      
      rho_matrix.change = dist_select[idx,][,idx]/dist.range_parameter
      opt = optimize_with_range_w_nugget(par.init, Y=observation_residual[,idx]/sqrt(cov.distort), 
                                         rho_matrix= rho_matrix.change, 
                                         upper.iter = 10, upper.convergence = 10^(-3),
                                         positive.check = TRUE, pd.check = TRUE, pd.adjust = TRUE,
                                         approx.option=2, approx.par=10^(-9), pos.coef=TRUE)
      opt$m = m
      opt$cov.distort = cov.distort
      opt$dist.range_parameter = dist.range_parameter
      
      current.m.end.time = Sys.time()
      opt$time = difftime(current.m.end.time,current.m.start.time,units = 'hours')
      
      #save results of current m
      opt_vector = c(opt$m,opt$time,opt$dist.range_parameter, opt$cov.distort,opt$est.range,opt$est.fvalue,opt$est.par)
      names(opt_vector) = c('m','time','dist.range_parameter','cov.distort','est.range','est.fvalue','nugget',paste0('coef',1:opt$m))
      write.table(opt_vector,file = paste0(result_path,'/Mymethod_m_',m,'.txt'),col.names = FALSE)
      
      time_seq[i] = opt$time
      print(paste0("##############  m=",m,"     #################"))
      print(current.m.end.time - current.m.start.time)
      print(paste0('fvalue:  ',opt$est.fvalue))
      all.m_seq.result[[i]] = opt
    }
  names(time_seq) = m_seq
  print(time_seq)
}






### Find the best m #####
if(TRUE){

Mymethod_all.est_fvalue = sapply(all.m_seq.result,function(x){x$est.fvalue})
Mymethod_all.est_nugget = sapply(all.m_seq.result,function(x){x$est.par[1]})
Mymethod_all.est_range = sapply(all.m_seq.result,function(x){x$est.range})
Mymethod_all.sum = sapply(all.m_seq.result,function(x){sum(x$est.par)})
Mymethod_all.dist.range_parameter = sapply(all.m_seq.result,function(x){x$dist.range_parameter})
Mymethod_all.cov.distort = sapply(all.m_seq.result,function(x){x$cov.distort})

fvalue_ratio = Mymethod_all.est_fvalue[2:length(Mymethod_all.est_fvalue)]/Mymethod_all.est_fvalue[1:(length(Mymethod_all.est_fvalue)-1)] - 1
names(fvalue_ratio) = paste0(m_seq[2:length(m_seq)],":",m_seq[1:(length(m_seq)-1)])


#relative ratio
if(Mymethod_all.est_fvalue[length(m_seq)] > 0){
  tmp = which(fvalue_ratio <0 & fvalue_ratio >- 1*10^(-3))
  if(length(tmp)!=0){
    index2 = min(which(fvalue_ratio <0 & fvalue_ratio >- 1*10^(-3)))+1
  }else{
    index2 = NA
  }
}else{
  tmp = which(fvalue_ratio>0 & fvalue_ratio  < 1*10^(-3))
  if(length(tmp)!=0){
    index2 = min(which(fvalue_ratio>0 & fvalue_ratio  < 1*10^(-3)))+1
  }else{
    index2 = NA
  }
}

print(index2)

#which min
index0 =  which.min(Mymethod_all.est_fvalue)

#which min not the last
index1 = ifelse(index0 == length(m_seq),NA,index0)

#if min is the last one: find the relative ratio one and if no relative ratio then the last one
#if min is the middle one; accept the middle one
index_combined = ifelse(!is.na(index1),index1,index2)
index_combined = ifelse(is.na(index_combined),index0,index_combined)
final_result = all.m_seq.result[[index_combined]]

print(paste0('index0: ',index0,'  index0 m:',m_seq[index0]))
print(paste0('index1: ',index1,'  index1 m:',m_seq[index1]))
print(paste0('index2: ',index2,'  index2 m:',m_seq[index2]))
print(paste0('best index: ',index_combined,'        best m:',m_seq[index_combined]))

}



### Parametric #####
time_Parametric_seq = rep(NA,4)

start.time = Sys.time()
optout_Matern <-  try(MLE_est_Matern_with_nugget(par.init=c(1,1,1),Y=observation_residual[,idx]/sqrt(cov.distort),
                                                 rho_matrix=dist_select[idx,][,idx]/dist.range_parameter,theta_3 = NULL),
                      silent = FALSE)
if('try-error' %in% class(optout_Matern)){
  optout_Matern = list(par = NA)
}
end.time = Sys.time()
time_Parametric_seq[1] = end.time - start.time

start.time = Sys.time()
optout_Cauchy <-  try(MLE_est_Cauchy_with_nugget(par.init=c(1,1,1),Y=observation_residual[,idx]/sqrt(cov.distort),rho_matrix=dist_select[idx,][,idx]/dist.range_parameter),
                   silent = FALSE)
if('try-error' %in% class(optout_Cauchy)){
  optout_Cauchy = list(par = NA)
}
end.time = Sys.time()
time_Parametric_seq[2] = end.time - start.time

start.time = Sys.time()
optout_Gaussian <-  try(MLE_est_Gaussian_with_nugget(par.init=c(1,1,1),Y=observation_residual[,idx]/sqrt(cov.distort),rho_matrix=dist_select[idx,][,idx]/dist.range_parameter),
                   silent = FALSE)
if('try-error' %in% class(optout_Gaussian)){
  optout_Gaussian = list(par = NA)
}
end.time = Sys.time()
time_Parametric_seq[3] = end.time - start.time

start.time = Sys.time()
optout_GenCauchy <-  try(MLE_est_GenCauchy_with_nugget(par.init=c(1,1,1),Y=observation_residual[,idx]/sqrt(cov.distort),rho_matrix=dist_select[idx,][,idx]/dist.range_parameter),
                        silent = FALSE)
if('try-error' %in% class(optout_GenCauchy)){
  optout_GenCauchy = list(par = NA)
}
end.time = Sys.time()
time_Parametric_seq[4] = end.time - start.time

## save par est
param_est_option = c('Matern','Cauchy','Gaussian','GenCauchy')
param_est_nugget = rep(NA,4)
param_est_vector = rep(NA,17)
param_est_vector[c(1,2)] = c(dist.range_parameter,cov.distort)
if(!is.na(optout_Matern$par)){
  param_est_vector[3:6] = optout_Matern$par
  param_est_nugget[1] = optout_Matern$par[1]
}
if(!is.na(optout_Cauchy$par)){
  param_est_vector[7:9] = optout_Cauchy$par
  param_est_nugget[2] = optout_Cauchy$par[1]
}
if(!is.na(optout_Gaussian$par)){
  param_est_vector[10:12] = optout_Gaussian$par
  param_est_nugget[3] = optout_Gaussian$par[1]
}
if(!is.na(optout_GenCauchy$par)){
  param_est_vector[13:17] = optout_GenCauchy$par
  param_est_nugget[4] = optout_GenCauchy$par[1]
}

names(param_est_vector) = c('dist.range_parameter','cov.distort',
  paste0('Matern.',1:4),paste0('Cauchy.',1:3),paste0('Gaussian.',1:3),paste0('GenCauchy.',1:5))
names(param_est_nugget) = param_est_option
names(time_Parametric_seq) = param_est_option
print(time_Parametric_seq)





######## Plot all NP #########
if(TRUE){
pdf(file=paste0(result_path,'/Mymethod_',format(Sys.time(), "%b%d_%X"),'_plot.pdf'))

#layout(matrix(1:2,1,2))
d_seq = seq(from=0,to=max(dist_select)*1.2,length.out = 100)

plot(NA,ylim=c(0,max(Mymethod_all.sum))*1.2*cov.distort,xlim=c(0,max(d_seq))/(10^4),xlab='dist/10^4',ylab='Cov',main="Cov without nugget:Mymethod")      
all_color  = 1:length(m_seq)
for(i in 1:length(m_seq)){
  opt = all.m_seq.result[[i]]
  m = length(opt$est.par[-1])
  
  f_cov_est_Mymethod <- function(x){
    tmp <- function(x){
      vector = sapply(1:m, function(k){
        return(-sum(log1p(rep(x^2/(dist.range_parameter*opt$est.range)^2,m-k+1)/seq(from=k,to=m,by=1))))}) 
      sum(exp(vector + log(opt$est.par[-1]))) 
    }
    sapply(x,tmp)*cov.distort
  }
  
  lines(d_seq/(10^4),sapply(d_seq, f_cov_est_Mymethod),col=all_color[i],lwd=2,lty='dashed') 
}
abline(v=max(dist_select)/10^4,col='red',lty='dashed')
legend('topright',legend = paste0('m=',m_seq,' nugget=',round(Mymethod_all.est_nugget*100,digits = 4)),bty='n',
       lty=rep('dashed',length(m_seq)),col=all_color)




plot(NA,ylim=c(0,max(Mymethod_all.sum))*1.2*cov.distort,xlim=c(0,max(d_seq))/10^4,main='semivariogram with nugget',xlab='dist/10^4',ylab='semi-variogram')      

all_color  = 1:length(m_seq) 
for(i in 1:length(m_seq)){
  opt = all.m_seq.result[[i]]
  m = length(opt$est.par[-1])
  
  f_cov_est_Mymethod <- function(x){
    tmp <- function(x){
      vector = sapply(1:m, function(k){
        return(-sum(log1p(rep(x^2/(dist.range_parameter*opt$est.range)^2,m-k+1)/seq(from=k,to=m,by=1))))}) 
      sum(exp(vector + log(opt$est.par[-1]))) 
    }
    sapply(x,tmp) * cov.distort
  }
  
  f_variog_est_Mymethod <- function(x){
    
    sapply(x,function(x){
      if(x==0){
        return(opt$est.par[1] * cov.distort)
      }else{
        cov.distort*sum(opt$est.par[-1])-f_cov_est_Mymethod(x) + opt$est.par[1]*cov.distort
      }
    })
    
  }
  lines(d_seq/(10^4),sapply(d_seq, f_variog_est_Mymethod),col=all_color[i],lwd=2,lty='dashed') 
}
abline(v=max(dist_select)/10^4,col='red',lty='dashed')
legend('topright',legend = paste0('m=',m_seq,' nugget=',round(Mymethod_all.est_nugget*cov.distort,digits = 3)),lty=rep('solid',length(m_seq)),col=all_color,
       bty='n')

dev.off()
}



#### Plot Parametric #####
if(TRUE){
  
pdf(file=paste0(result_path,'/Parametric_vs_best_NP_',format(Sys.time(), "%b%d_%X"),'_plot.pdf'))
## cov without nugget
plot(NA,ylim=c(0,max(Mymethod_all.sum)*cov.distort)*1.2,xlim=c(0,max(d_seq))/10^4,main='Cov of Parametric vs best Mymethod NP without nugget',xlab='dist/10^4',ylab='cov')      

if(optout_Matern$par){
  lines(d_seq/10^4, sapply(d_seq/dist.range_parameter, Cov_function,option='Matern',par.cov=optout_Matern$par[-1])*cov.distort,
        lty='dashed',col=1,lwd=2)
}
if(optout_Cauchy$par){
  lines(d_seq/10^4, sapply(d_seq/dist.range_parameter, Cov_function,option='Cauchy',par.cov=optout_Cauchy$par[-1])*cov.distort,
        lty='dashed',col=2,lwd=2)
}
if(optout_Gaussian$par){
  lines(d_seq/10^4, sapply(d_seq/dist.range_parameter, Cov_function,option='Gaussian',par.cov=optout_Gaussian$par[-1])*cov.distort,
        lty='dashed',col=3,lwd=2)
}
if(optout_GenCauchy$par){
  lines(d_seq/10^4, sapply(d_seq/dist.range_parameter, Cov_function,option='GenCauchy',par.cov=optout_GenCauchy$par[-1])*cov.distort,
        lty='dashed',col=4,lwd=2)
}

#best NP
opt = final_result
m = length(opt$est.par[-1])

best_cov_est_Mymethod <- function(x){
  tmp <- function(x){
    vector = sapply(1:m, function(k){
      return(-sum(log1p(rep(x^2/(dist.range_parameter*opt$est.range)^2,m-k+1)/seq(from=k,to=m,by=1))))}) 
    sum(exp(vector + log(opt$est.par[-1]))) 
  }
  sapply(x,tmp)*cov.distort
}

lines(d_seq/(10^4),sapply(d_seq, best_cov_est_Mymethod),col=5,lwd=2,lty='dashed') 
abline(v=max(dist_select)/10^4,col='red',lty='dashed')
legend('topright',legend = c(paste0(param_est_option,'   nugget:',round(param_est_nugget*cov.distort,digits = 3)),
                             paste0('best_NP m=',m_seq[index_combined],' nugget:',round(Mymethod_all.est_nugget[index_combined]*cov.distort,digits=3))),
       col = 1:5,lty=rep('dashed',5),bty='n')



## variogram with nugget
plot(NA,ylim=c(0,max(Mymethod_all.sum)*cov.distort)*1.2,xlim=c(0,max(d_seq))/10^4,main='SemiVariog of Parametric vs best Mymethod NP with nugget',xlab='dist/10^4',ylab='semivariogram')      

if(optout_Matern$par){
  lines(d_seq/10^4, 
        (optout_Matern$par[1] + optout_Matern$par[2] - 
           sapply(d_seq/dist.range_parameter, Cov_function,option='Matern',par.cov=optout_Matern$par[-1]))*cov.distort,
        lty='dashed',col=1,lwd=2)
}
if(optout_Cauchy$par){
  lines(d_seq/10^4, 
        (optout_Cauchy$par[1] + optout_Cauchy$par[2] - 
           sapply(d_seq/dist.range_parameter, Cov_function,option='Cauchy',par.cov=optout_Cauchy$par[-1]))*cov.distort,
        lty='dashed',col=2,lwd=2)
}
if(optout_Gaussian$par){
  lines(d_seq/10^4, 
        (optout_Gaussian$par[1] + optout_Gaussian$par[2] - 
           sapply(d_seq/dist.range_parameter, Cov_function,option='Gaussian',par.cov=optout_Gaussian$par[-1]))*cov.distort,
        lty='dashed',col=3,lwd=2)
}
if(optout_GenCauchy$par){
  lines(d_seq/10^4, 
        (optout_GenCauchy$par[1] + optout_GenCauchy$par[2] - 
           sapply(d_seq/dist.range_parameter, Cov_function,option='GenCauchy',par.cov=optout_GenCauchy$par[-1]))*cov.distort,
        lty='dashed',col=4,lwd=2)
}

best_variog_est_Mymethod <- function(x){
  sapply(x,function(x){
    if(x==0){
      return(opt$est.par[1] * cov.distort)
    }else{
      cov.distort*sum(opt$est.par[-1])-best_cov_est_Mymethod(x) + opt$est.par[1]*cov.distort
    }
  })
}

lines(d_seq/(10^4),sapply(d_seq, best_variog_est_Mymethod),col=5,lwd=2,lty='dashed') 

legend('topright',legend = c(paste0(param_est_option,'   nugget:',round(param_est_nugget*cov.distort,digits = 3)),
                             paste0('best_NP m=',m_seq[index_combined],' nugget:',round(Mymethod_all.est_nugget[index_combined]*cov.distort,digits=3))),
       col = 1:5,lty=rep('dashed',5),bty='n')


dev.off()
}



####### save results ##############
write_Mymethod_all_result_func <- function(all.m_seq.result,filename){
  
  Mymethod_all.est_fvalue = sapply(all.m_seq.result,function(x){x$est.fvalue})
  Mymethod_all.est_nugget = sapply(all.m_seq.result,function(x){x$est.par[1]})
  Mymethod_all.est_range = sapply(all.m_seq.result,function(x){x$est.range})
  Mymethod_all.sum = sapply(all.m_seq.result,function(x){sum(x$est.par)})
  Mymethod_all.dist.range_parameter = sapply(all.m_seq.result,function(x){x$dist.range_parameter})
  Mymethod_all.cov.distort = sapply(all.m_seq.result,function(x){x$cov.distort})
  
  names(Mymethod_all.est_fvalue) = m_seq
  names(Mymethod_all.est_nugget) = m_seq
  names(Mymethod_all.est_range) = m_seq
  names(Mymethod_all.sum) = m_seq
  names(Mymethod_all.dist.range_parameter) = m_seq
  names(Mymethod_all.cov.distort) = m_seq
  
  write.table(x = Mymethod_all.est_fvalue,file = paste0(filename,".est_fvalue.txt"),col.names = FALSE,row.names = TRUE)
  write.table(x = Mymethod_all.est_nugget,file = paste0(filename,".est_nugget.txt"),col.names = FALSE,row.names = TRUE)
  write.table(x = Mymethod_all.est_range,file = paste0(filename,".est_range.txt"),col.names = FALSE,row.names = TRUE)
  write.table(x = Mymethod_all.sum, file = paste0(filename,".sum.txt"),col.names = FALSE,row.names = TRUE)
  write.table(x = Mymethod_all.dist.range_parameter, file = paste0(filename,".dist.range_parameter.txt"),col.names = FALSE,row.names = TRUE)
  write.table(x = Mymethod_all.cov.distort, file = paste0(filename,".cov.distort.txt"),col.names = FALSE,row.names = TRUE)
  
}

write_Mymethod_all_result_func(all.m_seq.result, filename = paste0(result_path,'/','Mymethod_all_result'))
write.table(x = time_seq, file = paste0(result_path,"/","Mymethod_time_seq.txt"),col.names = FALSE,row.names = TRUE)


Parametric_result_list = list(optout_Matern = optout_Matern, optout_Cauchy = optout_Cauchy,
                              optout_Gaussian = optout_Gaussian, optout_GenCauchy = optout_GenCauchy)

save(Parametric_result_list,file=paste0(result_path,'/','Parametric_result_list.RData'))
save(all.m_seq.result,file=paste0(result_path,'/','Mymethod_all.m_seq.result.RData'))

write.table(param_est_vector, file=paste0(result_path,'/','Parametric_est_vector.txt'),col.names = FALSE,row.names = TRUE)
write.table(param_est_nugget, file=paste0(result_path,'/','Parametric_est_nugget.txt'),col.names = FALSE,row.names = TRUE)
write.table(time_Parametric_seq, file=paste0(result_path,'/','Parametric_est_time.txt'),col.names = FALSE,row.names = TRUE)


##### save NP results #####

#all m
m.vector = sapply(all.m_seq.result, function(x){return(x$m)})
#final range
transform.ratio.vector = sapply(all.m_seq.result, function(x){return(x$est.range)})
names(transform.ratio.vector) = paste("m=",m.vector)
#final w_tilde
est.par_matrix = matrix(NA, nrow = length(m.vector), ncol=max(m.vector)+1)
for(k in 1:length(m_seq)){
  est.par_matrix[k,1:(m.vector[k]+1)] = all.m_seq.result[[k]]$est.par
}
rownames(est.par_matrix)=paste("m=",m.vector)
colnames(est.par_matrix)=c('nugget',paste0('coef',1:max(m.vector)))

write.table(transform.ratio.vector, file=paste0(result_path,'/','Mymethod_est_transformed_ratio.txt'),col.names = FALSE,row.names = TRUE)
write.table(est.par_matrix, file=paste0(result_path,'/','Mymethod_est_par_matrix.txt'),col.names = TRUE,row.names = TRUE)

