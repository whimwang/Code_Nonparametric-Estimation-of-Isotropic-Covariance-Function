#Approximated Error
#2020 April 5: code from Code0920-0926/Solve.QP.Simulation.R

#1. Gaussian: m = 5,10,15,20




#function.path = "/mnt/home/ywang225/research_code/"
#function.path = "/Users/wangyiming/Documents/2019 Research Daily/2019Septemper/code0912-0919/"

#My Mac 
#function.path = "/Users/wangyiming/Documents/2019 Research Daily/research_code/Code1122-1129/"
function.path = '/Users/wangyiming/Documents/2019 Research Daily/research_code/Code_2020/Updated0620/Code0620/'

#server
#function.path = "/home/whimwang/research_code/Upload0425/Code1122-1129/"

source(paste0(function.path,"Code0912-0919/Part-I.Spectral.Cov.R"))
source(paste0(function.path,"Code0912-0919/Part-II.Data_Generation.R"))
source(paste0(function.path,"Code0912-0919/Part-III.Objective.functions.R"))
source(paste0(function.path,"Code0912-0919/Part-IIII.Check.Cov.Est_and_Plot.R"))
source(paste0(function.path,"Code1004-1010/Support.R"))


#source(paste0(function.path,"Code0920-0926/Function.Norm.R"))
#source(paste0(function.path,"Code0927-1003/RQ.Norm.Diff.R"))
source(paste0(function.path,"Code0927-1003/Update_Approximation_Error.R"))

source(paste0(function.path,"Code1018-1025/Additional_Covariance_Not_Inf_dim.R"))


library("quadprog")
library("Matrix")
library("geoR") #matern function in geoR


library(foreach)
library(doParallel)
library(Rsolnp)
#library(erer)


#option = "GneitMatern"; theta_2 = c(1,1,1); 
#option = 'Circular'; theta_2 = c(1,0.5,2)
#phi,kappa for Matern, where phi is range parameter
#phi,kappa1,kappa2 for GenCauchy, where phi is range parameter 
#option = 'Gaussian'; theta_2 =1
#option = 'Exp'; theta_2 = 1
option = 'Cauchy'; theta_2=1

if(option == 'Matern'){
  task.name = paste0("option_", option, "_phi_",theta_2[1],"_kappa_",theta_2[2])
}else if(option == 'GenCauchy'|option == 'GneitMatern'){
  task.name = paste0("option_", option, "_phi_",theta_2[1],"_kappa_",theta_2[2],"_",theta_2[3])
}else{
  task.name = paste0("option_", option, "_theta2_",theta_2)
}




#m_seq = c(seq(from=5,to=20, by=5),seq(from=40,to=80,by=20))
#m_seq = c(seq(from=1,to=20, by=1))
m_seq = c(seq(from=100,to=1000, by=100))
m_vector <- time.record.vector <- Diff.est.vector <- Norm2.vector <- numeric(length(m_seq))

sqrt(pi/2)*0.05


#m=150,Cauchy 1: norm2: 0.08304153



cores <- detectCores(logical=F)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)


#for(i in 1:length(m_seq)){
all.m_seq.list <- foreach(i = 1:length(m_seq),.packages=c("quadprog","Matrix","Rfast","geoR")) %dopar%
  {
    
    start.time = Sys.time()
    m = m_seq[i]
    
    Diff.matrix = Diff.Matrix.func(m, option, theta_2)
    #when m=30, not numerically positive definite
    PD.Diff.matrix = as.matrix(nearPD(Diff.matrix)$mat)
    
    #print(paste0("eigen is ",paste0(round(eigen(Diff.matrix)$value,4),collapse = " ")))
    
    opt.result = solve.QP(Dmat = PD.Diff.matrix, dvec = numeric(m),
                          Amat = t(rbind(rep(1,m),diag(rep(1,m),m))),
                          bvec = c(1, rep(0,m)),
                          meq = 1)#first 1 constraint is equality
    
    end.time = Sys.time()
    time.record.vector[i] = difftime(end.time, start.time, units = "mins")
    
    
    #print(paste0("m: ",m, " obj.value: ", Diff.est.vector[i]))
    #print(paste0("Time:", time.record[i]))
    
    
    current.seed.list <- list(m = m,
                              time = difftime(end.time, start.time, units = "mins"),
                              est.par = opt.result$solution,
                              est.fvalue = opt.result$value,
                              norm2 = sqrt(opt.result$value))
    
  }



#list to matrix
for(i in 1:length(m_seq)){
  m_vector[i] = all.m_seq.list[[i]]$m
  time.record.vector[i] = all.m_seq.list[[i]]$time
  Diff.est.vector[i] = all.m_seq.list[[i]]$est.fvalue
  Norm2.vector[i] = all.m_seq.list[[i]]$norm2
}

result.summary = data.frame(m_vector, Diff.est.vector,time.record.vector, Norm2.vector)

write.table(result.summary, file = paste0(task.name,"QP.result.txt"), quote = FALSE)


#write out pars
est.par.matrix = matrix(NA,nrow=length(m_seq),ncol=max(m_seq))
rownames(est.par.matrix) = paste0("m=",m_seq)
for(i in 1:length(m_seq)){
  est.par.matrix[i,1:m_seq[i]] = all.m_seq.list[[i]]$est.par 
}
write.table(est.par.matrix, file = paste0(task.name,"QP.est.par.result.txt"),
            col.names = FALSE,na="NA",quote = FALSE)


save.image(file=paste0(task.name, ".QP.RData"))


stopImplicitCluster()
stopCluster(cl)






jpeg(paste0(task.name,".jpg"))
Plot.Result.Check.Diff.Norm(all.m_seq.list, option, theta_2)
dev.off()


Plot.Result.Check.Diff.Norm(all.m_seq.list, option, theta_2)


