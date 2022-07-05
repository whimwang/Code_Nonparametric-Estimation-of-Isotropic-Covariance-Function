
#'This function generate positive definte and full rank basis
#'@return: Log.All.Basis.Cor_matrix: (n,n,m)

Log.Basis.Cor.Final <- function(rho_matrix, m, positive.check = TRUE, pd.check = TRUE,pd.adjust = TRUE){
  
  Log.All.Basis.Cor_matrix = Log.All.Basis.Cor(rho_matrix, m, 
                                               positive.check = positive.check, 
                                               pd.check = pd.check, pd.adjust=pd.adjust)
  
  All.Basis.Cor_matrix = exp(Log.All.Basis.Cor_matrix)
  
  #check positive definite after adjustment
  if(FALSE){
  if(pd.adjust){
    if(PD.Check.3d(All.Basis.Cor_matrix)){
      #print("After adjustment, pd")
    }else{
      #print("After adjustment, not pd")
    }
  }
  
  #check rank:do not check rank, updated on 20200621
  if(FALSE){
    rank.vector = numeric(m)
    for(t in seq(from=1,to=m,by=1)){
      rank.vector[t] = qr(All.Basis.Cor_matrix[,,t],LAPACK = TRUE)$rank # not set FALSE
    }
  }
  #print(rank.vector)
  #rank.file = paste0(task.name,"rank.txt")
  #write.table(rank.vector, file=rank.file,row.names = FALSE,quote=FALSE,col.names = FALSE)
  
  #if(min(rank.vector) != n){
  #  print("Not full rank!")
  #}else{
  #  print("fUll rank!")
  #}
  }
  return(Log.All.Basis.Cor_matrix)
}












#' This function return the best est par
#' @return a list, best.seed.list and all.seed.list



Optimization.Parallel <- function(Y, Log.All.Basis.Cor_matrix,n_seed = 2,doParallel,
                                  all.seed.list.return = TRUE,
                                  pos.coef = TRUE, approx.option = 2, approx.par = 10^(-8)){
  
  n_seed = n_seed;
  m = dim(Log.All.Basis.Cor_matrix)[3]
  
  #### setting parameters :DON NOT Change ######
  pos.coef = pos.coef; #whether estimation if positive or not
  approx.option = approx.option; approx.par = approx.par; 
  
  if(doParallel){
    
    cores <- detectCores(logical=T)
    cl <- makeCluster(cores)
    registerDoParallel(cl, cores=cores)
    
    
    all.seed.list <- foreach(i = 1:n_seed,
                             .packages=c("nloptr", "psych","mvtnorm","RcppZiggurat","Rfast","Rsolnp","pracma","penalized")) %dopar%
      {
        
        function.path = "//wolftech.ad.ncsu.edu/cos/stat/Redirect/ywang225/Desktop/Exp_lambda_basis/2019September/Code1026-1031/"
        
        source(paste0(function.path,"Code0912-0919/Part-I.Spectral.Cov.R"))
        source(paste0(function.path,"Code0912-0919/Part-II.Data_Generation.R"))
        source(paste0(function.path,"Code0912-0919/Part-III.Objective.functions.R"))
        source(paste0(function.path,"Code0912-0919/Part-IIII.Check.Cov.Est_and_Plot.R"))
        source(paste0(function.path,"Code1004-1010/Support.R"))
      
        
        set.seed(i*100)
        
        #if(i == 1){
        #  init_par = LSE.w_tilde(Y, Log.All.Basis.Cor_matrix, pos.coef)
        
        #}else{
        init_par = runif(m,0,1)
        init_par = init_par/sum(init_par)
        #}
        #print(init_par)
        if(pos.coef){
          sol = solnp(pars=init_par,fun = Two.Neg.Log.Likelihood.Use.Basis.Log,
                      eqfun = Equal.Constrain.Log, eqB = 0,
                      ineqfun = NULL, ineqLB = NULL, ineqUB = NULL,
                      LB = rep(0,m), UB = rep(1,m), control= NULL,Y = Y,
                      Log.All.Basis.Cor_matrix = Log.All.Basis.Cor_matrix,
                      approx.option = approx.option, approx.par = approx.par)
        }else{
          sol = solnp(pars=init_par,fun = Two.Neg.Log.Likelihood.Use.Basis.Log,
                      eqfun = Equal.Constrain.Log, eqB = 0,
                      ineqfun = NULL, ineqLB = NULL, ineqUB = NULL,
                      control= NULL,Y = Y,
                      Log.All.Basis.Cor_matrix = Log.All.Basis.Cor_matrix,
                      approx.option = approx.option, approx.par = approx.par)
          
        }
        
        current.seed.list <- list(seed = i,
                                  m = m,
                                  est.par = sol$pars,
                                  est.fvalue = sol$values[length(sol$values)],
                                  est.convergence = sol$convergence,
                                  init_par = init_par)
      }
    
    
    
    stopImplicitCluster()
    stopCluster(cl)
  }
  
  if(!doParallel){
    
   all.seed.list = list()
   for(i in 1:n_seed){
        
        set.seed(i*100)
        
        #if(i == 1){
        #  init_par = LSE.w_tilde(Y, Log.All.Basis.Cor_matrix, pos.coef)
        
        #}else{
        init_par = runif(m,0,1)
        init_par = init_par/sum(init_par)
        #}
        #print(init_par)
        if(pos.coef){
          sol = solnp(pars=init_par,fun = Two.Neg.Log.Likelihood.Use.Basis.Log,
                      eqfun = Equal.Constrain.Log, eqB = 0,
                      ineqfun = NULL, ineqLB = NULL, ineqUB = NULL,
                      LB = rep(0,m), UB = rep(1,m), control= NULL,Y = Y,
                      Log.All.Basis.Cor_matrix = Log.All.Basis.Cor_matrix,
                      approx.option = approx.option, approx.par = approx.par)
        }else{
          sol = solnp(pars=init_par,fun = Two.Neg.Log.Likelihood.Use.Basis.Log,
                      eqfun = Equal.Constrain.Log, eqB = 0,
                      ineqfun = NULL, ineqLB = NULL, ineqUB = NULL,
                      control= NULL,Y = Y,
                      Log.All.Basis.Cor_matrix = Log.All.Basis.Cor_matrix,
                      approx.option = approx.option, approx.par = approx.par)
          
        }
        
        #est.par_matrix[i,] = sol$pars
        #est.fvalue_matrix[i,] = sol$values[length(sol$values)]
        #est.convergence_matrix[i,] = sol$convergence
        #print(sol$convergence)
        
        current.seed.list <- list(seed = i,
                                  m = m,
                                  est.par = sol$pars,
                                  est.fvalue = sol$values[length(sol$values)],
                                  est.convergence = sol$convergence,
                                  init_par = init_par)
   
        all.seed.list[[length(all.seed.list)+1]] = current.seed.list
   }
    
  }
  
  
  sum.convergence = sum(sapply(all.seed.list, function(x){x$est.convergence}))
  if(sum.convergence !=0){
    warnings("some seed do not converge !")
  }
  
  best.idx = which.min(sapply(all.seed.list, function(x){x$est.fvalue}))
  
  best.seed.list = all.seed.list[[best.idx]]
  best.seed.list$Log.All.Basis.Cor_matrix = Log.All.Basis.Cor_matrix
  
  if(all.seed.list.return){
    return(list(best.seed.list = best.seed.list,
                all.seed.list = all.seed.list))
  }else{
    return(list(best.seed.list = best.seed.list))
  }
}
