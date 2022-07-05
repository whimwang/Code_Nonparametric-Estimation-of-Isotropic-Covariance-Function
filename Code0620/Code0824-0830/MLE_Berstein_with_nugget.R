#' ---
#' title: "optimize_with_range_w_with_nugget"
#' author: "Whim Wang"
#' date: "9/24/2020"
#' output: html_document
#' ---
#' 
#' Calculate -2MLE/(r*n) given nugget,w_tilde,range
#' 
#' @param par, a vector of length m+1, first element is nugget effect and the rest is w_tilde
#' @param range, value, positive
#' @param Y, matrix of r by n
#' @param rho_matrix, dist matrix n by n
#' @param approx.option and approx.par
#' 
#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------
Two.Neg.Log.BersteinSpline.given_pars.with.nugget <- function(par,range,Y,rho_matrix,positive.check = TRUE,
                                                              pd.check = TRUE,pd.adjust = TRUE,
                                                              approx.option=2,approx.par = 10^(-8)){
  nugget = par[1];
  w_tilde = par[2:length(par)]
  Log.All.Basis.Cor_matrix = Log.Basis.Cor.Final(rho_matrix/range, m, positive.check = positive.check, 
                                                 pd.check = pd.check,pd.adjust = pd.adjust)
  Sigma = Product.3d(w_tilde,Basis = exp(Log.All.Basis.Cor_matrix)) + diag(1,nrow=dim(Log.All.Basis.Cor_matrix)[1])*nugget
  Two.Neg.Log.Likelihood.Use.Sigma(Y,Sigma,approx.option = approx.option,approx.par = approx.par)
}

#' 
#' 
#' Calculate -2MLE/(r*n) given Log.All.Basis.Cor_matrix
#' 
#' @param par, vector of length m+1
#' @param Y, matrix r by n
#' @param All.Basis.Cor_matrix, matrix of n by n by m
#' @param approx.option, approx.par
#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------
Two.Neg.Log.BersteinSpline.given_pars.Basis.with.nugget <- function(par,Y,All.Basis.Cor_matrix,approx.option=2,approx.par = 10^(-8)){
  nugget = par[1];
  w_tilde = par[2:length(par)]
  Sigma = Product.3d(w_tilde,Basis = All.Basis.Cor_matrix) + diag(1,nrow=dim(All.Basis.Cor_matrix)[1])*nugget
  Two.Neg.Log.Likelihood.Use.Sigma(Y,Sigma,approx.option = approx.option,approx.par = approx.par)
}


#' 
#' 
#' Estimate par (nugget and w_tilde) using MLE given par.init
#' @param par.init, vector of length m+1
#' @param Y, matrix r by n
#' @param rho_matrix, dist, n by n
#' @return list:est.par, est.range,est.fvalue,fvalue.vector,etc
#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------
optimize_with_range_w_nugget <-function(par.init, Y, rho_matrix, 
                                 upper.iter = 10, upper.convergence = 10^(-3),
                                 positive.check = TRUE, pd.check = TRUE,pd.adjust = TRUE,
                                 approx.option, approx.par,pos.coef){
  #added myself
  #set.seed(0)
  
  time.1 = Sys.time()
  #====== Initialization =======
  m = length(par.init)-1
  iter.time = 0
  par.new = par.init
  range.new = range.current  = (min(rho_matrix[rho_matrix!=0]) + max(rho_matrix[rho_matrix!=0]))/2
  fvalue.new = 1000
  fvalue.current = 2 * 1000
  
  range.vector = NULL
  fvalue.vector = NULL
  par.matrix = NULL
  convergence.vector = NULL
  
  
 
  
  while(abs(fvalue.current - fvalue.new) > upper.convergence * abs(fvalue.current) && iter.time < upper.iter){
    iter.time = iter.time + 1
    print(paste0("This is the ",iter.time," iteration"))
    
    #===== Given w_tilde, optimize with respect to range parameter =========#
    
    par.current = par.new
    fvalue.current = fvalue.new
    
    #random search
    if(TRUE){
      #random search # 50 points 5-6mins
      time.start = Sys.time()
      
      
      
      #Idea 1: generate from range of rho_matrix
      #if(is.null(convergence.vector)){ #first time
        #range_seq = c(runif(n=max(30,floor(nrow(rho_matrix)/4)), min=min(rho_matrix[rho_matrix!=0]),max=max(rho_matrix[rho_matrix!=0])), range.current) # generate uniform from the range 
        range_seq = c(runif(n=max(30), min=min(rho_matrix[rho_matrix!=0]),max=max(rho_matrix[rho_matrix!=0])), range.current) # generate uniform from the range 
      #}else{
      #  range_seq = c(rnorm(n=30,mean=range.current,sd=max(abs(range_seq[order(fvalue_seq,decreasing = FALSE)]-range.current))/2),range.current)
      #  range_seq = range_seq[range_seq > 0]
      #}
      
      
      #Idea 2: generate range in  neighborhood
     
      
      #fvalue_seq = numeric(length(range_seq))
      fvalue_seq = rep(Inf,length(range_seq))
      start.time = Sys.time()
      
      for(i in 1:length(range_seq)){
        tryCatch({
          fvalue_seq[i] = Two.Neg.Log.BersteinSpline.given_pars.with.nugget(par=par.current, range = range_seq[i],Y=Y,rho_matrix=rho_matrix,
                                                                            positive.check = positive.check,
                                                                            pd.adjust = pd.adjust,
                                                                            pd.check = pd.check,
                                                                            approx.option = approx.option,approx.par = approx.par)
        },
        error=function(e){
          print(paste0(i,' :range: ',range_seq[i]," error!"))
        }
        )
      }
      
      end.time = Sys.time()
      print(end.time-start.time)
      idx = which.min(fvalue_seq)
      fvalue.new = fvalue_seq[idx]
      range.new = range_seq[idx]
      
      #save new results
      fvalue.vector = c(fvalue.vector, fvalue.new)
      range.vector = c(range.vector, range.new) #update range
      par.matrix = rbind(par.matrix, par.current)#fixed w_tilde
      
      time.end = Sys.time()
    }
    
    
    
    #===== Given range parameter, optimize with respect to w_tilde =========
    
    range.current = range.new
    par.current = par.new
    fvalue.current = fvalue.new
    
    Log.All.Basis.Cor_matrix = Log.Basis.Cor.Final(rho_matrix/range.current, m, positive.check=positive.check, 
                                                   pd.check = pd.check, pd.adjust = pd.adjust)
    
    
    sol = solnp(pars=par.current,fun = Two.Neg.Log.BersteinSpline.given_pars.Basis.with.nugget,
                  eqfun = NULL, eqB = 0,
                  ineqfun = NULL, ineqLB = NULL, ineqUB = NULL,
                  LB = rep(0,m+1), UB = NULL, control= NULL,Y = Y,
                  All.Basis.Cor_matrix = exp(Log.All.Basis.Cor_matrix),
                  approx.option = approx.option, approx.par = approx.par)
    
    par.new = sol$pars
    fvalue.new = sol$values[length(sol$values)]
    
    convergence.vector = c(convergence.vector,sol$convergence)
    fvalue.vector = c(fvalue.vector, fvalue.new)
    range.vector = c(range.vector, range.current) #fixed range
    par.matrix = rbind(par.matrix, par.new)#update w_tilde
    
  }
  
  #plot(a1$range_seq,a1$fvalue_seq,ylim=range(c(a1$fvalue_seq,a2$fvalue_seq)))
  #points(a2$range_seq,a2$fvalue_seq,col='red')
  
  
  time.2 = Sys.time()
  
  
  return(list(est.par = par.new,
              est.range = range.new,
              est.fvalue = fvalue.new,
              fvalue.vector = fvalue.vector,
              range.vector = range.vector,
              convergence.vector = convergence.vector,
              par.matrix = par.matrix))
}

#' 
#' 
#' 
#' 
#' 
#' 
## ----eval=FALSE----------------------------------------------------------------------------------------------------------------------------------
## 
## m_seq = c(2,3,5,7,11,17)
## all.m_seq.result = list()
## for(i in 1:length(m_seq)){
##   m = m_seq[i]
##   par.init = c(0.5,rep(1/m,m))
## 
##   opt = optimize_with_range_w_nugget(par.init, Y, rho_matrix,
##                                  upper.iter = 10, upper.convergence = 10^(-3),
##                                  positive.check = TRUE, pd.check = TRUE,pd.adjust = TRUE,
##                                  approx.option, approx.par,pos.coef)
## 
##   all.m_seq.result[[i]] = opt
## }
## 
## 
## 
## all.est_fvalue = sapply(all.m_seq.result,function(x){x$est.fvalue})
## 
## all.est_nugget = sapply(all.m_seq.result,function(x){x$est.par[1]})
## 
## all.est_range = sapply(all.m_seq.result,function(x){x$est.range})
## 
## 
## fvalue_ratio = all.est_fvalue[2:length(all.est_fvalue)]/all.est_fvalue[1:(length(all.est_fvalue)-1)] - 1
## names(fvalue_ratio) = paste0(m_seq[2:length(m_seq)],":",m_seq[1:(length(m_seq)-1)])
## 
## 
## #relative ratio
## index2 = if(all.est_fvalue[length(m_seq)] > 0){
##   min(which(fvalue_ratio <0 & fvalue_ratio >- 1*10^(-3)))+1
## }else{
##   min(which(fvalue_ratio>0 & fvalue_ratio  < 1*10^(-3)))+1
## }
## 
## #which min
## index0 =  which.min(all.est_fvalue)
## 
## #which min not the last
## index1 = ifelse(index0 == length(m_seq),NA,index0)
## 
## #if min is the last one: find the relative ratio one and if no relative ratio then the last one
## #if min is the middle one; accept the middle one
## index_combined = ifelse(!is.na(index1),index1,index2)
## index_combined = ifelse(is.na(index_combined),index0,index_combined)
## 
## 
## final_result = all.m_seq.result[[index_combined]]
## 
## final_result$est.par[1]
## sum(final_result$est.par[-1])
## final_result$est.fvalue
## 
## 
## 
## #plot true vs est matrix
## for(i in 1:length(m_seq)){
##  opt = all.m_seq.result[[i]]
##  transformed.result = Transform.Distance(rho_matrix.orig, ratio = opt$est.range )# rho_matrix/ratio
##  rho_matrix = transformed.result$transformed_rho_matrix
##  transform.ratio = transformed.result$ratio #max of rho_matrix.orig
## 
## 
##  m = length(opt$est.par[-1])
##  Log.All.Basis.Cor_matrix = Log.Basis.Cor.Final(rho_matrix, m=m, positive.check = TRUE, pd.check = TRUE,pd.adjust = TRUE)
## 
##  est_Cov_matrix = Product.3d(w=opt$est.par[-1],Basis=exp(Log.All.Basis.Cor_matrix))
## 
##  if(i==1){
##    plot(Data.Preparation.result$Hidden.list$True.Sigma_matrix,
##       est_Cov_matrix)
##      abline(a=0,b=1,col='red')
##   }else{
##     points(Data.Preparation.result$Hidden.list$True.Sigma_matrix,est_Cov_matrix,pch=20)
##   }
## }
## 
## 
## 
## #plot d vs true, d vs est
## d_seq = seq(from=0,to=24,length.out = 100)
## f_cov_true <- function(x){
##         Cov_function(x,option = option,par.cov=par.cov)
## }
## 
##   ## to calculate metric
## plot(d_seq,sapply(d_seq, f_cov_true),type='l',lwd=5)
## 
## all_color  = 1:length(m_seq)+1
## for(i in 1:length(m_seq)){
##   opt = all.m_seq.result[[i]]
##   m = length(opt$est.par[-1])
## 
##   f_cov_est <- function(x){
##         tmp <- function(x){
##           vector = sapply(1:m, function(k){
##             return(-sum(log1p(rep(x^2/opt$est.range^2,m-k+1)/seq(from=k,to=m,by=1))))})
##           sum(exp(vector + log(opt$est.par[-1])))
##         }
##         sapply(x,tmp)
##   }
## 
##   lines(d_seq,sapply(d_seq, f_cov_est),col=all_color[i],lwd=2,lty='dashed')
## 
## }
## 
## legend('topright',legend = c('true',paste0('m=',m_seq)),lty=c('solid',rep('dashed',length(m_seq))),col=c('black',all_color))
## 

#' 
#' 
#' 
## ----eval=FALSE----------------------------------------------------------------------------------------------------------------------------------
## nugget = 1
## True_Sigma = Data.Preparation.result$Hidden.list$True.Sigma_matrix +
##   nugget * diag(x=1,nrow=nrow(Data.Preparation.result$Hidden.list$True.Sigma_matrix),ncol=ncol(Data.Preparation.result$Hidden.list$True.Sigma_matrix))
## Two.Neg.Log.Likelihood.Use.Sigma(Y,Sigma=True_Sigma,approx.option = approx.option,approx.par = approx.par)
## 

#' 
#' 
## ----eval=FALSE----------------------------------------------------------------------------------------------------------------------------------
## knitr::purl('MLE_Berstein_with_nugget.Rmd',documentation = 2L)

