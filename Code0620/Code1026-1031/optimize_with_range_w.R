#' ---
#' title: "range_parameter.Rmd"
#' author: "Yiming Wang"
#' date: "10/25/2019"
#' output: html_document
#' ---
#' 
#' #' Function
#' #' This function calculate the MLE of correlation $\Sigma_R^m(D)=\sum_{k=1}^mA_m^k(\frac{\rho}{\theta_{range}})$
#' @parameter: range: a scalar; range parameter
#' @parameter: rho_matrix: distance matrix; (n,n)
#' @parameter: w_tilde: a vector of length m; sum equal to 1
#' @parameter: Y: observation matrix, a (r,n) matrix
#' 
#' @return: a scalar
#' #'Step 1: calculate basis 
#' $R^m(\frac{D}{\theta_{range}})=\sum_{k=1}^mA_m^k(\frac{\rho}{\theta_{range}})$
#' #'Step 2: Use Two.Neg.Log.Likelihood.Use.Basis.Log to calculate 
#' $nrlog(\sum_{i=1}^r y_i^T R^{-1}y_i)+rlogdet(R)$
#' 
#' 
#' 
## ------------------------------------------------------------------------
Two.Neg.Log.Likelihood.Range.func <- function(range, w_tilde,Y, rho_matrix,
                                              positive.check = TRUE, pd.check = TRUE,pd.adjust = TRUE,
                                              approx.option, approx.par){
  m = length(w_tilde)
  Log.All.Basis.Cor_matrix = Log.Basis.Cor.Final(rho_matrix/range, m, positive.check = positive.check, pd.check = pd.check,pd.adjust = pd.adjust)
  
  return(Two.Neg.Log.Likelihood.Use.Basis.Log(w_tilde,Y,Log.All.Basis.Cor_matrix,approx.option,approx.par))
}

#'
#' Function
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #' Function
#' #' This function find the best range parameter so that
#' $\theta_{range}=argmin_{\theta_{range}}[nrlog(\sum_{i=1}^r y_i^T R^{-1}y_i)+rlogdet(R)]$, where $R=\sum_{k=1}^mA_m^k(\frac{\rho}{\theta_{range}})$
#' 
#' @parameter: w_tilde.init, a vector of length m
#' @parameter: Y, observation matrix, (r,n)
#' @parameter: rho_matrix, a distance matrix, (n,n)
#' @parameter: upper.iter
#' @parameter: upper.convergence
#' @return: a list
#' 
#' #' Step 1: Initialize 
#' $\tilde_{w}=(\frac{1}{m},\cdots,\frac{1}{m})$
#' #' Iteration Starts:
#' #'              Given $\tilde{w}$, solve range with mimum objective value
#' #'              Given range, solve best $\tilde{w}$
#' #' Stop if iter.time > upper.iter or $|f_{new}-f_{old}|<|f_{new}|*upper.convergence$
#' 
## ------------------------------------------------------------------------
optimize_with_range_w <-function(w_tilde.init, Y, rho_matrix, 
                                 upper.iter = 10, upper.convergence = 10^(-3),
                                 positive.check = TRUE, pd.check = TRUE,pd.adjust = TRUE,
                                 approx.option, approx.par,pos.coef){
  #added myself
  #set.seed(0)
  
  time.1 = Sys.time()
  #====== Initialization =======
  m = length(w_tilde.init)
  iter.time = 0
  w_tilde.new = w_tilde.init
  range.new = range.current  = (min(rho_matrix[rho_matrix!=0]) + max(rho_matrix[rho_matrix!=0]))/2
  fvalue.new = 1000
  fvalue.current = 2 * 1000
  
  range.vector = NULL
  fvalue.vector = NULL
  w_tilde.matrix = NULL
  convergence.vector = NULL
  
  
  #w_tilde.init = c(1,rep(0,m-1))
  #w_tilde.init = c(rep(0,m-1),1)
  
  while(abs(fvalue.current - fvalue.new) > upper.convergence * abs(fvalue.current) && iter.time < upper.iter){
    iter.time = iter.time + 1
    print(paste0("This is the ",iter.time," iteration"))
    
    #===== Given w_tilde, optimize with respect to range parameter =========#
    
    w_tilde.current = w_tilde.new
    fvalue.current = fvalue.new
    
    
    #optimize: require unimodal: golden section search
    #temp.result = optimise(Two.Neg.Log.Likelihood.Range.func, lower = min(rho_matrix[rho_matrix!=0]), upper = max(rho_matrix[rho_matrix!=0]),
    #        maximum = FALSE, w_tilde = w_tilde.current, Y = Y, rho_matrix=rho_matrix, positive.check = positive.check,
    #        pd.check = pd.check, approx.option = approx.option, approx.par = approx.par)
    
    
    
    #grid search
    if(FALSE){
      range_seq = seq(from=min(rho_matrix[rho_matrix!=0]),to=max(rho_matrix[rho_matrix!=0]),length.out = 20)
      fvalue_seq = numeric(length(range_seq))
      
      start.time = Sys.time()
      for(i in 1:length(range_seq)){
        print(i)
        fvalue_seq[i] = Two.Neg.Log.Likelihood.Range.func(range_seq[i], w_tilde = w_tilde.current,Y,rho_matrix,
                                                          positive.check, pd.check, pd.adjust,
                                                          approx.option, approx.par)
        
      }
      end.time = Sys.time()
      print(end.time-start.time)
      idx = which.min(fvalue_seq)
      fvalue.new = fvalue_seq[idx]
      range.new = range_seq[idx]
      fvalue.vector = c(fvalue.vector, fvalue.new)
      range.vector = c(range.vector, range.new)
    }
    
    #GA
    if(FALSE){
      #GA can somehow get a value closer to global maximum but cannot tell which is better among GA and grid search but may take longer time.
      time.start = Sys.time()
      GA = ga(type = "real-valued",fitness = function(range){
        return(-Two.Neg.Log.Likelihood.Range.func(range, w_tilde = w_tilde.current,Y, rho_matrix,
                                                  positive.check = TRUE, pd.check = TRUE,pd.adjust = TRUE,
                                                  approx.option, approx.par))},
        lower = min(rho_matrix[rho_matrix!=0]),upper = max(rho_matrix[rho_matrix!=0]), 
        maxiter = 5, popSize = 20, pcrossover = 0.2, pmutation = 0.1)
      
      time.end = Sys.time()
      
      plot(range_seq, fvalue_seq)
      points(GA@solution,-GA@fitnessValue,col="red")
    }
    
    #random search
    if(TRUE){
      #random search # 50 points 5-6mins
      time.start = Sys.time()
      
    

      #Idea 1: generate from range of rho_matrix
      if(is.null(convergence.vector)){ #first time
        range_seq = c(runif(n=max(30,nrow(rho_matrix)/4), min=min(rho_matrix[rho_matrix!=0]),max=max(rho_matrix[rho_matrix!=0])), range.current) # generate uniform from the range 
      }else{
        range_seq = c(rnorm(n=30,mean=range.current,sd=max(abs(range_seq[order(fvalue_seq,decreasing = FALSE)]-range.current))/2),range.current)
        range_seq = range_seq[range_seq > 0]
      }
     
      
      #Idea 2: generate range in  neighborhood
      if(FALSE){
        if(iter.time == 1){
          range_seq = c(runif(n=max(50,nrow(rho_matrix)/4), min=min(rho_matrix[rho_matrix!=0]),max=max(rho_matrix[rho_matrix!=0])), range.current) # generate uniform from the range 
        }else{
          generate_range_orig = range(rho_matrix[rho_matrix!=0])
          range_sd = sd(rho_matrix[rho_matrix!=0])
          generate_range_new = c(range.current - 2*range_sd, range.current + 2*range_sd)
          range_seq = c(runif(n=max(30,nrow(rho_matrix)/4), min=max(generate_range_orig[1], generate_range_new[1]),
                              max=min(generate_range_orig[2], generate_range_new[2])), range.current) # generate uniform from the range 
        }
      }
      
      #fvalue_seq = numeric(length(range_seq))
      fvalue_seq = rep(Inf,length(range_seq))
      start.time = Sys.time()
    
      for(i in 1:length(range_seq)){
        tryCatch({
          fvalue_seq[i] = Two.Neg.Log.Likelihood.Range.func(range_seq[i], w_tilde = w_tilde.current,Y,rho_matrix,
                                                                    positive.check, pd.check, pd.adjust,
                                                                    approx.option, approx.par)
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
      w_tilde.matrix = rbind(w_tilde.matrix, w_tilde.current)#fixed w_tilde
      
      time.end = Sys.time()
    }
    
   
    
    #===== Given range parameter, optimize with respect to w_tilde =========
   
    range.current = range.new
    w_tilde.current = w_tilde.new
    fvalue.current = fvalue.new
    
    Log.All.Basis.Cor_matrix = Log.Basis.Cor.Final(rho_matrix/range.current, m, positive.check, pd.check, pd.adjust)
    
    if(pos.coef){
      #w_tilde.current = c(1,rep(0.0001,19))
      #w_tilde.current = w_tilde.current/sum(w_tilde.current)
      sol = solnp(pars=w_tilde.current,fun = Two.Neg.Log.Likelihood.Use.Basis.Log,
                  eqfun = Equal.Constrain.Log, eqB = 0,
                  ineqfun = NULL, ineqLB = NULL, ineqUB = NULL,
                  LB = rep(0,m), UB = rep(1,m), control= NULL,Y = Y,
                  Log.All.Basis.Cor_matrix = Log.All.Basis.Cor_matrix,
                  approx.option = approx.option, approx.par = approx.par)
    }else{
      sol = solnp(pars=w_tilde.current,fun = Two.Neg.Log.Likelihood.Use.Basis.Log,
                  eqfun = Equal.Constrain.Log, eqB = 0,
                  ineqfun = NULL, ineqLB = NULL, ineqUB = NULL,
                  control= NULL,Y = Y,
                  Log.All.Basis.Cor_matrix = Log.All.Basis.Cor_matrix,
                  approx.option = approx.option, approx.par = approx.par)
      
    }
    
    
    w_tilde.new = sol$pars
    fvalue.new = sol$values[length(sol$values)]
  
    convergence.vector = c(convergence.vector,sol$convergence)
    fvalue.vector = c(fvalue.vector, fvalue.new)
    range.vector = c(range.vector, range.current) #fixed range
    w_tilde.matrix = rbind(w_tilde.matrix, w_tilde.new)#update w_tilde
  
  }
  
  #plot(a1$range_seq,a1$fvalue_seq,ylim=range(c(a1$fvalue_seq,a2$fvalue_seq)))
  #points(a2$range_seq,a2$fvalue_seq,col='red')
  
  
  time.2 = Sys.time()
  
  
  return(list(est.w_tilde = w_tilde.new,
              est.range = range.new,
              est.fvalue = fvalue.new,
              fvalue.vector = fvalue.vector,
              range.vector = range.vector,
              convergence.vector = convergence.vector,
              w_tilde.matrix = w_tilde.matrix))
}


#' 
#' 
#' Function
#' This function 
#' Step 1. Generate w_tilde.init
#' Step 2. Given range, update w_tilde and Given w_tilde, update range

#'
#' Function
#' 
#' 
#' 
#' 
#' 
#' 

