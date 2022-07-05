
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #############################################################
#' 
#' **Part IIII: check the correctness.**
#' 
#' 
#' *All functions:*
#' 
#' #### 1. Est.Sigma_and_Est.v and Est.Sigma_and_Est.v.Approx: Given est.par_matrix and Basis and Y, return est.Sigma_R and est.V_vector
#' 
#' #### 2. Plot.Results: Given est.par_matrix, est.V_matrix, use Basis, rho_matrix, True.Sigma_matrix,
#' ####                  plot true correlation function vs estimated correlation function
#' ####                  plot true covariance function vs estimated covariance function
#' 
#' 
#' 
#' 
#' 
#' 
#' #### 1. Est.Sigma_and_Est.v_Log and Est.Sigma_and_Est.v.Approx_Log: Given est.par_matrix and Log.Basis and Y, return est.Sigma_R and est.V_vector
#' 
#' #### 2. Plot.Results.Log 
#' ####                  : Given est.par_matrix, est.V_matrix, use Log.Basis, rho_matrix, True.Sigma_matrix,
#' ####                  plot true correlation function vs estimated correlation function
#' ####                  plot true covariance function vs estimated covariance function
#' 
#' 
#' 
#' 
#' 
#' #### Function 1
#' #### Description --- This function computes estimated Correlation matrix and estimated variance value.
#' ####                 Given a r by n matrix, calculate est.par_matrix, a n_seed by m matrix
#' 
#' #### Input -- est.par_matrix =  a n_seed by m matrix
#' ####       -- Y = a r by n matrix;
#' 
#' 
#' #### Return -- a list = est.Sigma_R_matrix: a n by n by n_seed matrix;
#' ####                  and est.V_vector: a vector of length n_seed
#' 
## ------------------------------------------------------------------------
Est.Sigma_and_Est.v <- function(est.par_matrix,All.Basis.Cor_matrix,Y,
                                approx.option,approx.par){
  
  n_seed = nrow(est.par_matrix)
  n = dim(All.Basis.Cor_matrix)[1]
  est.Sigma_R_matrix = array(0,dim=c(n,n,n_seed))
  est.V_vector = numeric(n_seed)
  
  for(i in 1:n_seed){
      est.Sigma_R_matrix[,,i] = Product.3d(w = est.par_matrix[i,],
                                           Basis = All.Basis.Cor_matrix)
      
   if(approx.option == 1){#not approx
     est.V_vector[i] = Cov.Zero(Y,est.Sigma_R_matrix[,,i])
   }
   if(approx.option == 2){#approx by using tol
      est.V_vector[i] = Cov.Zero.Approx(Y,est.Sigma_R_matrix[,,i],tol=approx.par) 
   }
   if(approx.option == 3){#approx by using tol.ratio
      est.V_vector[i] = Cov.Zero.Approx.ratio(Y,est.Sigma_R_matrix[,,i],
                                              tol.ratio=approx.par)
   }
}
  return(list(est.Sigma_R_matrix = est.Sigma_R_matrix,
              est.V_vector = est.V_vector))
}



#' 
#' 
#' 
#' 
## ------------------------------------------------------------------------
Est.Sigma_and_Est.v_Log <- function(est.par_matrix,Log.All.Basis.Cor_matrix,
                                    Y,approx.option,approx.par){
  
  n_seed = nrow(est.par_matrix)
  n = dim(Log.All.Basis.Cor_matrix)[1]
  est.Sigma_R_matrix = array(0,dim=c(n,n,n_seed))
  est.V_vector = numeric(n_seed)
  
  for(i in 1:n_seed){
    if(all(est.par_matrix[i,]>0)){
      est.Sigma_R_matrix[,,i] = Product.3d.Log(w = est.par_matrix[i,], 
                                             Log.Basis = Log.All.Basis.Cor_matrix)
    }else{
      est.Sigma_R_matrix[,,i] = Product.3d(w = est.par_matrix[i,], 
                                           Basis = exp(Log.All.Basis.Cor_matrix))
    }
    
    if(approx.option == 1){#not approx
      est.V_vector[i] = Cov.Zero(Y,est.Sigma_R_matrix[,,i])
    }
    if(approx.option == 2){#approx by using tol
      est.V_vector[i] = Cov.Zero.Approx(Y,est.Sigma_R_matrix[,,i],tol=approx.par) 
    }
    if(approx.option == 3){#approx by using tol.ratio
      est.V_vector[i] = Cov.Zero.Approx.ratio(Y,est.Sigma_R_matrix[,,i],
                                              tol.ratio=approx.par)
    }
  }
  return(list(est.Sigma_R_matrix = est.Sigma_R_matrix,
              est.V_vector = est.V_vector))
}




#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #### Function 2
#' #### Description -- This function shows some figures.
#' #### Use the same rho_matrix to generate different solutions est.par saved in est.par_matrix
#' 
#' #### Input -- est.par_matrix = a nrow by m matrix; each row is an estimation of m parameters
#' ####          est_V_vector = a nrow by m matrix;  use each row of Y.matrix to get each row of ####    est.par_matrix
#' ####          All.Basis.Cor_matrix = Basis matrix (a n by n by 5 array);
#' ####          calculated by function All.Basis.Cor()
#' ####   
#' $$Basis[i,j,k] =\frac{Beta(k+\rho_{[i,j]}^2,m-k+1)}{Beta(k,m-k+1)}$$
#' ####         True.Sigma_matrix = a n by n matrix; True Sigma Covariance matrix, calculated by Data_Generation.
#' ####         rho_matrix = a n by n distance matrix
#' 
#' 
#' ####         Cov_function = the function used for generating data, with option and par.cov as parameters.
#' 
#' #### Return -- Two figures
#' Figure 1: observations's true covariance and estimated covariance
#' 
#' Figure 2: true covariance function and estimated covariance function
#' 
#' 
## ------------------------------------------------------------------------

Plot.Results <- function(est.par_matrix,est.V_matrix,
                         All.Basis.Cor_matrix,
                         rho_matrix,True.Sigma_matrix,
                         Cov_function=Cov_function,option=option,par.cov=par.cov,
                         text.legend="est"){
  
#####
### True
#######
   
  true_Sigma_Cov = True.Sigma_matrix
  true_V = as.numeric(true_Sigma_Cov[1,1])
  true_Sigma_R = true_Sigma_Cov/true_V

  n = dim(All.Basis.Cor_matrix)[1]
  #m = dim(All.Basis.Cor_matrix)[3]
  

######
### Estimated
#######
  
  est_Sigma_Cov_matrix = est_Sigma_R_matrix = matrix(0,nrow = nrow(est.par_matrix),ncol=n*n)
 
  
  for(i in 1:nrow(est.par_matrix)){
  est_par = as.vector(est.par_matrix[i,])
  est_Sigma_R = Product.3d(w = est_par,Basis = All.Basis.Cor_matrix)
  
  est_Sigma_R_matrix[i,] = as.vector(est_Sigma_R)
  est_Sigma_Cov_matrix[i,] = est.V_matrix[i] * est_Sigma_R_matrix[i,]
  
}



#layout(matrix(1:2,1,2))
 

## Fig 1:covariance
#--------------------
  #true  
  
  rho_seq = seq(from=0,to=max(rho_matrix)+1,by=0.1)
  plot(as.vector(rho_matrix),as.vector(true_Sigma_Cov),type="p",
     ylim=c(0,max(max(est_Sigma_Cov_matrix),max(true_Sigma_Cov))+0.2),
     main="true Cov vs est Cov",
     xlab="dist",ylab="cov",cex=.2)
  lines(rho_seq,sapply(rho_seq,Cov_function,option,par.cov),lty="dashed")
  text(x=0+0.1,y=true_V,labels=round(true_V,4),pos=4)

  
  #est
  for(i in 1:nrow(est.par_matrix)){
    est_par = est.par_matrix[i,];est_V = est.V_matrix[i];
    #use rho_matrix
    points(as.vector(rho_matrix),est_V * est_Sigma_R_matrix[i,],col="red",cex=.2)
    #use rho_seq
    lines(rho_seq,
          sapply(rho_seq,Cor.Approx.2d, w_tilde = est_par)*est_V,
          col="red",lty="dashed")
    text(x=0+0.1,y=est_V,labels=round(est_V,4),pos=4)
    
  }
  
  #plot true again: for clarity
  lines(rho_seq,sapply(rho_seq,Cov_function,option,par.cov),lty="dashed")
  
  legend("topright",legend = c("true",paste(text.legend,1:nrow(est.par_matrix))),
         col=c("black",rep("red",nrow(est.par_matrix))),
         lty = c("dashed",rep("dashed",nrow(est.par_matrix))))

  
# Fig 2: Correlation
#----------------------------
  #true
  plot(as.vector(rho_matrix),as.vector(true_Sigma_R),type="p",
     ylim=c(0,max(max(est_Sigma_R),max(true_Sigma_R))+0.2),
     main="true Cor vs est Cor",col="black",
     xlab="dist",ylab="cor",cex=.2)
  lines(rho_seq,sapply(rho_seq,Cov_function,option,par.cov)/true_V,lty="dashed")

  #est
  for(i in 1:nrow(est.par_matrix)){
    est_par = est.par_matrix[i,];
    #use rho_matrix
    points(as.vector(rho_matrix),est_Sigma_R_matrix[i,],col="red",cex=.2)
    #use rho_seq
    lines(rho_seq,
          sapply(rho_seq,Cor.Approx.2d,w_tilde = est_par),
          col="red",lty="dashed")
  }
  
  #plot true again: clarity
  lines(rho_seq,sapply(rho_seq,Cov_function,option,par.cov)/true_V,lty="dashed")
  
   legend("topright",legend = c("true",paste(text.legend,1:nrow(est.par_matrix))),
         col=c("black",rep("red",nrow(est.par_matrix))),
         lty = c("dashed",rep("dashed",nrow(est.par_matrix))))
  
}




#' 
#' 
#' 
#' 
#' 
#' 
## ------------------------------------------------------------------------
Plot.Results.Log <- function(est.par_matrix,
                             est.V_matrix,
                         Log.All.Basis.Cor_matrix,
                         rho_matrix,True.Sigma_matrix,
                         Cov_function,option,par.cov,
                         text.legend="est",bound){
  
#####
### True
#######
   
  true_Sigma_Cov = True.Sigma_matrix
  true_V = as.numeric(true_Sigma_Cov[1,1])
  true_Sigma_R = true_Sigma_Cov/true_V

  n = dim(Log.All.Basis.Cor_matrix)[1]
  #m = dim(All.Basis.Cor_matrix)[3]
  

######
### Estimated
#######
  
  est_Sigma_Cov_matrix = est_Sigma_R_matrix = matrix(0,nrow = nrow(est.par_matrix),ncol=n*n)
 
  
  for(i in 1:nrow(est.par_matrix)){
  est_par = as.vector(est.par_matrix[i,])
  est_Sigma_R = Product.3d.Log(w = est_par,Log.Basis = Log.All.Basis.Cor_matrix)
  
  est_Sigma_R_matrix[i,] = as.vector(est_Sigma_R)
  est_Sigma_Cov_matrix[i,] = est.V_matrix[i] * est_Sigma_R_matrix[i,]
  
}



#layout(matrix(1:2,1,2))
 

## Fig 1:covariance
#--------------------
  #true  
  
  rho_seq = seq(from=0,to=bound,length.out = 1000)
  plot(as.vector(rho_matrix),as.vector(true_Sigma_Cov),type="p",
       ylim=c(0,max(max(est_Sigma_Cov_matrix),max(true_Sigma_Cov))+0.2),
     main="true Cov vs est Cov",xlim=c(0,bound),
     xlab="dist",ylab="cov",cex=.2)
  lines(rho_seq,sapply(rho_seq,Cov_function,option,par.cov),lty="dashed")
  text(x=0+0.1,y=true_V,labels=round(true_V,4),pos=4)

  
  #est
  for(i in 1:nrow(est.par_matrix)){
    est_par = est.par_matrix[i,];est_V = est.V_matrix[i];
    #use rho_matrix
    points(as.vector(rho_matrix),est_V * est_Sigma_R_matrix[i,],col="red",cex=.2)
    #use rho_seq
    lines(rho_seq,
          sapply(rho_seq,Cor.Approx.2d.Log, w_tilde = est_par)*est_V,
          col="red",lty="dashed")
    text(x=0+0.1,y=est_V,labels=round(est_V,4),pos=4)
    
  }
  
  #plot true again: for clarity, last layer is true
  lines(rho_seq,sapply(rho_seq,Cov_function,option,par.cov),lty="dashed", lwd=2)
  
  legend("topright",legend = c("true",paste(text.legend,1:nrow(est.par_matrix))),
         col=c("black",rep("red",nrow(est.par_matrix))),
         lty = c("dashed",rep("dashed",nrow(est.par_matrix))))

  
# Fig 2: Correlation
#----------------------------
  #true
  plot(as.vector(rho_matrix),as.vector(true_Sigma_R),type="p",
       ylim=c(0,max(max(est_Sigma_R),max(true_Sigma_R))+0.2),
     main="true Cor vs est Cor",col="black",xlim=c(0,bound), 
     xlab="dist",ylab="cor",cex=.2)
  lines(rho_seq,sapply(rho_seq,Cov_function,option,par.cov)/true_V,lty="dashed")

  #est
  for(i in 1:nrow(est.par_matrix)){
    est_par = est.par_matrix[i,];
    #use rho_matrix
    points(as.vector(rho_matrix),est_Sigma_R_matrix[i,],col="red",cex=.2)
    #use rho_seq
    lines(rho_seq,
          sapply(rho_seq,Cor.Approx.2d.Log,w_tilde = est_par),
          col="red",lty="dashed")
  }
  
  #plot true again: for clarity, last layer is true
  lines(rho_seq,sapply(rho_seq,Cov_function,option,par.cov)/true_V,lty="dashed", lwd=2)
  
   legend("topright",legend = c("true",paste(text.legend,1:nrow(est.par_matrix))),
         col=c("black",rep("red",nrow(est.par_matrix))),
         lty = c("dashed",rep("dashed",nrow(est.par_matrix))))
  
}






Plot.Results.Log.Transform <- function(est.par_matrix,
                                       est.V_matrix,
                                       Log.All.Basis.Cor_matrix,
                                       rho_matrix.orig,True.Sigma_matrix,
                                       Cov_function,option,par.cov,
                                       text.legend="est",bound,transform.ratio){
  
  #####
  ### True
  #######
  
  true_Sigma_Cov = True.Sigma_matrix
  true_V = as.numeric(true_Sigma_Cov[1,1])
  true_Sigma_R = true_Sigma_Cov/true_V
  
  n = dim(Log.All.Basis.Cor_matrix)[1] #basis for transformed rho
  #m = dim(All.Basis.Cor_matrix)[3]
  
  
  ######
  ### Estimated
  #######
  
  est_Sigma_Cov_matrix = est_Sigma_R_matrix = matrix(0,nrow = nrow(est.par_matrix),ncol=n*n)
  
  
  for(i in 1:nrow(est.par_matrix)){
    est_par = as.vector(est.par_matrix[i,])
    est_Sigma_R = Product.3d.Log(w = est_par,Log.Basis = Log.All.Basis.Cor_matrix)
    
    est_Sigma_R_matrix[i,] = as.vector(est_Sigma_R)
    est_Sigma_Cov_matrix[i,] = est.V_matrix[i] * est_Sigma_R_matrix[i,]
    
  }
  
  
  
  #layout(matrix(1:2,1,2))
  
  
  ## Fig 1:covariance
  #--------------------
  #true  
  
  #orig vs true covariance
  rho_seq.orig = seq(from=0,to=bound,length.out = 1000)
  rho_seq = rho_seq.orig / transform.ratio
  
  plot(as.vector(rho_matrix.orig),as.vector(true_Sigma_Cov),type="p",
       ylim=c(0,max(max(est_Sigma_Cov_matrix),max(true_Sigma_Cov))+0.2),
       main="true Cov vs est Cov",xlim=c(0,bound),
       xlab="dist",ylab="cov",cex=.2)
  
  lines(rho_seq.orig,sapply(rho_seq.orig,Cov_function,option,par.cov),lty="dashed")
  text(x=0+0.1,y=true_V,labels=round(true_V,4),pos=4)
  
  
  #est
  #orig vs true
  for(i in 1:nrow(est.par_matrix)){
    est_par = est.par_matrix[i,];est_V = est.V_matrix[i];
    #use rho_matrix
    points(as.vector(rho_matrix.orig),est_V * est_Sigma_R_matrix[i,],col="red",cex=.2)
    #use rho_seq
    lines(rho_seq.orig,
          sapply(rho_seq,Cor.Approx.2d.Log, w_tilde = est_par)*est_V,
          col="red",lty="dashed")
    text(x=0+0.1,y=est_V,labels=round(est_V,4),pos=4)
    
  }
  
  #plot true again: for clarity, last layer is true
  lines(rho_seq.orig,sapply(rho_seq.orig,Cov_function,option,par.cov),lty="dashed", lwd=2)
  
  legend("topright",legend = c("true",paste(text.legend,1:nrow(est.par_matrix))),
         col=c("black",rep("red",nrow(est.par_matrix))),
         lty = c("dashed",rep("dashed",nrow(est.par_matrix))))
  
  
  # Fig 2: Correlation
  #----------------------------
  #true
  plot(as.vector(rho_matrix.orig),as.vector(true_Sigma_R),type="p",
       ylim=c(0,max(max(est_Sigma_R),max(true_Sigma_R))+0.2),
       main="true Cor vs est Cor",col="black",xlim=c(0,bound), 
       xlab="dist",ylab="cor",cex=.2)
  lines(rho_seq.orig,sapply(rho_seq.orig,Cov_function,option,par.cov)/true_V,lty="dashed")
  
  #est
  for(i in 1:nrow(est.par_matrix)){
    est_par = est.par_matrix[i,];
    #use rho_matrix
    points(as.vector(rho_matrix.orig),est_Sigma_R_matrix[i,],col="red",cex=.2)
    #use rho_seq
    lines(rho_seq.orig,
          sapply(rho_seq,Cor.Approx.2d.Log,w_tilde = est_par),
          col="red",lty="dashed")
  }
  
  #plot true again: for clarity, last layer is true
  lines(rho_seq.orig,sapply(rho_seq.orig,Cov_function,option,par.cov)/true_V,lty="dashed", lwd=2)
  
  legend("topright",legend = c("true",paste(text.legend,1:nrow(est.par_matrix))),
         col=c("black",rep("red",nrow(est.par_matrix))),
         lty = c("dashed",rep("dashed",nrow(est.par_matrix))))
  
}




