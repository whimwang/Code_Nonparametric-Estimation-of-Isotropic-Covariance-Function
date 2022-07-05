
Plot.List.Result <- function(True.Sigma_matrix, rho_support, par.cov, option,
                             rho_matrix.orig,
                             all.best.seed.list,
                             epsilon,
                             text.legend = "best "){

  
  #### bound #####
  m.vector = sapply(all.best.seed.list, function(x){return(x$m)})
  ratio.vector = sapply(all.best.seed.list, function(x){return(x$transform.ratio)})
  bound.basis = max(sapply(m.vector, Support.given.m.m, epsilon) * ratio.vector)
  bound.true = max(rho_matrix.orig)
  bound = max(rho_support, bound.true)*1.5
  bound = ifelse(bound.basis/bound.true >2, bound, max(bound, bound.basis))
   
  
  #### True V and True Cor
  true_V = True.Sigma_matrix[1,1]
  true.Sigma_R = True.Sigma_matrix/true_V
  
  
  
  ## Fig 1:True Covariance Points
  #-------------------
  max.est.V = max(sapply(all.best.seed.list, function(x){x$est.V}))
  plot(NA,
       ylim=c(0,max(max.est.V,max(True.Sigma_matrix))+0.2),
       main="true Sample Cov",xlim=c(0,bound),
       xlab="dist",ylab="cov")
  
  points(as.vector(rho_matrix.orig),as.vector(True.Sigma_matrix),cex=.8,pch=16)
  
  rho_seq.orig = seq(from=0,to=bound,length.out = 1000)
  lines(rho_seq.orig,sapply(rho_seq.orig,Cov_function,option,par.cov),lty="dashed")
  text(x=0+0.1,y=true_V,labels=round(true_V,4),pos=4)
  
  
  
  
  ## Fig 2: True Cov vs Est Cov lines
  plot(NA,
       ylim=c(0,max(max.est.V,max(True.Sigma_matrix))+0.2),
       main="true Cov vs est Cov",xlim=c(0,bound),
       xlab="dist",ylab="cov")
  
  #est
  #orig vs true
  N= length(all.best.seed.list)
  for(i in 1:N){
    rho_seq = rho_seq.orig / all.best.seed.list[[i]]$transform.ratio
    est_R = all.best.seed.list[[i]]$est.R;
    est_V = all.best.seed.list[[i]]$est.V;
    #use rho_matrix
    points(as.vector(rho_matrix.orig),est_V * est_R,col="red",cex=.3,pch=16)
    #use rho_seq
    lines(rho_seq.orig,
          sapply(rho_seq,Cor.Approx.2d.Log, w_tilde = all.best.seed.list[[i]]$est.par)*est_V,
          col="red",lty="dashed")
    text(x=0+0.1,y=est_V,labels=round(est_V,4),pos=4)
    
  }
  
  #plot true again: for clarity, last layer is true
  lines(rho_seq.orig,sapply(rho_seq.orig,Cov_function,option,par.cov),lty="dashed", lwd=2)
  
  legend("topright",legend = c("true",paste(text.legend,1:N,rep(":",N),m.vector)),
         col=c("black",rep("red",N)),
         lty = c("dashed",rep("dashed",N)))
  
  
  # Fig 2: Correlation
  #----------------------------
  #true
  plot(NA,
       ylim=c(0,1.2),
       main="true Cor vs est Cor",xlim=c(0,bound),
       xlab="dist",ylab="cor")
  
  points(as.vector(rho_matrix.orig),as.vector(true.Sigma_R),cex=.5)
  lines(rho_seq.orig,sapply(rho_seq.orig,Cov_function,option,par.cov)/true_V,lty="dashed")
  
  #est
  for(i in 1:N){
    est_par = all.best.seed.list[[i]]$est.par;
    #use rho_matrix
    points(as.vector(rho_matrix.orig),as.vector(all.best.seed.list[[i]]$est.R),col="red",cex=.2)
    #use rho_seq
    lines(rho_seq.orig,
          sapply(rho_seq,Cor.Approx.2d.Log,w_tilde = est_par),
          col="red",lty="dashed")
  }
  
  #plot true again: for clarity, last layer is true
  lines(rho_seq.orig,sapply(rho_seq.orig,Cov_function,option,par.cov)/true_V,lty="dashed", lwd=2)
  
  legend("topright",legend = c("true",paste(text.legend,1:N,rep(":",N),m.vector)),
         col=c("black",rep("red",N)),
         lty = c("dashed",rep("dashed",N)))
  
}

