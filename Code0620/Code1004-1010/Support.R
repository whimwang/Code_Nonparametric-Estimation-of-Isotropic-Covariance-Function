
#' 
#' Function: Support.Basis
#' Description: This function calculate the upper bound of numerical support of m basis
#' @param m: total number of basis function
#' @param epsilon: accuracy tol
#' 
#' The solution of $A^m_m(\rho)=\epsilon$
#' $$(1+\frac{\rho^2}{m})^{-1}=\epsilon$$
#' $$\frac{\rho^2}{m}=\epsilon^{-1}-1$$
#' 
#' The solution of $A^k_m(\rho)=\epsilon$
#' 
## ------------------------------------------------------------------------
Support.given.m.m <- function(m,epsilon=10^(-3)){#
  sqrt(m * (epsilon^(-1) - 1))
}

Support.given.k.m <- function(k,m,epsilon=10^(-3)){
  if(k == m){
    return(Support.given.m.m(m,epsilon))
  }
  target.func <- function(rho){
    Log.Basis.Cor(rho,k,m)-log(epsilon)
  }
  upper.bound = Support.given.m.m(m,epsilon)
  return(uniroot(target.func, interval = c(0,upper.bound))$root)
  
}

#' 
#' 
## ----eval=FALSE----------------------------------------------------------
## epsilon = 10^(-5)
## #support of A_m^m
## m_seq = c(seq(from=5,to=40, by=5),seq(from=50,80,by=10))
## support.seq = sapply(m_seq, Support.given.m.m,epsilon = epsilon)
## plot(support.seq, m_seq, xlim=c(0,max(support.seq)),xlab="m",main="numerical support")
## 
## 
## #support of A^1_m to A^m_m
## 
## for(m in m_seq){
##   k_seq = seq(from=1,to=m,by=1)
##   support_current_m = sapply(k_seq, Support.given.k.m, m = m, epsilon = epsilon)
##   result = list(m = m, support_m = support_current_m)
## 
##   assign(x = paste0("result_",m), result)
## }
## 
## 
## for(m in m_seq){
##   points(get(paste0("result_",m))$support_m, rep(m,m),col=m,cex=0.5)
## }
## 
## 
## #Plot density
## 
## labels = NULL
## total.support = NULL
## for(m in m_seq){
##   labels = c(labels, rep(m,m))
##   total.support = c(total.support,get(paste0("result_",m))$support_m)
## }
## labels = as.factor(labels)
## all.result = data.frame(total.support = total.support, labels = labels)
## 
## library("ggplot2")
## ggplot(all.result, aes(x= total.support, group=labels,color=labels),alpha=0.5) + geom_density(adjust=2)
## 
## 
## 
## library("ggplot2")
## ggplot(all.result, aes(x= total.support, group=labels,color=labels),alpha=0.5) + geom_density(adjust=2) + xlim(0,100)
## 

#' 
#' If we only consider to approximate using 50 basis or less to estimate covariance function, since most support are below 10( or below 5), so we have to make sure we can transform max distance to be 10 or 5.
#' 
#' 
#' Function:
#' Description: 
#' @param bound: 
#' @return: basis m
#' 
#' $$m = \frac{\rho^2}{\epsilon^{-1}-1}$$
#' 
## ------------------------------------------------------------------------

m.given.support <- function(bound, epsilon=10^(-3)){
    return(floor(exp(2 * log(bound) - log(epsilon^(-1)-1))))
  
}



#' 
#' 
#' Function: Plot Basis
#' Description:
#' 
#' 
## ------------------------------------------------------------------------
Plot.Basis <- function(m,epsilon=0.01){
  bound = Support.given.m.m(m, epsilon)
  
  plot(NA,xlab="distance",ylab="basis values",xlim=c(0,bound),ylim=c(0,1))
  for(k in 1:m){
    g=function(r){exp(Log.Basis.Cor(r,k,m))}
    curve(g,0,bound,add=T,col=k,lty=k)
  }
  legend("topright",legend=paste("k=",1:m),col=1:m,lty=1:m,bty="n",cex=0.5)
  title(paste("Basis functions with m=",m))
  
}


#' 
#' 
#' Test Plot.Basis function
## ----eval=FALSE----------------------------------------------------------
## layout(matrix(1:4,2,2))
## for(m in c(2,4,6,8)){
##   Plot.Basis(m)
## }
## 

#' 
#' 
#' 
#' 
#' Function:
#' Description: This function calculate the numerical support of Exp Covariance function
#' @param par.cov
#' @return rho
#' 
#' $$exp(-\frac{\rho}{\theta_2})=\epsilon$$
#' $$\rho = -\theta_2 * log(\epsilon)$$
#' 
## ------------------------------------------------------------------------
Support.Exp <- function(par.cov,epsilon=0.01){
  return(-par.cov[2] * log(epsilon))
}

#' 
#' 
#' 
#' 
#' Function
#' Description: This function computes the support of Gaussian
#' $$\theta_1 exp(-(\frac{\rho}{\theta_2})^2)=\epsilon$$
#' $$\rho = \theta_2\sqrt{log(\theta_1)-log(\epsilon)}$$
#' 
## ------------------------------------------------------------------------
Support.Gaussian <- function(par.cov,epsilon=0.01){
   return(par.cov[2] * sqrt(-log(epsilon)))
}

#' 
#' 
#' $$\theta_1\frac{\theta_2}{\rho^2+\theta_2}=\epsilon$$
#' $$\rho = \sqrt{\frac{\theta_1\theta_2}{\epsilon}-\theta_2}$$
## ------------------------------------------------------------------------
Support.RQ <- function(par.cov, epsilon=0.01){
  return(par.cov[2] * sqrt(1/epsilon-1))
  
}


# support of Matern
Support.Matern <- function(par.cov, epsilon){
  if(length(par.cov)!=3){
    stop("par.cov of Matern should be length 3!")
  }
  temp.func <- function(rho){
    return(Cov_Matern(rho, theta_1 = 1, theta_2 = par.cov[2], theta_3 = par.cov[3])-epsilon)
  }
  return(uniroot(temp.func, lower=0, upper=100000, extendInt = "yes")$root)  #extentint can extend the interval
}

# support of GenCauchy
Support.GenCauchy <- function(par.cov, epsilon){
  if(length(par.cov)!=4){
    stop("par.cov of Matern should be length 3!")
  }
  
 return((epsilon^(-par.cov[4]/par.cov[3])-1)^(1/par.cov[4]) * par.cov[2])
}


# support of Circular
Support.Circular <- function(par.cov, epsilon){
  temp.func <- function(rho){
    return(Cov_Circular(rho, theta_1 = 1, theta_2 = par.cov[2]) - epsilon)
  }
  return(uniroot(temp.func,lower = 0,upper = par.cov[2])$root)
  
}

# support of Spherical
Support.Spherical <- function(par.cov, epsilon){
  temp.func <- function(rho){
    return(Cov_Spherical(rho, theta_1 = 1, theta_2 = par.cov[2]) - epsilon)
  }
  return(uniroot(temp.func,lower = 0,upper = par.cov[2])$root)
}

# support of Cauchy
Support.Cauchy <- function(par.cov, epsilon){
  sqrt(epsilon^(-2)-1) * par.cov[2]
}

#support of linearmatern
Support.LinearMatern <- function(par.cov, epsilon){
  temp.func <- function(rho){
    return(Cov_LinearMatern(rho,theta_1=1,theta_2=par.cov[2],theta_3 = par.cov[3],
                            theta_4 = par.cov[4], theta_5 = par.cov[5])-epsilon)
  }
  uniroot(temp.func,lower = 0,upper = 100000000,extendInt = "yes")$root
}



#' 
#' 
## ------------------------------------------------------------------------
Support.function <- function(par.cov,option, epsilon){
  if(option == "Gaussian"){
    return(Support.Gaussian(par.cov,epsilon))
  }
  if(option == "Exp"){
    return(Support.Exp(par.cov, epsilon))
  }
  if(option == "RQ"){
    return(Support.RQ(par.cov, epsilon))
  }
  if(option == "Matern"){
    return(Support.Matern(par.cov, epsilon))
  }
  if(option == 'Cauchy'){
    return(Support.Cauchy(par.cov, epsilon))
  }
  if(option == 'LinearMatern'){
    return(Support.LinearMatern(par.cov, epsilon))
  }
  if(option == 'GenCauchy'){
    return(Support.GenCauchy(par.cov,epsilon))
  }
  if(option == "Circular"){
    return(Support.Circular(par.cov, epsilon))
  }
  if(option == "Spherical"){
    return(Support.Spherical(par.cov, epsilon))
  }
  stop('Invalid option!')
}


l.range.function <- function(rho_support,d,option.range){
  if(option.range != "max" && option.range != "mean"){
    stop("option range is not valid !")
  }
  # mean:
  if(option.range == "mean"){
   l = round(rho_support * sqrt(3/2/d),3)
  }
  
  if(option.range == "max"){
    l = round(rho_support/2/sqrt(d),3)
  }
  return(l)
}

#' 
#' 
## ------------------------------------------------------------------------
Plot.Basis.Option <- function(m,option, par.cov, epsilon=0.01){
  bound = Support.given.m.m(m, epsilon)
  
  #basis
  plot(NA,xlab="distance",ylab="basis values",xlim=c(0,bound),ylim=c(0,1))
  for(k in 1:m){
    g=function(r){exp(Log.Basis.Cor(r,k,m))}
    curve(g,0,bound,add=T,col=k+1,lty=k+1)
  }
  
  #true
  g = function(r){Cov_function(r, option = option, par.cov = par.cov)}
  curve(g,0,bound,add=T,col=1,lty=1,lwd=2)
  
  legend("topright",
         legend=c(paste(option,":", par.cov[1],par.cov[2]),paste("k=",1:m)),col=1:(m+1),lty=1:(m+1),bty="n",cex=0.5)
  title(paste(option,"par.cov:(",par.cov[1],",",par.cov[2],")","m=",m))
  
}



















##Distance of target Correlation function and Basis function

##$$\sqrt{\int_0^\infty |exp(-\frac{\rho}{\theta_2})-A_k^m(\rho)|^2d\rho}$$
  

Exp.Cor.Basis.Diff.Norm2 <- function(k,m,theta_2){
  
  target.func <- function(rho){
    (Cov_Exp(rho,theta_1 = 1, theta_2 = theta_2) - exp(Log.Basis.Cor(rho,k,m)))^2
  }
  
  return(sqrt(integrate(target.func, lower = 0, upper = Inf)$value))
  
}



Gaussian.Cor.Basis.Diff.Norm2 <- function(k,m,theta_2){
  
  target.func <- function(rho){
    (Cov_Gaussian(rho,theta_1 = 1, theta_2 = theta_2) - exp(Log.Basis.Cor(rho,k,m)))^2
  }
  
  return(sqrt(integrate(target.func, lower = 0, upper = Inf)$value))
  
}


RQ.Cor.Basis.Diff.Norm2 <- function(k,m,theta_2){
  
  target.func <- function(rho){
    (Cov_RQ(rho,theta_1 = 1, theta_2 = theta_2) - exp(Log.Basis.Cor(rho,k,m)))^2
  }
  
  return(sqrt(integrate(target.func, lower = 0, upper = Inf)$value))
  
}


Cor.Basis.Diff.Norm2 <- function(k,m,theta_2, option){
  if(option == "Gaussian")
    return(Gaussian.Cor.Basis.Diff.Norm2(k,m,theta_2))
  
  if(option == "Exp")
    return(Exp.Cor.Basis.Diff.Norm2(k,m,theta_2))
  
  if(option == "RQ")
    return(RQ.Cor.Basis.Diff.Norm2(k,m,theta_2))
  
}



