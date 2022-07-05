#' ---
#' title: "Data_Generation_new.log.Rmd"
#' author: "Whim Wang"
#' date: "9/13/2019"
#' output: html_document
#' ---
#' 
#' 
#' 
#' 
#' 
#' #### Part III
#' $$C(h)=\int_0^\infty exp(-u^2h^2)P(u)du =\int_0^1 s^{h^2}\frac{1}{2s\sqrt{-log s}}P(\sqrt{-log s})ds=\int_0^1 s^{h^2}g(s)ds$$.
#' $g(s)$ is a function defined in $(0,1]$. And we use Berstein sum $g_m(s)$ to approximate $g(s)$.
#' 
#' $$g_m(s) = \sum_{k=1}^m w_k {m-1\choose k-1}s^{k-1}(1-s)^{m-k}ds$$
#' Covariance function
#' 
#' $$C(\rho)=\int_0^1 s^{\rho^2}g_m(s)ds
#'      =\sum_{k=1}^m w_k {m-1\choose k-1}\int_0^1 s^{\rho^2+k-1}(1-s)^{m-k}ds
#'       =\sum_{k=1}^m w_k {m-1\choose k-1}Beta(\rho^2+k,m-k+1)
#'       =\sum_{k=1}^m w_k \frac{1}{m}\frac{Beta(k+\rho^2,m-k+1)}{Beta(k,m-k+1)}$$
#'       
#' 
#' 
#' $$C(0)=\sum_{k=1}^m w_k\frac{1}{m}$$
#' 
#' Correlation function
#' 
#' $$R(\rho) = \frac{\sum_{k=1}^m w_k \frac{1}{m}\frac{Beta(k+\rho^2,m-k+1)}{Beta(k,m-k+1)}}{\sum_{k=1}^m w_k\frac{1}{m}}=\sum_{k=1}^m \frac{w_k}{\sum_{k=1}^m w_k}\frac{Beta(k+\rho^2,m-k+1)}{Beta(k,m-k+1)}=\sum_{k=1}^m \tilde{w}_k\frac{Beta(k+\rho^2,m-k+1)}{Beta(k,m-k+1)}$$.
#' 
#' And we have 
#' $$\sum_{k=1}^m \tilde{w}_k=1$$
#' 
#' Therefore we have two ways to approximate covariance.
#' $$C(\rho)=\sum_{k=1}^m w_k \frac{1}{m}\frac{Beta(k+\rho^2,m-k+1)}{Beta(k,m-k+1)}$$
#' 
#' $$C(\rho)=v_m\sum_{k=1}^m \tilde{w}_k \frac{Beta(k+\rho^2,m-k+1)}{Beta(k,m-k+1)}$$.
#' 
#' $$C(0)=v_m=\sum_{k=1}^m w_k\frac{1}{m}$$
#' Therefore, we have one way to approximate correlation.
#' $$R(\rho) = \sum_{k=1}^m \tilde{w}_k\frac{Beta(k+\rho^2,m-k+1)}{Beta(k,m-k+1)}$$.
#' 
#' 
#' 
#' We use l to denote log likelihood. 
#' When there is no replicates, $Y_{1 \times n}\sim N(0_n,\Sigma)$
#' $$l = -\frac{1}{2}log|\Sigma^{(m)}|-\frac{1}{2}y^T[\Sigma^{(m)}]^{-1}y$$
#' Therefore, there are two ways to calculate objective function.
#' $$\Sigma^{(m)}=v_m\Sigma_R^{(m)}$$
#' $$-2l = log|\Sigma^{(m)}|+y^T[\Sigma^{(m)}]^{-1}y$$
#' $$-2l = nlog(v_m) + log|\Sigma_R^{(m)}|+\frac{1}{v_m}y^T[\Sigma_R^{(m)}]^{-1}y$$
#' 
#' 
#' When there is r replicates, $Y_1,Y_2\cdots,Y_r iid \sim N(0_n,\Sigma)$
#' $$l = -\frac{r}{2}log|\Sigma^{(m)}|-\frac{1}{2}\sum_{i=1}^ry_i^T[\Sigma^{(m)}]^{-1}y_i$$
#' Therefore, there are two ways to calculate objective function.
#' $$\Sigma^{(m)}=v_m\Sigma_R^{(m)}$$
#' $$-2l = rlog|\Sigma^{(m)}|+\sum_{i=1}^r y_i^T[\Sigma^{(m)}]^{-1}y_i$$
#' $$-2l = rnlog(v_m) + rlog|\Sigma_R^{(m)}|+\frac{1}{v_m}\sum_{i=1}^r y_i^T[\Sigma_R^{(m)}]^{-1}y_i$$
#' 
#' 
#' To minimize the objective function,
#' $$min -2l(v_m,\tilde{w}_{1:m}) = rnlog(v_m) + rlog|\Sigma_R^{(m)}|+\frac{1}{v_m}\sum_{i=1}^r y_i^T[\Sigma_R^{(m)}]^{-1}y_i$$
#' 
#' $$\hat{v}_m=\frac{\sum_{i=1}^ry_i^T[\Sigma_R^{(m)}]^{-1}y_i}{rn}$$
#' Therefore, we substitute $\hat{v}_m$ to get
#' $$-2l(\tilde{w}_{1:m},\hat{v}_m)=rnlog(\sum_{i=1}^r y_i^T[\Sigma_R^{(m)}]^{-1}y_i)+rlog|\Sigma_R^{(m)}|-rnlog(rn)+ rn$$
#' 
#' Thus we only need to minimize
#' $$-2l(\tilde{w}_{1:m})=rnlog(\sum_{i=1}^r y_i^T[\Sigma_R^{(m)}]^{-1}y_i)+rlog|\Sigma_R^{(m)}|-rnlog(rn)+ rn$$
#' 
#' 
#' 
#' 
#' #### All Functions.
#' #### 1. Basis.Cor
#' #### 2. Cor.Approx.2d 
#' #### 3. Positive.Check
#' #### 4. PD.Check.2d and PD.Check.3d
#' ####    PD.Adjust.2d and PD.Adjust.3d
#' #### 5. All.Basis.Cor: 
#' #### 6. Cov.Zero and Cov.Zero.Approx and Cov.Zero.Approx.ratio: return sum_{1 to r} [y^T R^[-1] y]/(rn)
#' #### 7. Product.3d: return a n by n matrix; Basis[i,j,] %*% w
#' #### 8. LogDet.Inv
#' 
#' #### 9. Two.Neg.Log.Likelihood.Use.Basis: with 3 approx.option 
#' #### Use Basis to get R then calculate det(R) and v.hat
#' 
#' #### 10. Two.Neg.Log.Likelihood.Use.Basis.minus.one.par:  use a vector of  m-1 with 3 approx.option
#' 
#' #### 11. Two.Neg.Log.Likelihood.Use.Sigma_R and Two.Neg.Log.Likelihood.Use.Sigma
#' 
#' #### 12. Equal.Constrain: sum of a vector of length m minus 1 equal to zero
#' 
#' 
#' #### All Function for Log version
#' #### 1. Log.Basis.Cor : return log(Basis)
#' #### 2. Cor.Approx.2d.Log: return sum_k [exp(log(Basis_k) + log(w_k))]
#' #### 3. same
#' #### 4: same
#' #### 5: Log.All.Basis.Cor.Log: return log(Basis[,,]), a n by n by m array
#' #### 6: same
#' #### 7: Product.3d.Log: return a n by n matrix; sum(exp(log(Basis[i,j,]) + log(w)))
#' #### 8: same
#' 
#' #### 9:Two.Neg.Log.Likelihood.Use.Basis.Log: with 3 approx option 
#' #### Use log(Basis) to get R then calculate det(R) and v.hat
#' 
#' #### 10. Two.Neg.Log.Likelihood.Use.Basis.Log.minus.one.par:  use a length of m-1 with 3 approx.option
#' 
#' #### 11. same
#' 
#' #### 12. Equal.Constrain.Log: sum of a vector of length m minus 1 equal to zero.
#' 
#' 
#' 
#' #### Function 1
#' #### Description -- This function calculates the basis value of Correlation
#' #### Input - rho = distance   
#' ####         k = an integer in 1,2 to m
#' ####         m = an integer greater or equal to k and 1      
#' 
#' 
#' #### Return = Basis value
#' $$\frac{Beta(k+\rho^2,m-k+1)}{Beta(k,m-k+1)}$$
#' 
#' 
## ------------------------------------------------------------------------
Basis.Cor <- function(rho,k,m){
  return(beta(k+rho^2,m-k+1)/beta(k,m-k+1))
}

#' 
#' 
#' 
#' ##### Function 1
#' #### Description -- This function calculates the log basis value of correlation
#' #### Input - rho=distance
#' ####         method="Lbeta" or "Lgamma" 
#' ####         (LBeta and LGamma in Rfast, more effection than lbeta and lgamma in MASS)
#' 
#' Log Beta
#' $$LBeta(k+\rho^2,m-k+1)-LBeta(k,m-k+1)$$
#' 
#' Log Gamma
#' $$log(\frac{\frac{\Gamma(k+\rho^2)\Gamma(m-k+1)}{\Gamma(m+\rho^2+1)}}{\frac{\Gamma(k)\Gamma(m-k+1)}{\Gamma(m+1)}})=[LGamma(k+\rho^2)-LGamma(m+\rho^2+1)]-[LGamma(k)-LGamma(m+1)]$$
#' 
#' 
## ------------------------------------------------------------------------
#library("Rfast")
#Log.Basis.Cor <- function(rho,k,m,method){
#  if(method == "Lgamma"){
#    return(Lgamma(k+rho^2) - Lgamma(m+rho^2+1) - Lgamma(k) + Lgamma(m+1))
#  }
#  if(method == "Lbeta"){
#    return(Lbeta(k+rho^2,m-k+1) - Lbeta(k,m-k+1))
#  }
#}
#' Function: Log.Basis.Cor
#' Description:  This function calculates the log basis value of correlation 
#' 
#' @param rho: distance
#' @param k: current basis
#' @param m: total basis
#' 
#' Note: log1p(x) can calculate log(1+x) accurately when |x|<1
#' 
#' $$\frac{\frac{\Gamma(k+\rho^2)\Gamma(m-k+1)}{\Gamma(m+\rho^2+1)}}{\frac{\Gamma(k)\Gamma(m-k+1)}{\Gamma(m+1)}}=\prod_{j=k}^m(1+\frac{\rho^2}{j})^{-1}$$
#' $$log\frac{\frac{\Gamma(k+\rho^2)\Gamma(m-k+1)}{\Gamma(m+\rho^2+1)}}{\frac{\Gamma(k)\Gamma(m-k+1)}{\Gamma(m+1)}}=-\sum_{j=k}^mlog(1+\frac{\rho^2}{j})$$
#' 
#' 
## ------------------------------------------------------------------------
Log.Basis.Cor <- function(rho,k,m){
  
  vector = sapply(rho, function(rho.scalar){
    return(-sum(log1p(rho.scalar^2/seq(from=k,to=m,by=1))))}) # sapply so that we can use curve function
 # vector = -sum(log1p(rho^2/seq(from=k,to=m,by=1)))
  
  return(vector)
}


#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' Test Basis.Cor and Log.Basis.Cor
## ----eval=FALSE----------------------------------------------------------
## rho=sqrt(0)
## lgamma.ratio = exp(Log.Basis.Cor(rho,k=2,m=5,method = "Lgamma"))/Basis.Cor(rho,k=2,m=5)-1
## lbeta.ratio = exp(Log.Basis.Cor(rho,k=2,m=5,method = "Lbeta"))/Basis.Cor(rho,k=2,m=5)-1
## lgamma.ratio
## lbeta.ratio
## lgamma.ratio - lbeta.ratio
## 

#' 
#' 
#' 
#' 
#' 
#' #### Function 2
#' #### Description -- This function calculates the estimated Correlation.
#' #### Input - rho = distance;  w_tilde =  weights;
#' ####         Basis.Cor = function which is used for Basis.Cor
#' $$R(\rho) =\sum_{k=1}^m \tilde{w}_k\frac{Beta(k+\rho^2,m-k+1)}{Beta(k,m-k+1)}$$.
#' 
#' #### Return --  an estimated correlation value.
#' 
## ------------------------------------------------------------------------
## when Sigma
Cor.Approx.2d <- function(rho,w_tilde){
  m = length(w_tilde)
  basis = numeric(m)
  
  basis_func <- function(k){ 
    return(Basis.Cor(rho=rho,k,m))
  }
  
  basis = sapply(1:m,basis_func) 
  return(as.numeric(basis %*% w_tilde))
}


#' 
#' 
#' 
#' $$R(\rho) =\sum_{k=1}^m exp(log(\tilde{w}_k)+log(\frac{Beta(k+\rho^2,m-k+1)}{Beta(k,m-k+1)}))$$.
Cor.Basis.Log.rho_Rcpp_code = "

  Rcpp::NumericVector xcpp(rho_vector); //length n
  Rcpp::NumericVector w_cpp(w_tilde);   // length m
  Rcpp::NumericMatrix result_cpp(result); // matrix n by m
  
  int n = xcpp.length();
  int m_cpp = w_cpp.length();

  for (int i = 0; i < n; i++){
    for (int k = m_cpp; k >= 1; k--){
      if(k == m_cpp){
        result_cpp(i,k-1) = -log1p(double(pow(xcpp(i),2)/ k));
      }else{
        result_cpp(i,k-1) = -log1p(double(pow(xcpp(i),2)/ k)) + result_cpp(i,k);
      }
    }
  }
  
  double a=0;
  for(int k = 0; k < m_cpp; k++){
    a = log(w_cpp(k));
    for(int i = 0; i < n; i++){
      result_cpp(i,k) = exp(result_cpp(i,k) + a);
    }
  }

  //return sum(result_cpp(1,) * w_cpp)
  //return xcpp;
  "

Cor.Basis.Log.rho_Rcpp_func <- cxxfunction(signature(rho_vector = "numeric",w_tilde = "numeric",result = "numeric"), 
                                           body=Cor.Basis.Log.rho_Rcpp_code, 
                                           plugin="Rcpp")

## ------------------------------------------------------------------------
Cor.Approx.2d.Log <- function(rho,w_tilde){
  m = length(w_tilde)
 
  #Log.basis_func <- function(k){ 
  #  return(Log.Basis.Cor(rho=rho,k,m))
  #}
  tmp = matrix(0,nrow=length(rho),ncol=length(w_tilde))
  Cor.Basis.Log.rho_Rcpp_func(rho_vector = rho,w_tilde = w_tilde,tmp)
  
  return(apply(tmp,1,sum))
  #Log.basis = sapply(seq(from=1,to=m,by=1),Log.basis_func) 
  # Log.basis = sapply(1:5,Log.basis_func) can not be used here
  
  
  #return(sum(exp(Log.basis + log(w_tilde))))
}






#' 
#' 
#' test Cor.Approx.2d and Cor.Approx.2d.Log
## ----eval=FALSE----------------------------------------------------------
## rho = sqrt(200);m=5
## w_tilde <- runif(m,0,1)
## w_tilde <- w_tilde/sum(w_tilde)
## 
## Cor.Approx.2d(rho,w_tilde)
## Cor.Approx.2d.Log(rho,w_tilde,method="Lgamma")
## Cor.Approx.2d.Log(rho,w_tilde,method="Lbeta")
## 
## Cor.Approx.2d(rho,w_tilde)/Cor.Approx.2d.Log(rho,w_tilde,method="Lgamma")-1
## Cor.Approx.2d(rho,w_tilde)/Cor.Approx.2d.Log(rho,w_tilde,method="Lbeta")-1
## 
## Cor.Approx.2d.Log(rho,w_tilde,method="Lbeta")/Cor.Approx.2d.Log(rho,w_tilde,method="Lgamma")-1

#' 
#' 
#' 
#' 
#' 
#' 
#' #### Function 3
#' #### Description -- This function check whether all elements are positive.
#' #### Input --- target_matrix = the matrix you want to check
#' #### Return -- FALSE if non-positive element exists.
#' 
## ------------------------------------------------------------------------
Positive.Check <- function(target_matrix){
  if(length(which(target_matrix<=0))==0){
    return(TRUE)
  }else{
    return(FALSE)
  }
}


#' 
#' 
#' 
#' #### Function 4
#' #### Description -- This function checks whether positive definite
#' #### Return  -- TRUE if positive definite
#' 
## ------------------------------------------------------------------------
PD.Check.2d <- function(target_matrix){
  if(any(eigen(target_matrix)$values<=0))#check positive
     return(FALSE)
  
  return(TRUE)
}

#check whether [,,k] is non-negative when k is from 1 to n
PD.Check.3d <- function(target_matrix){
  v = dim(target_matrix)[3]
  for(k in 1:v){
    if(PD.Check.2d(target_matrix[,,k])!=TRUE){
      print(paste0(k," dim is not positive definite "))
      return(FALSE)
    }
  }
  return(TRUE)
}



PD.Adjust.2d <- function(target_matrix){
  if(any(eigen(target_matrix)$values<=0)){#not positive definite, change
    new_matrix = as.matrix(nearPD(target_matrix,corr = TRUE,doSym = TRUE, eig.tol = 1e-8/eigen(target_matrix)$values[1])$mat)
    return(new_matrix)
  }#if positive definite, not change
  return(target_matrix)
}



PD.Adjust.3d <- function(target_matrix){
  v = dim(target_matrix)[3]; 
  new_matrix = array(0,dim=dim(target_matrix))
  for(k in 1:v){
    new_matrix[,,k] = PD.Adjust.2d(target_matrix[,,k])
  }
  return(new_matrix)
}



#' 
#' 
#' 
#' #### Function 5 
#' #### Description -- This function uses Basis.Cor to calculates Correlation basis when k=1 to m given rho_matrix
#' #### Input --   rho_matrix = a n by n matrix
#' ####            Basis.Cor = a function defined previously, given rho,k and m return Basis value
#' ####            pd.check = a logical variable, if true then check whether positive definite
#' ####            pos.check = a logical variable, if true then check whether all elements are all positive
#' 
#' $$Basis[i,j,k] =\frac{Beta(k+\rho_{[i,j]}^2,m-k+1)}{Beta(k,m-k+1)}$$
#' 
## ------------------------------------------------------------------------

All.Basis.Cor <- function(rho_matrix,m,positive.check, pd.check, pd.adjust){
  n = nrow(rho_matrix)
  
  # a n by n by m array (lower dimension at start)
  All.Basis.Cor_matrix <- array(1,dim=c(n,n,m))
  
  for(k in 1:m){
    All.Basis.Cor_matrix[,,k] <- sapply(rho_matrix,function(x){Basis.Cor(x,k,m)})
  }
   
  # check whether positive
  if(positive.check)
    print(Positive.Check(All.Basis.Cor_matrix))
  
  # check whether positive definite
  if(pd.check)
    print(PD.Check.3d(All.Basis.Cor_matrix))
  
  # Adjust to make sure Basis matrix are all positive definite
  if(pd.adjust){
    new.All.Basis.Cor_matrix = PD.Adjust.3d(All.Basis.Cor_matrix)
    return(new.All.Basis.Cor_matrix)
  }
  
  
 return(All.Basis.Cor_matrix)
}


#' 

Log.All.Basis.Cor_Rcpp_code = "

  Rcpp::NumericMatrix xcpp(x);
  Rcpp::NumericMatrix result_cpp(result);
  
  int n = xcpp.nrow();
  int m_cpp = Rcpp::as<int>(m);

  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      for (int k = m_cpp; k >= 1; k--){
        /*xcpp[nr * j + i] *= kcpp;*/
        /*xcpp(i,j) *=kcpp;*/
        if(k == m_cpp){
          result_cpp(n * j + i,k-1) = -log1p(double(pow(xcpp(i,j),2)/ k));
        }else{
          result_cpp(n * j + i,k-1) = -log1p(double(pow(xcpp(i,j),2)/ k)) + result_cpp(n * j + i,k);
        }
      }
    }
  }

  //return xcpp;
  "

Log.All.Basis.Cor_Rcpp_func <- cxxfunction(signature(x = "numeric",result = "numeric",m = "numeric"), 
                     body=Log.All.Basis.Cor_Rcpp_code, 
                     plugin="Rcpp")

#' 
## ------------------------------------------------------------------------

Log.All.Basis.Cor <- function(rho_matrix,m,positive.check, pd.check, pd.adjust){
  n = nrow(rho_matrix)
  
  # a n by n by m array (lower dimension at start)
  if(TRUE){
    Log.All.Basis.Cor_matrix <- array(1,dim=c(n,n,m))
    
    
    result <- matrix(0,n^2,m)
    Log.All.Basis.Cor_Rcpp_func(rho_matrix,result,m)
    
    for(k in seq(from=1,to=m,by=1)){
      #Log.All.Basis.Cor_matrix[,,k] <- sapply(rho_matrix,function(x){Log.Basis.Cor(x,k,m)})
      Log.All.Basis.Cor_matrix[,,k] <- matrix(result[,k],nrow = n,ncol=n)
    }
    
    All.Basis.Cor_matrix = exp(Log.All.Basis.Cor_matrix)
    
    # check whether positive
    if(positive.check){
      #print(Positive.Check(All.Basis.Cor_matrix))
    }
    
    # check whether positive definite
    if(pd.check){
      #print(PD.Check.3d(All.Basis.Cor_matrix))
    }
    
    # whether to adjust to make sure pd
    if(pd.adjust){
      new.Log.All.Basis.Cor_matrix = log(PD.Adjust.3d(All.Basis.Cor_matrix))
      return(new.Log.All.Basis.Cor_matrix)
    }
    
    
    
    return(Log.All.Basis.Cor_matrix)
  }
  
  
  
  if(FALSE){
    Log.All.Basis.Cor_matrix <- array(1,dim=c(n,n,m))
    All.Basis.Cor_matrix <- array(1,dim=c(n,n,m))
    
    
    result <- matrix(0,n^2,m)
    Log.All.Basis.Cor_Rcpp_func(rho_matrix,result,m)
    
    for(k in seq(from=1,to=m,by=1)){
      #print(k)
      #Log.All.Basis.Cor_matrix[,,k] <- sapply(rho_matrix,function(x){Log.Basis.Cor(x,k,m)})
      Log.All.Basis.Cor_matrix[,,k] <- matrix(result[,k],nrow = n,ncol=n)
      
      ## make it stable
      tmp = exp(Log.All.Basis.Cor_matrix[,,k])
      eigen.decompose =  eigen(tmp)
      eigen.approx = sapply(eigen.decompose$values, function(x){ifelse(x<10^(-4),10^(-4),x)})
      tmp.approx = eigen.decompose$vectors %*% diag(eigen.approx) %*% t(eigen.decompose$vectors)
      
      All.Basis.Cor_matrix[,,k] = tmp.approx
    }
    
    # check whether positive
    #if(positive.check)
    #  print(Positive.Check(All.Basis.Cor_matrix))
    
    # check whether positive definite
    #if(pd.check)
    #  print(PD.Check.3d(All.Basis.Cor_matrix))
    
    # whether to adjust to make sure pd
    #if(pd.adjust){
    #  new.Log.All.Basis.Cor_matrix = log(PD.Adjust.3d(All.Basis.Cor_matrix))
    #  return(new.Log.All.Basis.Cor_matrix)
    #}
    Log.All.Basis.Cor_matrix = log(All.Basis.Cor_matrix)
    
    return(Log.All.Basis.Cor_matrix)
  }
}





#' 
#' 
#' 
#' 
#' #### Function 6
#' #### Description -- This function computes $v_m$ or say $C(0)$ given Correlation matrix $\Sigma_R$
#' #### Input -- Y = a r by n matrix; Sigma_R = a n by n matrix
#' #### Return 
#' 
#' $$\hat{v}_m=\frac{\sum_{i=1}^ry_i^T[\Sigma_R^{(m)}]^{-1}y_i}{rn}$$
#' 
#' 
#' eigen decomposition : $\Sigma_R^{(m)}=VDV^{-1}$, where V is a orthogonal matrix satisfying $V^{-1}=V^T$ and $D$ is a diagnol matrix with all positive elements
#' 
#' $$y^T[\Sigma_R^{(m)}]^{-1}y = y^T[VDV^{-1}]^{-1}y=y^T[VD^{-1}V^{-1}]y=y^T[VD^{-1}V^T]y=(V^Ty)^TD^{-1}(V^Ty)$$
#' We use $D^{-1,Approx}$ to approximate $D^{-1}$.
#' $$D=[\lambda_1,\lambda_2,\cdots,\lambda_n]$$
#' $$D^{-1,Approx}=[I_{(\lambda_1\geq tol)}\lambda_1^{-1},I_{(\lambda_s\geq tol)}\lambda_2^{-1},\cdots,I_{(\lambda_n\geq tol)} \lambda_n^{-1}]$$
#' 
#' 
#' $$[\Sigma_R^{(m)}]^{-1,Approx} = VD^{-1,Approx}V^T$$
#' $$\hat{v}_m^{Approx}=\frac{\sum_{i=1}^ry_i^T[\Sigma_R^{(m)}]^{-1,Approx}y_i}{rn}$$
#' 
#' 
## ------------------------------------------------------------------------
Cov.Zero <- function(Y,Sigma_R){
     r = nrow(Y); n = ncol(Y)
    
     Sigma_R.inverse = spdinv(Sigma_R)#inverse of symmetric positive definite matrix
  
     sum_vector = apply(Y, 1, function(x){x %*% Sigma_R.inverse %*% x})

     return(mean(sum_vector)/n)
}


Cov.Zero.Approx <- function(Y,Sigma_R,tol){
     r = nrow(Y); n = ncol(Y)
     eigen.decompose =  eigen(Sigma_R)
     #inverse.eigen.approx = sapply(eigen.decompose$values, function(x){ifelse(abs(x)<tol,0,abs(x)^(-1))})
     inverse.eigen.approx = sapply(eigen.decompose$values, function(x){ifelse(x<tol,0,x^(-1))})
     
     Sigma_R.inverse.approx = eigen.decompose$vectors %*% diag(inverse.eigen.approx) %*% t(eigen.decompose$vectors)
     
     sum_vector = apply(Y, 1, function(x){as.numeric(x %*% Sigma_R.inverse.approx %*% x)})
     
     return(mean(sum_vector)/n)
}


Cov.Zero.Approx.ratio <- function(Y,Sigma_R,tol.ratio){
     r = nrow(Y); n = ncol(Y)
     eigen.decompose =  eigen(Sigma_R)
     inverse.eigen.approx = sapply(eigen.decompose$values, 
                                   function(x){
                                     ifelse(abs(x) < abs(max(eigen.decompose$values)*tol.ratio),0,abs(x)^(-1))
                                     })
     
     Sigma_R.inverse.approx = eigen.decompose$vectors %*% 
       diag(inverse.eigen.approx) %*% t(eigen.decompose$vectors)
     
     sum_vector = apply(Y, 1, function(x){as.numeric(x %*% Sigma_R.inverse.approx %*% x)})
     
     return(mean(sum_vector)/n)
}


#' 
#' 
#' 
#' 
#' 
#' #### Function 7
#' #### Description -- This function computes Correlation_matrix based on a symmetric basis
#' #### Input -- Basis = a n by n by m array; w = a 1 by m matrix
#' ####
#' #### Return -- a n by n array
#' $$Basis[i,j,]=Basis[j,i,]$$
#' $$A[i,j]=Basis[i,j,]* w$$
#' 
## ------------------------------------------------------------------------
Product.3d<- function(w,Basis){
  
  Basis_dim = dim(Basis)
  n = Basis_dim[1]; m = Basis_dim[3]
  result_matrix = array(0,dim=c(n,n))
  
  
  result_matrix = apply(Basis, c(1,2), function(x){x %*% w})
  
  return(result_matrix)
}

#' 
#' 
#' 
#' 
## ------------------------------------------------------------------------
Product.3d.Log<- function(w,Log.Basis){
  
  Basis_dim = dim(Log.Basis)
  n = Basis_dim[1]; m = Basis_dim[3]
  result_matrix = array(0,dim=c(n,n))
  
  
  result_matrix = apply(Log.Basis, c(1,2), function(x){sum(exp(x + log(w)))})
  
  return(result_matrix)
}

#' 
#' 
#' 
#' Test Product.3d and Log.Product.3d
## ----eval=FALSE----------------------------------------------------------
## m=20
## set.seed(0)
## w_tilde = runif(m,0,1)
## w_tilde = w_tilde/sum(w_tilde)
## 
## #old
## old.All.Basis.Cor_matrix = All.Basis.Cor(rho_matrix,m,
##                                               positive.check = TRUE, pd.check = TRUE, pd.adjust = TRUE)
## old.Log.All.Basis.Cor_matrix = log(old.All.Basis.Cor_matrix)
## old.result = Product.3d(w_tilde,
##                         Basis = old.All.Basis.Cor_matrix)
## 
## #new
## new.Log.All.Basis.Cor_matrix = Log.All.Basis.Cor(rho_matrix,m,method="Lgamma",
##                                                           positive.check = TRUE, pd.check = TRUE, pd.adjust = TRUE)
## 
## new.All.Basis.Cor_matrix = exp(new.Log.All.Basis.Cor_matrix)
## new.result = Product.3d.Log(w_tilde,
##                             Log.Basis = new.Log.All.Basis.Cor_matrix)
## 
## 
## #
## max(abs(old.All.Basis.Cor_matrix - exp(new.Log.All.Basis.Cor_matrix)))
## 
## eigen(old.result)$values
## eigen(new.result)$values

#' 
#' 
#' 
#' 
#' 
#' #### Function 8
#' #### Description -- This function computes the log det and inverse of a symmetric positive definite matrix
#' #### Input -- A = a symmetric positive definite matrix
#' 
#' $$LogDet(A) = \sum_i(log(\lambda_i))$$
#' 
#' #### return -- a list = LogDet and inv_matrix
## ------------------------------------------------------------------------
library("RcppZiggurat")
library("Rfast")
LogDet.Inv <- function(A){
  
  #===== old =====#
  #logdet = sum(log(eigen(A,symmetric=TRUE)$values)) #
  
  #==== new ===
  logdet = logdet(A) 
  #library("msos")
  return(list(LogDet = logdet))
}

#' 
#' 
#' 
#' #### Function 9
#' #### Description -- This function computes objective function -2l given correlation basis matrix using m parameters  with m-1 degree of freedom
#' #### Input -- w_tilde = weights(length m) whose sum is 1; Y a r by n matrix
#' #### return -- 2 neg log likelihood.
#' 
#' 
#' 
#' $$Basis[i,j,k] =\frac{Beta(k+\rho_{[i,j]}^2,m-k+1)}{Beta(k,m-k+1)}$$
#' 
#' $$\Sigma_R[i,j]=\tilde{w} \times Basis[i,j,]$$
#' $$-2l(\tilde{w}_{1:m})=rnlog(\sum_{i=1}^r y_i^T[\Sigma_R^{(m)}]^{-1}y_i)+rlog|\Sigma_R^{(m)}|-rnlog(rn)+ rn$$
#' 
#' 
#' 
## ------------------------------------------------------------------------

Two.Neg.Log.Likelihood.Use.Basis <- function(w_tilde,Y,All.Basis.Cor_matrix, approx.option,approx.par){
  # if(sum(w_tilde)!=1) stop("sum of w_tilde should be 1")
  
   r = nrow(Y); n = ncol(Y)
   #check sum of w_tilde
   #if(all.equal(sum(w_tilde),1)!=TRUE)
     #stop("the sum of w_tilde is not equal to 1")
   
   #Sigma_R
   Sigma_R = Product.3d(w = w_tilde,Basis = All.Basis.Cor_matrix)
   
   #logdet and inverse
   #temp = LogDet.Inv(Sigma_R)
   #logdet = temp$LogDet
   eigen = eigen(Sigma_R)$values
   logdet = sum(log(eigen[eigen>approx.par & eigen>0]))
   
   
   
   #quad_sum = apply(Y, 1, function(x){x %*% inv_matrix %*% x})
   if(approx.option == 1){#not approx
     quad_sum = Cov.Zero(Y,Sigma_R) * r * n
   }
   if(approx.option == 2){#approx by using tol
     quad_sum = Cov.Zero.Approx(Y,Sigma_R,tol=approx.par) * r * n
   }
   if(approx.option == 3){#approx by using tol
     quad_sum = Cov.Zero.Approx.ratio(Y,Sigma_R,tol.ratio=approx.par) * r * n
   }
   
   
   #result = r*n*log(sum(quad_sum)) + r*logdet
   #result = log(sum(quad_sum)) + logdet/n
   #updated 0902
   result = log(sum(quad_sum)) + logdet/n + 1 - log(r*n)#-2MLE/(r*n)
  return(result)
}
  

  

#' 
#' 
#' 
## ------------------------------------------------------------------------

Two.Neg.Log.Likelihood.Use.Basis.Log <- function(w_tilde,Y,Log.All.Basis.Cor_matrix,approx.option,approx.par){
  # if(sum(w_tilde)!=1) stop("sum of w_tilde should be 1")
  
   r = nrow(Y); n = ncol(Y)
   #check sum of w_tilde
   #if(all.equal(sum(w_tilde),1)!=TRUE)
     #stop("the sum of w_tilde is not equal to 1")
   
   #Sigma_R
   if(all(w_tilde>0)){ #all positive
     Sigma_R = Product.3d.Log(w = w_tilde,Log.Basis = Log.All.Basis.Cor_matrix)
   }else{ #exist negative coef
     Sigma_R = Product.3d(w_tilde,Basis = exp(Log.All.Basis.Cor_matrix))
   }
   
   #logdet and inverse
   #temp = LogDet.Inv(Sigma_R) #only calculate sum(log(eigen))
   #logdet = temp$LogDet
   eigen = eigen(Sigma_R)$values
   logdet = sum(log(eigen[eigen>approx.par & eigen>0]))
   

   #quad_sum = apply(Y, 1, function(x){x %*% inv_matrix %*% x})
   if(approx.option == 1){#not approx
     quad_sum = Cov.Zero(Y,Sigma_R) * r * n
   }
   if(approx.option == 2){#approx by using tol
     quad_sum = Cov.Zero.Approx(Y,Sigma_R,tol=approx.par) * r * n
   }
   if(approx.option == 3){#approx by using tol.ratio
     quad_sum = Cov.Zero.Approx.ratio(Y,Sigma_R,tol.ratio=approx.par) * r * n
   }
   
   
   #result = r*n*log(sum(quad_sum)) + r*logdet - r*n*log(r*n) + r*n
   #result = r*n*log(sum(quad_sum)) + r*logdet 
   #result = log(sum(quad_sum)) + logdet/n #divide by r*n
   #updated 0902
   result = log(sum(quad_sum)) + logdet/n + 1 - log(r*n) #-2MLE/(r*n)
  return(result)
}
  


#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' **Test Part III**
#' 
#' Test objective function
## ----eval=FALSE----------------------------------------------------------
## 
## m=5
## 
## 
## # generate a random vector whose sum is 1
## set.seed(0)
## w_tilde = runif(m,min=0,max=1)
## w_tilde = w_tilde/sum(w_tilde)
## 
## 
## #PD.Check  - new basis
## PD.Check.3d(new.All.Basis.Cor_matrix)
## 
## #PD.Check - old basis
## PD.Check.3d(old.All.Basis.Cor_matrix)
## 
## 
## #objective - new basis
## #Two.Neg.Log.Likelihood.Use.Basis.Log(w_tilde,Y,Log.All.Basis.Cor_matrix = new.Log.All.Basis.Cor_matrix)
## #Two.Neg.Log.Likelihood.Use.Basis(w_tilde, Y, All.Basis.Cor_matrix = new.All.Basis.Cor_matrix)
## 
## #Two.Neg.Log.Likelihood.Use.Basis.Approx.Log(w_tilde, Y, Log.All.Basis.Cor_matrix = new.Log.All.Basis.Cor_matrix,
## #                                            tol = 10^(-8))
## #Two.Neg.Log.Likelihood.Use.Basis.Approx.ratio.Log(w_tilde, Y, Log.All.Basis.Cor_matrix = new.Log.All.Basis.Cor_matrix,
## #                                                  tol.ratio = 10^(-8))
## 
## 
## #objective - old basis
## #Two.Neg.Log.Likelihood.Use.Basis(w_tilde, Y, All.Basis.Cor_matrix = old.All.Basis.Cor_matrix)
## #Two.Neg.Log.Likelihood.Use.Basis.Log(w_tilde, Y,
## #                                     Log.All.Basis.Cor_matrix = old.Log.All.Basis.Cor_matrix)
## 
## 
## 
## #derivative
## #Derivative.Use.Basis(w_tilde,Y, All.Basis.Cor_matrix = All.Basis.Cor_matrix)

#' 
#' 
#' #### Function 10
#' #### Description -- Compute -2l by using a vector of length m-1
#' #### Input -- w_tilde.minus.one = a vector of length m-1
#' 
## ------------------------------------------------------------------------
Two.Neg.Log.Likelihood.Use.Basis.minus.one.par <- function(w_tilde.minus.one,Y,
                                                           All.Basis.Cor_matrix,                      approx.option,approx.par){
  
  w_tilde = c(w_tilde.minus.one,1-sum(w_tilde.minus.one))
  
  result = Two.Neg.Log.Likelihood.Use.Basis(w_tilde,Y, All.Basis.Cor_matrix,
                                            approx.option,approx.par)
   
  return(result)
}

#' 
#' 
## ------------------------------------------------------------------------
Two.Neg.Log.Likelihood.Use.Basis.Log.minus.one.par <- function(w_tilde.minus.one,Y,
                                                               Log.All.Basis.Cor_matrix,approx.option,approx.par){
  
  w_tilde = c(w_tilde.minus.one,1-sum(w_tilde.minus.one))
  
  result = Two.Neg.Log.Likelihood.Use.Basis.Log(w_tilde,Y, Log.All.Basis.Cor_matrix,
                                                approx.option,approx.par)
   
  return(result)
}

#' 
#' 
#' 
#' 
#' 
#' 
#' #### Function 11
#' #### Description -- This function computes the -2l (-2 negative log likelihood)
#' #### Input -- Y = a r by n matrix, Sigma_R = Correlation matrix, Sigma = Covariance matrix
#' $$-2l = rlog|\Sigma^{(m)}|+\sum_{i=1}^r y_i^T[\Sigma^{(m)}]^{-1}y_i$$
#' $$-2l(\tilde{w}_{1:m})=rnlog(\sum_{i=1}^r y_i^T[\Sigma_R^{(m)}]^{-1}y_i)+rlog|\Sigma_R^{(m)}|-rnlog(rn)+ rn$$
#' 
#' 
#' 
## ------------------------------------------------------------------------
Two.Neg.Log.Likelihood.Use.Sigma_R <- function(Y,Sigma_R,approx.option,approx.par){
    
    r = nrow(Y); n = ncol(Y)
   
   #logdet and inverse
   #temp = LogDet.Inv(Sigma_R)
   #logdet = temp$LogDet

   eigen = eigen(Sigma_R)$values
   logdet = sum(log(eigen[eigen>approx.par & eigen>0]))
   
   #quad_sum = apply(Y, 1, function(x){x %*% inv_matrix %*% x})
   if(approx.option == 1){#not approx
     quad_sum = Cov.Zero(Y,Sigma_R) * r * n
   }
   if(approx.option == 2){#approx by using tol
     quad_sum = Cov.Zero.Approx(Y,Sigma_R,tol=approx.par) * r * n
   }
   if(approx.option == 3){#approx by using tol.ratio
     quad_sum = Cov.Zero.Approx.ratio(Y,Sigma_R,tol.ratio=approx.par) * r * n
   }
  
   
   #result = r*n*log(sum(quad_sum)) + r*logdet 
   #result = log(sum(quad_sum)) + logdet/n 
   #updated 0902
   result = log(sum(quad_sum)) + logdet/n + 1 - log(r*n) #-2MLE/(r*n)
  return(result)
  
}


Two.Neg.Log.Likelihood.Use.Sigma <- function(Y,Sigma,approx.option,approx.par){
  
  r = nrow(Y); n = ncol(Y)
   
   #logdet and inverse
  #temp = LogDet.Inv(Sigma)
  #logdet = temp$LogDet
  eigen = eigen(Sigma)$values
  logdet = sum(log(eigen[eigen>approx.par & eigen>0]))
 
   
   #quad_sum = apply(Y,1,function(y){y %*% inv_matrix %*% y})
  
   #quad_sum = apply(Y, 1, function(x){x %*% inv_matrix %*% x})
   if(approx.option == 1){#not approx
     quad_sum = Cov.Zero(Y,Sigma) * r * n
   }
   if(approx.option == 2){#approx by using tol
     quad_sum = Cov.Zero.Approx(Y,Sigma,tol=approx.par) * r * n
   }
   if(approx.option == 3){#approx by using tol.ratio
     quad_sum = Cov.Zero.Approx.ratio(Y,Sigma,tol.ratio=approx.par) * r * n
   }
  #result = r*logdet + log(sum(quad_sum)) 
  #result = logdet/n + log(sum(quad_sum)) 
  #updated 0902
  result = sum(quad_sum/r/n) + logdet/n   #-2MLE/(r*n)
  return(result)
  
}  
  

#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #### Function 12
#' #### Description -- This function calculate $\sum w_tilde - 1$
#' #### Input -- w-tilde is a vector of length m
#' ###        -- approx.option = 1 or 2 or 3
#' 
#' 
#' 
#' 
#' 
#' Define eqfun
## ----eval=FALSE----------------------------------------------------------
## #parameters same as objective function
## Equal.Constrain <- function(w_tilde, Y, All.Basis.Cor_matrix,approx.option,approx.par){
##   return(sum(w_tilde)-1)
## }
## 
## 

#' 
#' 
#' 
## ------------------------------------------------------------------------
#parameters same as objective function
Equal.Constrain.Log <- function(w_tilde, Y, Log.All.Basis.Cor_matrix,
                                approx.option,approx.par){
  return(sum(w_tilde)-1)
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
#' 
#' 
#' #### Description -- This function computes the derivative of objective function with respect to w_tilde_k
#' #### Input -- 
#' $$
#'    \frac{dL}{d\tilde{w}_k}  = -\frac{rn}{\sum_{i=1}^m y_k [\Sigma_R^{(m)}]^{-1}y_k}tr([\Sigma_R^{(m)}]^{-1}(\sum_{i=1}^r y_iy_i^T) [\Sigma_R^{(m)}]^{-1}A_k)+r tr([\Sigma_R^{(m)}]^{-1}A_k)=rtr([\Sigma_R^{(m)}]^{-1}A_k(-\frac{n\Sigma_R^{-1}(\sum_{i=1}^ry_i y_i^T)}{\sum_{i=1}^r y_i^T [\Sigma_R]^{-1}y_i}+I_{n\times n}))
#' $$
#' $$A_k = Basis[,,k]$$
#' $$\Sigma_R^{(m)}=\sum_{k=1}^m\tilde{w}_k A_k$$
#' 
#' Denote 
#' $$\frac{dL}{d\tilde{w}_k}  = rtr(A_k(-\frac{n\Sigma_R^{-1}(\sum_{i=1}^ry_i y_i^T)}{tr(\Sigma_R^{-1}(\sum_{i=1}^ry_i y_i^T))}+I_{n\times n})[\Sigma_R^{(m)}]^{-1})$$
#' #### Return = a vector of length m
#' 
## ------------------------------------------------------------------------
library("psych") # for tr function

Derivative.Use.Basis <- function(w_tilde,Y,All.Basis.Cor_matrix){
    
   r = nrow(Y); n = ncol(Y); m = dim(All.Basis.Cor_matrix)[3]
  
   #check sum of w_tilde
   if(all.equal(sum(w_tilde),1)!=TRUE)
     stop("the sum of w_tilde is not equal to 1")
   
   #Sigma_R
   Sigma_R = Product.3d(w_tilde,All.Basis.Cor_matrix)
   
   #inverse
   inv_matrix = LogDet.Inv(Sigma_R)$inv_matrix
   
   prod_sum = matrix(apply(apply(Y, 1, function(x){x %*% t(x)}),1,sum),nrow=n,ncol=n)
   
   #a = matrix(0,nrow=n,ncol=n)
   #for(i in 1:r){
   #   a = a + Y[i,] %*%  t(Y[i,])
   #}
   
   temp1 = inv_matrix %*% prod_sum
   temp2 = (-n*temp1/tr(temp1) + diag(n)) %*% inv_matrix
   
   result = numeric(m)
   for(k in 1:m){
     result[k] = tr(All.Basis.Cor_matrix[,,k] %*% temp2)
   }
  return(result)
   
}
  
  

#' 
#' 
#' 
#' #### Input -- w_tilde.minus.one = weights of length m-1
#' #### Return -- 
#' 
## ------------------------------------------------------------------------
Derivative.Use.Basis.minus.one <- function(w_tilde.minus.one,Y,All.Basis.Cor_matrix){
  m = length(w_tilde.minus.one)+1
  w_tilde = c(w_tilde.minus.one,1-sum(w_tilde.minus.one))
  result = Derivative.Use.Basis(w_tilde,Y,All.Basis.Cor_matrix)#length m
  result.new = result[1:(m-1)] - result[m]
  return(result.new)
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
#' 
#' To estimate $\tilde{w}_k$
#' $$min -2l(\tilde{w}_{1:m})=rnlog(\sum_{i=1}^r y_i^T[\Sigma_R^{(m)}]^{-1}y_i)+rlog|\Sigma_R^{(m)}|-rnlog(rn)+ rn$$
#' 
#' $$0\leq \tilde{w}_k \leq 1$$
#' 
#' $$\sum_{k=1}^m \tilde{w}_k=1$$
## ------------------------------------------------------------------------
library(Rsolnp) #solnp

#' 
#' 
#' 
#' **Method 1**
#' 
#' Use solnp function: solnp(pars,fun,eqfun=,eqB=,LB,UB)
#' converge=0 means sucess.
#' the number of parameters is m.
#' 
#' Comments: In updating parameters, the constraints may not be met while the result is sure to meet the constraints.
#' 
#' 
#' 
#' 
#' 
#' 
## ------------------------------------------------------------------------
#parameters same as objective function
Inequal.Constrain.minus.one <- function(w_tilde.minus.one, Y, All.Basis.Cor_matrix){
 
  return(sum(w_tilde.minus.one)-1) # sum of m-1 weights <=1
}


Derivative.Inequal.Constrain.minus.one <- function(w_tilde.minus.one, Y, All.Basis.Cor_matrix){
  return(rep(1,length(w_tilde.minus.one)))
}


#' 
#' 
#' 
#' 
#' 

