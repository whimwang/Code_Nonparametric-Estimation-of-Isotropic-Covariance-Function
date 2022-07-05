#' ---
#' title: "NP_Choi_Estimation"
#' author: "Whim Wang"
#' date: "7/24/2020"
#' output: html_document
#' ---
#' 
#' First,calculate empirical estimator of $C(h_i)$
#' Second, calculate $f_j^q(h_i^2)$ for $j=1,2\cdots,m+p$
#' 
#' 
#' 
#' Calculate statistics for observations for r independent replications of n locations
#' @param coord matrix,n by d
#' @param observation Y,r by n
#' @return list,sum_covhat $\hat{C}(h_i)\cdot N(h_i)$;n_h $N(h_i)$;dist $h_i$
#' 
## -------------------------------------------------------------------------------------------------
Choi_Empirical_Cov <- function(Y,coord_matrix,dist_type='norm2'){
  
  r = nrow(Y);n=ncol(Y)
  if(nrow(coord_matrix)!=n){
    stop('col of Y is not equal to row of coord matrix')
  }
  
  
  #assume constant mean
  #Y_mean = mean(Y)
  #Y = Y - Y_mean
  
  #nonzero dist=n*(n-1)/2 + zero dist=n
  result = matrix(0,nrow=n*(n-1)/2+1,ncol=3)
  
  ## nonzero dist
  start_idx = 1
  for(i in 1:(n-1)){
    #ith sample with (i+1):n sample
    end_idx = start_idx + (n-i)-1 #save result for (i,i+1) to (i,n),n-i in total
    #dist
    if(dist_type=='norm2'){
      dist = sapply(seq(from=i+1,to=n,by=1), function(idx){norm(coord_matrix[i,]-coord_matrix[idx,],type='2')}) 
    }else if(dist_type == 'distCosine'){
      dist = sapply(seq(from=i+1,to=n,by=1), function(idx){distm(coord_matrix[i,],coord_matrix[idx,],fun=distCosine)})
    }
    result[start_idx:end_idx,1] = dist
    #sum of cov_hat
    cov_hat = sapply(seq(from=i+1,to=n,by=1), function(idx){sum(Y[,i] * Y[,idx])})
    result[start_idx:end_idx,2] = cov_hat
    #N_h
    num = rep(r,end_idx-start_idx+1)
    result[start_idx:end_idx,3] = num
    #
    #print(paste0('start: ',start_idx,' end: ',end_idx))
    start_idx = end_idx+1
  }
  
  ## zero dist
  result[start_idx,1] = 0
  result[start_idx,2] = sum(sapply(seq(from=1,to=n,by=1),function(idx){sum(Y[,idx]^2)}))
  result[start_idx,3] = n * r
  
  ## check unique dist
  if(length(unique(result[,1])) == nrow(result)){
    return(list(dist = result[,1],
              sum_covhat = result[,2],
              n_h = result[,3]))
  }
  
  ## If not unique dist: need a summary
  unique_dist = unique(result[,1])
  summary_result = matrix(0,nrow=length(unique_dist),ncol=3)
  for(i in 1:length(unique_dist)){
    idx = which(result[,1]==unique_dist[i])
    if(length(idx)!=1){
      summary_result[i,] = c(unique_dist[i],apply(result[idx,2:3],2,sum))
    }else{
      summary_result[i,] = result[idx,]
    }
  }
  
  
  return(list(dist = summary_result[,1],
              sum_covhat = summary_result[,2],
              n_h = summary_result[,3]))
}


#' 
#' 
#' 
#' 
## -------------------------------------------------------------------------------------------------
##knitr::purl('NP_Choi_Estimation.Rmd',documentation = 2L)

#' 
