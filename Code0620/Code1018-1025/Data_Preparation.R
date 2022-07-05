







################ Data Preparation ################
#'
#'

Data_Preparation <- function(option, par.cov, 
                             d,N,r,n,epsilon,option.range,seed=0){
  
  rho_support = Support.function(par.cov=c(1,par.cov[2:length(par.cov)]), option, epsilon = epsilon)
  l_max = l.range.function(rho_support,d, option.range)
  
  
  coord_matrix = Coordinates_Generation(d,n,-l_max,l_max,seed=seed)
  
  
  rho_matrix.orig = Distance(coord_matrix)
  
  Y_matrix = array(0,dim=c(N,r,n))
  for(i in 1:N)
  {
    temp = Data_Generation(Cov_function,rho_matrix.orig,r,option,par.cov,seed=i*1000)
    Y_matrix[i,,] = temp$Y
    if(i == 1){
      True.Sigma_matrix = temp$Sigma_matrix
    }
  }
  
  Output.list <- list(Y_matrix = Y_matrix,
                      rho_matrix.orig = rho_matrix.orig,
                      coord_matrix = coord_matrix)
  
  Hidden.list <- list(True.Sigma_matrix = True.Sigma_matrix,
                      rho_support = rho_support,
                      l_max = l_max)
  
  return(list(Output.list = Output.list, Hidden.list = Hidden.list))
}





Data_Preparation_Fixed_range <- function(option, par.cov, 
                             d,N,r,n,l_max,seed=0){
  
  coord_matrix = Coordinates_Generation(d,n,0,l_max,seed=seed)
  
  
  rho_matrix.orig = Distance(coord_matrix)
  
  Y_matrix = array(0,dim=c(N,r,n))
  for(i in 1:N)
  {
    temp = Data_Generation(Cov_function,rho_matrix.orig,r,option,par.cov,seed=i*1000)
    Y_matrix[i,,] = temp$Y
    if(i == 1){
      True.Sigma_matrix = temp$Sigma_matrix
    }
  }
  
  Output.list <- list(Y_matrix = Y_matrix,
                      rho_matrix.orig = rho_matrix.orig,
                      coord_matrix = coord_matrix)
  
  Hidden.list <- list(True.Sigma_matrix = True.Sigma_matrix,
                      l_max = l_max)
  
  return(list(Output.list = Output.list, Hidden.list = Hidden.list))
}

