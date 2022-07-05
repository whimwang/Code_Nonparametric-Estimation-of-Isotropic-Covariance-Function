Write.List.result <- function(all.best.seed.list, task.name, whether.best,text=NULL){
  
  N = length(all.best.seed.list)
 
  est.fvalue.vector = sapply(all.best.seed.list, function(x){return(x$est.V)})
  transform.ratio.vector = sapply(all.best.seed.list, function(x){return(x$transform.ratio)})
  
  m.vector = sapply(all.best.seed.list, function(x){return(x$m)})
  
  est.par_matrix = matrix(NA, nrow = N, ncol=max(m.vector))
  for(k in 1:N){
    est.par_matrix[k,1:m.vector[k]] = all.best.seed.list[[k]]$est.par
  }
  rownames(est.par_matrix)=paste("m=",m.vector)
  
  
  
  
  if(whether.best){
  #ratio
  ratio.file = paste0(task.name,".","best",".transform.ratio.txt")
  write.table(matrix(transform.ratio.vector,nrow = N,ncol=1),
              file=ratio.file,
              quote = FALSE,col.names = FALSE,
              row.names = paste0("Simulation",1:N))
  }
  
  #est.V
  est.V.vector = sapply(all.best.seed.list, function(x){return(x$est.V)})
  V.file = paste0(task.name,".",text,".V.txt")
  write.table(matrix(est.V.vector,nrow = N,ncol=1),
              file=V.file,
              quote = FALSE,col.names = FALSE,
              row.names = paste0("Simulation",1:N))
  #m
  m.file = paste0(task.name,".",text,".m.txt")
  write.table(matrix(m.vector,nrow = N,ncol=1),
              file=m.file,
              quote = FALSE,col.names = FALSE,
              row.names = paste0("Simulation",1:N))
  
  #est.par_matrix
  est.par.file = paste0(task.name,".",text,".est.par.txt")
  write.table(est.par_matrix,
              file=est.par.file,
              quote = FALSE,col.names = FALSE)
  
  
}