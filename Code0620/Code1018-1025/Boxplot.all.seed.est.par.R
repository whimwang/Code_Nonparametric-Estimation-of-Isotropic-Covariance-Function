
Boxplot.all.seed.est.par <- function(all.est.par_matrix,fig.name){
  colnames(all.est.par_matrix) = paste0("m",1:ncol(all.est.par_matrix))
  
  temp.matrix = melt(all.est.par_matrix)[,-1]
  colnames(temp.matrix) = c("variable","value")
  
  sp <- ggplot(data = temp.matrix, aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))+ylim(0,1)+
    theme(plot.title = element_text(size = 10), axis.text.x = element_text(angle = 90,size=5),legend.position = "none")
  sp + labs(x = "estimated coef",
            title = paste0("d:",d," n:",n," m:",m," n_seed:",n_seed," option:",option," N:",N,
                           " par.cov:","(",paste0(par.cov,collapse = ","),")"))
  ggsave(filename = fig.name,width=7)
}