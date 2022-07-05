fig_prefix = 'Year40_Station_189'

input_data.path = "/Users/wangyiming/Documents/2019 Research Daily/research_code/Code_2020/Updated0620/Code0620/input_real_data/Year40_Station189_46_30_80_100/"
load(paste0(input_data.path,"Year40_Station189_46_30_80_100.RData"))
lon_lat_data = input_data$lon_lat_data
x_y_scale = input_data$x_y_scale
observation_residual = input_data$observation_residual
dist_select = input_data$dist_select
cov.distort = 16
dist.range_parameter = 10^5

setwd("/Users/wangyiming/Documents/2019 Research Daily/research_code/Code_2020/Result_1022/Choi_real_data/m_7_p_3_cov.distort_16_dist.range_parameter_1e+05_Oct22_10AM")

load('Choi_NP_change_weight_result.RData')
load('Choi_result_orig_converted.RData')
load('Choi_parametric_change_weight.RData')

paste0(round(sapply(Choi_parametric_change_weight,function(x){x$deviance}),digits = 3),collapse = ' & ')
paste0(round(sapply(Choi_parametric_change_weight,function(x){x$time}),digits = 2),collapse = ' & ')
paste0(round(sapply(Choi_parametric_change_weight,function(x){x$par[1]} * cov.distort),digits = 2),collapse = ' & ')
paste0(round(sapply(Choi_parametric_change_weight,function(x){x$par[2]} * cov.distort),digits = 2),collapse = ' & ')
paste0(round(sapply(Choi_parametric_change_weight,function(x){x$par[3]} * dist.range_parameter/10^5),digits = 2),collapse = ' & ')

sapply(Choi_parametric_change_weight,function(x){x$par})



optout_BsplineChoi_change_weight = list(par = Choi_NP_change_weight_result$par)
ChoiFitted_arg = Choi_NP_change_weight_result$ChoiFitted_arg

f_cov_est_Choi_change_weight <- function(rho){
  Cov_BsplineChoi(rho/dist.range_parameter,beta=optout_BsplineChoi_change_weight$par,p=ChoiFitted_arg$p,m=ChoiFitted_arg$m) * cov.distort
}

f_cor_est_Choi_change_weight <- function(rho){
  f_cov_est_Choi_change_weight(rho)/f_cov_est_Choi_change_weight(0)
}

f_semivariog_est_Choi_change_weight <- function(rho){
  f_cov_est_Choi_change_weight(0) - f_cov_est_Choi_change_weight(rho) + Choi_NP_change_weight_result$nugget * cov.distort
}





d_seq = seq(from=0,to=max(dist_select)*1.2,length.out = 100)
Choi_change_seq_plot = list(d_seq = d_seq,
                         cov_seq = sapply(d_seq, f_cov_est_Choi_change_weight),
                         cor_seq = sapply(d_seq, f_cor_est_Choi_change_weight),
                         semivariog_seq = sapply(d_seq, f_semivariog_est_Choi_change_weight))

save(Choi_change_seq_plot,file='Choi_change_seq_plot.RData')

plot(Choi_change_seq_plot$d_seq/10^5,Choi_change_seq_plot$cov_seq)
plot(Choi_change_seq_plot$d_seq/10^5,Choi_change_seq_plot$semivariog_seq)


if(TRUE){
#Matern
Matern_change_optout = list(d_seq = d_seq,
                     cov_seq = sapply(d_seq/dist.range_parameter, Cov_function,option='Matern',par.cov = Choi_parametric_change_weight$optout_Matern_change_weight$par[-1]) * cov.distort)

Matern_change_optout$cor_seq = Matern_change_optout$cov_seq/Choi_parametric_change_weight$optout_Matern_change_weight$par[2]/cov.distort
Matern_change_optout$semivariog_seq = Choi_parametric_change_weight$optout_Matern_change_weight$par[2] * cov.distort - Matern_change_optout$cov_seq + Choi_parametric_change_weight$optout_Matern_change_weight$par[1] * cov.distort

#Cauchy
Cauchy_change_optout = list(d_seq = d_seq,
                            cov_seq = sapply(d_seq/dist.range_parameter, Cov_function,option='Cauchy',par.cov = Choi_parametric_change_weight$optout_Cauchy_change_weight$par[-1]) * cov.distort)

Cauchy_change_optout$cor_seq = Cauchy_change_optout$cov_seq/Choi_parametric_change_weight$optout_Cauchy_change_weight$par[2]/cov.distort
Cauchy_change_optout$semivariog_seq = Choi_parametric_change_weight$optout_Cauchy_change_weight$par[2] * cov.distort - Cauchy_change_optout$cov_seq + Choi_parametric_change_weight$optout_Cauchy_change_weight$par[1] * cov.distort

#Gaussian
Gaussian_change_optout = list(d_seq = d_seq,
                            cov_seq = sapply(d_seq/dist.range_parameter, Cov_function,option='Gaussian',par.cov = Choi_parametric_change_weight$optout_Gaussian_change_weight$par[-1]) * cov.distort)

Gaussian_change_optout$cor_seq = Gaussian_change_optout$cov_seq/Choi_parametric_change_weight$optout_Gaussian_change_weight$par[2]/cov.distort
Gaussian_change_optout$semivariog_seq = Choi_parametric_change_weight$optout_Gaussian_change_weight$par[2] * cov.distort - Gaussian_change_optout$cov_seq + Choi_parametric_change_weight$optout_Gaussian_change_weight$par[1] * cov.distort

#GenCauchy
GenCauchy_change_optout = list(d_seq = d_seq,
                              cov_seq = sapply(d_seq/dist.range_parameter, Cov_function,option='GenCauchy',par.cov = Choi_parametric_change_weight$optout_GenCauchy_change_weight$par[-1]) * cov.distort)

GenCauchy_change_optout$cor_seq = GenCauchy_change_optout$cov_seq/Choi_parametric_change_weight$optout_GenCauchy_change_weight$par[2]/cov.distort
GenCauchy_change_optout$semivariog_seq = Choi_parametric_change_weight$optout_GenCauchy_change_weight$par[2] * cov.distort - GenCauchy_change_optout$cov_seq + Choi_parametric_change_weight$optout_GenCauchy_change_weight$par[1] * cov.distort



Choi_change_parametric_seq_plot = list(Matern_change_optout = Matern_change_optout,
                                       Cauchy_change_optout = Cauchy_change_optout,
                                       Gaussian_change_optout = Gaussian_change_optout,
                                       GenCauchy_change_optout = GenCauchy_change_optout)
save(Choi_change_parametric_seq_plot,file='Choi_change_parametric_seq_plot.RData')

}

