source('validation/R/theory.R')
stopifnot(requireNamespace('data.table',quietly=TRUE),
          requireNamespace('ggplot2',quietly=TRUE))
library(data.table)
library(ggplot2)

sim_file <- 'validation/data/stochastic_scan_results.csv'
sim <- fread(sim_file)

dir.create('validation/figures',FALSE,TRUE)

# Use the same R0-1 curves, K facets, stochastic points, and Wilson intervals
# as Figures 9--11.
analytic <- fread('validation/data/stochastic_scan_analytic.csv')
analytic[,BI:=vapply(seq_len(.N),function(i)
  bi_quantities(R0[i],rho[i],theta[i],K[i])$P_conditional,0.)]
draw_bi_page <- function(probability=c('unconditional','conditional'),rho_value,theta_value,
                         curve_data=analytic,simulation_data=sim,large_K=FALSE) {
  probability <- match.arg(probability)
  a <- copy(curve_data[rho==rho_value & theta==theta_value,.(rho,theta,K,R0,value=BI)])
  s <- copy(simulation_data[rho==rho_value & theta==theta_value])
  if(probability=='unconditional') {
    a[,value:=value*(1-1/R0)]
    s[,`:=`(estimate=P_unconditional,low=uncond_low,high=uncond_high)]
    title <- 'Unconditional persistence probability'
    ylabel <- 'Unconditional persistence probability'
  } else {
    s[,`:=`(estimate=P_conditional,low=cond_low,high=cond_high)]
    title <- 'Conditional persistence probability'
    ylabel <- 'Persistence probability, conditional on not fizzling'
  }
  q <- ggplot()+
    geom_line(data=a,aes(R0-1,value),colour='#009E73',linewidth=.9)+
    geom_errorbar(data=s,aes(R0-1,ymin=low,ymax=high),width=0,
      colour='grey25',linewidth=.42)+
    geom_point(data=s,aes(R0-1,estimate),colour='#2C7FB8',size=2)+
    facet_wrap(~K,ncol=2,labeller=labeller(K=function(x) paste0('K = ',
      if(large_K) format(as.numeric(x),scientific=TRUE) else x)))+
    scale_x_log10(breaks=c(.03,.05,.1,.2,.5,1,2,5),
      labels=c('0.03','0.05','0.10','0.20','0.50','1','2','5'))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),
      expand=expansion(mult=c(.01,.03)))+
    labs(title=title,
      subtitle=sprintf('Boundary-layer-independent leading order; rho = %.2f, theta = %g; I0 = 1',rho_value,theta_value),
      x=expression(R[0]-1~'(log scale)'),y=ylabel,
      caption=paste0('Blue points and grey bars are the existing stochastic estimates and 95% Wilson intervals; ',
        'the green curve is BI. BI remains defined through old NO_ENTRY regions.'))+
    theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),
      strip.background=element_rect(fill='grey88'),plot.title=element_text(size=15),
      plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0))
  if(probability=='unconditional') {
    ref <- unique(a[,.(K,R0)]);ref[,value:=1-1/R0]
    q <- q+geom_line(data=ref,aes(R0-1,value,group=K),colour='grey60',
      linewidth=.55,linetype='dashed')
  }
  q
}

large_analytic <- unique(fread('validation/data/fig11_scan_analytic.csv')[,
  .(rho,theta,K,R0)])
large_analytic[,BI:=vapply(seq_len(.N),function(i)
  bi_quantities(R0[i],rho[i],theta[i],K[i])$P_conditional,0.)]
large_sim <- fread('validation/data/fig11_scan_results.csv')
stopifnot(nrow(large_analytic)==2172L,nrow(large_sim)==168L)

write_bi_pdf <- function(probability,file) {
  cairo_pdf(file,width=10.5,height=8.4,onefile=TRUE)
  for(rho_value in c(.01,.02,.05,.10)) for(theta_value in c(0,.5,1))
    print(draw_bi_page(probability,rho_value,theta_value))
  for(theta_value in c(0,.5,1))
    print(draw_bi_page(probability,.01,theta_value,large_analytic,large_sim,TRUE))
  dev.off()
}
write_bi_pdf('unconditional','validation/figures/fig12_BI_unconditional_stochastic_validation.pdf')
write_bi_pdf('conditional','validation/figures/fig13_BI_conditional_stochastic_validation.pdf')

cat('BI conditional and unconditional validation PDFs complete\n')
