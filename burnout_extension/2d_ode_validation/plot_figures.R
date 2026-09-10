library(data.table);library(ggplot2)
sim<-fread('psi_validation/data/psi05_stochastic_results.csv')
curves<-fread('2d_ode_validation/data/psi05_2d_ode_curves.csv')[status=='OK']
dir.create('2d_ode_validation/figures',FALSE,TRUE)
theme_2d<-theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),
  strip.background=element_rect(fill='grey88'),plot.title=element_text(size=15),
  plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0),
  legend.position='bottom')
cols_method<-c('Numerical finite-horizon Kendall'='#009E73',
  'Closed trough saddle'='#CC79A7')

draw_probability<-function(kind,rho_value,theta_value){
  if(kind=='conditional'){
    measures<-c('P_kendall','P_trough');yl<-'Persistence probability, conditional on establishment'
    ttl<-'Conditional persistence probability'
    s<-copy(sim[rho==rho_value&theta==theta_value])
    s[,`:=`(estimate=P_conditional,low=cond_low,high=cond_high)]
  }else{
    measures<-c('P_kendall_unconditional','P_trough_unconditional');yl<-'Unconditional persistence probability'
    ttl<-'Unconditional persistence probability'
    s<-copy(sim[rho==rho_value&theta==theta_value])
    s[,`:=`(estimate=P_unconditional,low=uncond_low,high=uncond_high)]
  }
  a<-melt(copy(curves[rho==rho_value&theta==theta_value]),
    id.vars=c('rho','theta','K','R0'),measure.vars=measures,
    variable.name='method',value.name='value')
  a[,method:=factor(method,measures,c('Numerical finite-horizon Kendall','Closed trough saddle'))]
  p<-ggplot()+geom_line(data=a,aes(R0-1,value,colour=method),linewidth=.9,na.rm=TRUE)+
    geom_errorbar(data=s,aes(R0-1,ymin=low,ymax=high),width=0,colour='grey25',linewidth=.42)+
    geom_point(data=s,aes(R0-1,estimate),colour='#2C7FB8',size=2)+
    facet_wrap(~K,ncol=2,labeller=labeller(K=function(x)paste0('K = ',x)))+
    scale_x_log10(breaks=c(.05,.1,.2,.5,1,2,5))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),expand=expansion(mult=c(.01,.03)))+
    scale_colour_manual(values=cols_method)+
    labs(title=ttl,subtitle=sprintf('Full 2D ODE; psi = 0.5, rho = %.2f, theta = %g; I0 = 1',rho_value,theta_value),
      x=expression(R[0]-1~'(log scale)'),y=yl,colour=NULL,
      caption='Blue points and grey bars: reused CTMC estimates and 95% Wilson intervals. Green: direct finite-horizon Kendall quadrature. Purple: closed physical-trough saddle formula.')+theme_2d
  if(kind=='unconditional'){
    ref<-unique(a[,.(K,R0)]);ref[,value:=1-1/R0]
    p<-p+geom_line(data=ref,aes(R0-1,value,group=K),colour='grey60',linetype='dashed',linewidth=.55)
  }
  p
}
write_pages<-function(kind,file){
  cairo_pdf(file,10.5,8.4,onefile=TRUE)
  for(rho in c(.01,.02,.05,.10))for(theta in c(0,.5,1))print(draw_probability(kind,rho,theta))
  dev.off()
}
write_pages('unconditional','2d_ode_validation/figures/fig12_psi05_2d_ode_unconditional.pdf')
write_pages('conditional','2d_ode_validation/figures/fig13_psi05_2d_ode_conditional.pdf')

cmp<-curves[K==1000&R0%in%c(1.5,2,3,6)]
mae<-mean(abs(cmp$P_trough-cmp$P_kendall))
p<-ggplot(cmp,aes(P_kendall,P_trough,colour=factor(rho),shape=factor(theta)))+
  geom_abline(slope=1,intercept=0,colour='grey55',linetype='dashed')+
  geom_point(size=2.3)+coord_equal(xlim=c(0,1),ylim=c(0,1))+
  scale_colour_manual(values=c('0.01'='#0072B2','0.02'='#56B4E9','0.05'='#D55E00','0.1'='#009E73'),name=expression(rho))+
  scale_shape_manual(values=c('0'=16,'.5'=17,'0.5'=17,'1'=15),name=expression(theta))+
  labs(title='Closed trough formula versus numerical Kendall quadrature',
    subtitle=sprintf('48-point theory grid; psi = 0.5, K = 1000; mean absolute difference = %.5f',mae),
    x='Numerical finite-horizon Kendall probability',y='Closed trough-saddle probability',
    caption='The diagonal is equality. This comparison isolates the Gaussian trough reduction\non the same full nonlinear deterministic background.')+
  theme_2d
cairo_pdf('2d_ode_validation/figures/fig15_psi05_2d_ode_trough_vs_kendall.pdf',7.4,7.0,onefile=FALSE)
print(p);dev.off()
cat('2D ODE Figures 12, 13, and 15 complete\n')
