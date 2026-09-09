library(data.table);library(ggplot2)
sim<-fread('psi_validation/data/psi05_stochastic_results.csv')
curves<-fread('psi_validation/data/psi05_nonlinear_tail_curves.csv')
dir.create('psi_validation/figures',FALSE,TRUE)
theme_nl<-theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),
  strip.background=element_rect(fill='grey88'),plot.title=element_text(size=15),
  plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0),legend.position='bottom')
draw_leading<-function(probability,rho_value,theta_value){
  a<-copy(curves[rho==rho_value&theta==theta_value]);s<-copy(sim[rho==rho_value&theta==theta_value])
  if(probability=='unconditional'){
    a[,value:=leading*(1-1/R0)];s[,`:=`(estimate=P_unconditional,low=uncond_low,high=uncond_high)]
    ttl<-yl<-'Unconditional persistence probability'
  }else{
    a[,value:=leading];s[,`:=`(estimate=P_conditional,low=cond_low,high=cond_high)]
    ttl<-'Conditional persistence probability';yl<-'Persistence probability, conditional on not fizzling'
  }
  p<-ggplot()+geom_line(data=a,aes(R0-1,value),colour='#009E73',linewidth=.9,na.rm=TRUE)+
    geom_errorbar(data=s,aes(R0-1,ymin=low,ymax=high),width=0,colour='grey25',linewidth=.42,na.rm=TRUE)+
    geom_point(data=s,aes(R0-1,estimate),colour='#2C7FB8',size=2,na.rm=TRUE)+
    facet_wrap(~K,ncol=2,labeller=labeller(K=function(x)paste0('K = ',x)))+
    scale_x_log10(breaks=c(.05,.1,.2,.5,1,2,5))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),expand=expansion(mult=c(.01,.03)))+
    labs(title=ttl,subtitle=sprintf('Nonlinear-tail action amplitude; psi = 0.5, rho = %.2f, theta = %g; I0 = 1',rho_value,theta_value),
      x=expression(R[0]-1~'(log scale)'),y=yl,
      caption='Blue points and grey bars: reused CTMC estimates and 95% Wilson intervals; green: nonlinear-tail Laplace prediction. Curves stop where no admissible descending-tail initialization exists.')+theme_nl
  if(probability=='unconditional'){
    ref<-unique(a[,.(K,R0)]);ref[,value:=1-1/R0]
    p<-p+geom_line(data=ref,aes(R0-1,value,group=K),colour='grey60',linetype='dashed',linewidth=.55)
  };p
}
write_leading<-function(probability,file){cairo_pdf(file,10.5,8.4,onefile=TRUE)
  for(rho in c(.01,.02,.05,.10))for(theta in c(0,.5,1))print(draw_leading(probability,rho,theta));dev.off()}
write_leading('unconditional','psi_validation/figures/fig12_psi05_nonlinear_tail_unconditional.pdf')
write_leading('conditional','psi_validation/figures/fig13_psi05_nonlinear_tail_conditional.pdf')

draw_next<-function(rho_value,theta_value){
  a<-melt(curves[rho==rho_value&theta==theta_value],
    id.vars=c('rho','theta','K','R0'),measure.vars=c('leading','next_order'),
    variable.name='method',value.name='value')
  a[,method:=factor(method,c('leading','next_order'),c('Leading nonlinear tail','Next Laplace correction'))]
  s<-sim[rho==rho_value&theta==theta_value]
  ggplot()+geom_line(data=a,aes(R0-1,value,colour=method),linewidth=.9,na.rm=TRUE)+
    geom_errorbar(data=s,aes(R0-1,ymin=cond_low,ymax=cond_high),width=0,colour='grey25',linewidth=.42,na.rm=TRUE)+
    geom_point(data=s,aes(R0-1,P_conditional),colour='#2C7FB8',size=2,na.rm=TRUE)+
    facet_wrap(~K,ncol=2,labeller=labeller(K=function(x)paste0('K = ',x)))+
    scale_x_log10(breaks=c(.05,.1,.2,.5,1,2,5))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),expand=expansion(mult=c(.01,.03)))+
    scale_colour_manual(values=c('Leading nonlinear tail'='#009E73','Next Laplace correction'='#D55E00'))+
    labs(title='Conditional persistence probability',subtitle=sprintf('Nonlinear-tail comparison; psi = 0.5, rho = %.2f, theta = %g; I0 = 1',rho_value,theta_value),
      x=expression(R[0]-1~'(log scale)'),y='Persistence probability, conditional on not fizzling',colour=NULL,
      caption='Blue points and grey bars: reused CTMC estimates and 95% Wilson intervals. Green: leading nonlinear tail; orange: exp(-rho*cL) correction.')+theme_nl
}
cairo_pdf('psi_validation/figures/fig15_psi05_nonlinear_tail_next_stochastic_curves.pdf',10.5,8.4,onefile=TRUE)
for(rho in c(.01,.02,.05,.10))for(theta in c(0,.5,1))print(draw_next(rho,theta));dev.off()
cat('nonlinear-tail Figures 12, 13, and 15 complete\n')
