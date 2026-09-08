source('psi_validation/R/theory_psi.R')
stopifnot(requireNamespace('data.table',quietly=TRUE),requireNamespace('ggplot2',quietly=TRUE))
library(data.table);library(ggplot2)
dir.create('psi_validation/figures',FALSE,TRUE)
sim<-fread('psi_validation/data/psi05_stochastic_results.csv')
stopifnot(all(sim$attempts>=3000),all(sim$unresolved==0))
psi<-0.5
dense_for<-function(scan_value){
  template<-unique(sim[scan==scan_value,.(rho,theta,K)])
  rv<-if(scan_value=='exact')1+exp(seq(log(.05),log(5),length.out=101)) else
    1+exp(seq(log(.03),log(5),length.out=101))
  g<-template[,.(R0=rv),by=.(rho,theta,K)];keys<-unique(g[,.(R0,theta)])
  keys[,`:=`(C=NA_real_,D=NA_real_)]
  for(i in seq_len(nrow(keys))){cd<-CD_psi(keys$R0[i],keys$theta[i],psi);keys[i,`:=`(C=cd['C'],D=cd['D'])]}
  g<-merge(g,keys,by=c('R0','theta'))
  v<-t(vapply(seq_len(nrow(g)),function(i){
    q1<-bi_quantities_psi(g$R0[i],g$rho[i],g$theta[i],psi,g$K[i],next_order=FALSE)
    q2<-bi_quantities_psi(g$R0[i],g$rho[i],g$theta[i],psi,g$K[i],next_order=TRUE)
    c(leading=q1['P_conditional'],combined=q2['P_conditional'])},numeric(2)))
  g[,`:=`(leading=v[,1],combined=v[,2],scan=scan_value)];g
}
cache<-'psi_validation/data/psi05_BI_curves.csv'
if(file.exists(cache))curves<-fread(cache)else{curves<-dense_for('exact');fwrite(curves,cache)}

theme_psi<-theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),
  strip.background=element_rect(fill='grey88'),plot.title=element_text(size=15),
  plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0),legend.position='bottom')
draw_leading<-function(probability,rho_value,theta_value,scan_value){
  a<-copy(curves[scan==scan_value&rho==rho_value&theta==theta_value]);s<-copy(sim[scan==scan_value&rho==rho_value&theta==theta_value])
  if(probability=='unconditional'){a[,value:=leading*(1-1/R0)];s[,`:=`(estimate=P_unconditional,low=uncond_low,high=uncond_high)]
    ttl<-'Unconditional persistence probability';yl<-ttl}else{a[,value:=leading];s[,`:=`(estimate=P_conditional,low=cond_low,high=cond_high)]
    ttl<-'Conditional persistence probability';yl<-'Persistence probability, conditional on not fizzling'}
  p<-ggplot()+geom_line(data=a,aes(R0-1,value),colour='#009E73',linewidth=.9)+
    geom_errorbar(data=s,aes(R0-1,ymin=low,ymax=high),width=0,colour='grey25',linewidth=.42)+
    geom_point(data=s,aes(R0-1,estimate),colour='#2C7FB8',size=2)+facet_wrap(~K,ncol=2,
      labeller=labeller(K=function(x)paste0('K = ',if(scan_value=='tau')format(as.numeric(x),scientific=TRUE)else x)))+
    scale_x_log10(breaks=c(.03,.05,.1,.2,.5,1,2,5),labels=c('.03','.05','.10','.20','.50','1','2','5'))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),expand=expansion(mult=c(.01,.03)))+
    labs(title=ttl,subtitle=sprintf('Boundary-layer-independent leading order; psi = 0.5, rho = %.2f, theta = %g; I0 = 1',rho_value,theta_value),
      x=expression(R[0]-1~'(log scale)'),y=yl,caption='Blue points and grey bars: stochastic estimates and 95% Wilson intervals; green curve: BI theory.')+theme_psi
  if(probability=='unconditional'){ref<-unique(a[,.(K,R0)]);ref[,value:=1-1/R0];p<-p+geom_line(data=ref,aes(R0-1,value,group=K),colour='grey60',linetype='dashed',linewidth=.55)}
  p
}
write_leading<-function(probability,file){cairo_pdf(file,10.5,8.4,onefile=TRUE)
  for(rho in c(.01,.02,.05,.10))for(theta in c(0,.5,1))print(draw_leading(probability,rho,theta,'exact'));dev.off()}
write_leading('unconditional','psi_validation/figures/fig12_psi05_BI_unconditional_stochastic_validation.pdf')
write_leading('conditional','psi_validation/figures/fig13_psi05_BI_conditional_stochastic_validation.pdf')

draw_next<-function(rho_value,theta_value){
  a<-melt(curves[scan=='exact'&rho==rho_value&theta==theta_value],id.vars=c('rho','theta','K','R0'),
    measure.vars=c('leading','combined'),variable.name='method',value.name='value')
  a[,method:=factor(method,c('leading','combined'),c('Leading BI','Combined next-order BI'))]
  s<-sim[scan=='exact'&rho==rho_value&theta==theta_value]
  ggplot()+geom_line(data=a,aes(R0-1,value,colour=method),linewidth=.9)+
    geom_errorbar(data=s,aes(R0-1,ymin=cond_low,ymax=cond_high),width=0,colour='grey25',linewidth=.42)+
    geom_point(data=s,aes(R0-1,P_conditional),colour='#2C7FB8',size=2)+facet_wrap(~K,ncol=2,labeller=labeller(K=function(x)paste0('K = ',x)))+
    scale_x_log10(breaks=c(.05,.1,.2,.5,1,2,5))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),expand=expansion(mult=c(.01,.03)))+
    scale_colour_manual(values=c('Leading BI'='#009E73','Combined next-order BI'='#D55E00'))+
    labs(title='Conditional persistence probability',subtitle=sprintf('Boundary-independent comparison; psi = 0.5, rho = %.2f, theta = %g; I0 = 1',rho_value,theta_value),
      x=expression(R[0]-1~'(log scale)'),y='Persistence probability, conditional on not fizzling',colour=NULL,
      caption='Blue points and grey bars: stochastic estimates and 95% Wilson intervals. Green: leading BI; orange: combined next-order BI.')+theme_psi
}
cairo_pdf('psi_validation/figures/fig15_psi05_BI_next_order_stochastic_curves.pdf',10.5,8.4,onefile=TRUE)
for(rho in c(.01,.02,.05,.10))for(theta in c(0,.5,1))print(draw_next(rho,theta));dev.off()
cat('psi=0.5 Figures 12, 13, and 15 complete\n')
