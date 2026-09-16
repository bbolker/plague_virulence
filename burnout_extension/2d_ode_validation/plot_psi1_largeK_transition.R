library(data.table)
library(ggplot2)

s<-fread('2d_ode_validation/data/psi1_largeK_transition_pilot.csv')
t<-fread('2d_ode_validation/data/psi1_largeK_transition_2d_ode.csv')
d<-merge(s,t[,.(rho,K,q,P_kendall,status)],by=c('rho','K','q'),all.x=TRUE)
z<-qnorm(.975);n<-d$established;x<-d$persistent;den<-1+z^2/n
ctr<-(x/n+z^2/(2*n))/den
hw<-z*sqrt((x/n)*(1-x/n)/n+z^2/(4*n^2))/den
d[,`:=`(cond_low=ctr-hw,cond_high=ctr+hw,
  theory=ifelse(status=='OK',P_kendall,NA_real_),
  theory_status=ifelse(status=='OK','available','unresolved'))]
fwrite(d,'2d_ode_validation/data/psi1_largeK_transition_comparison.csv')

rho_lab<-function(x)sprintf('rho = %.2f',as.numeric(x))
k_lab<-function(x)paste0('K = ',format(as.numeric(x),scientific=TRUE))
p<-ggplot(d,aes(q,P_conditional,group=interaction(rho,K)))+
  annotate('rect',xmin=-Inf,xmax=Inf,ymin=0,ymax=.05,fill='grey93')+
  annotate('rect',xmin=-Inf,xmax=Inf,ymin=.95,ymax=1,fill='grey93')+
  geom_hline(yintercept=c(.05,.95),colour='grey75',linetype='dotted')+
  geom_line(colour='#2C7FB8',linewidth=.65)+
  geom_errorbar(aes(ymin=cond_low,ymax=cond_high),colour='#2C7FB8',width=.035,linewidth=.42)+
  geom_point(colour='#2C7FB8',size=2)+
  geom_point(aes(y=theory),colour='#009E73',shape=17,size=2.3,na.rm=TRUE)+
  facet_grid(K~rho,labeller=labeller(K=k_lab,rho=rho_lab))+
  scale_x_continuous(breaks=c(2.5,3,3.5,4,4.5))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),expand=expansion(mult=c(.01,.02)))+
  labs(title='Transition-focused conditional persistence: psi = 1, theta = 0.5',
    subtitle='The moving grid follows R0 = rho*sqrt(K)/q; grey bands mark near-0 and near-1 outcomes',
    x=expression(q==rho*sqrt(K)/R[0]),y='Persistence probability, conditional on establishment',
    caption=paste0('Blue points and 95% Wilson intervals: adaptive tau-leaping pilot (epsilon = 0.01; ',
      '50-200 paths per point). Green triangles: numerical finite-horizon Kendall when the 2D ODE event solve is resolved.\n',
      'Missing green triangles indicate BAD_CROSSING_DIRECTION or NO_RECOVERY_PEAK, not zero probability.'))+
  theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),strip.background=element_rect(fill='grey88'),
    plot.title=element_text(size=15),plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0),
    panel.spacing=unit(4,'pt'))

dir.create('2d_ode_validation/figures',FALSE,TRUE)
cairo_pdf('2d_ode_validation/figures/fig16_psi1_largeK_transition_focus.pdf',12,10,onefile=FALSE)
print(p);dev.off()
