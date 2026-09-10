library(data.table);library(ggplot2)
sim<-fread('psi_validation/data/psi05_stochastic_results.csv');curves<-fread('psi_validation/data/psi05_BI_curves.csv')
rho_show<-c(.01,.05,.10);sim<-sim[rho%in%rho_show];curves<-curves[rho%in%rho_show]
curves[,value:=pmin(.9999,pmax(.0001,leading))];sim[,`:=`(estimate=pmin(.9999,pmax(.0001,P_conditional)),low=pmin(.9999,pmax(.0001,cond_low)),high=pmin(.9999,pmax(.0001,cond_high)))]
cols<-c('0.01'='#0072B2','0.05'='#D55E00','0.1'='#009E73')
p<-ggplot(curves,aes(R0-1,value,colour=factor(rho),linetype=factor(theta),group=interaction(factor(rho),factor(theta),K)))+
 geom_line(linewidth=.85,na.rm=TRUE,alpha=.85)+geom_pointrange(data=sim,aes(R0-1,y=estimate,ymin=low,ymax=high,colour=factor(rho),shape=factor(theta)),inherit.aes=FALSE,linewidth=.38,size=.65,fatten=1.8)+
 facet_wrap(~K,ncol=2,labeller=labeller(K=function(x)paste0('K = ',x)))+
 scale_colour_manual(values=cols,name=expression(rho))+scale_linetype_manual(values=c('0'='solid','0.5'='22','1'='44'),name=expression(theta))+scale_shape_manual(values=c('0'=16,'0.5'=17,'1'=15),name=expression(theta))+
 scale_x_log10(breaks=c(.05,.1,.2,.5,1,2,5),labels=c('.05','.10','.20','.50','1','2','5'))+scale_y_continuous(trans='logit',limits=c(.0001,.9999),breaks=c(.001,.01,.05,.1,.25,.5,.75,.9,.99,.999),labels=function(x)format(x,trim=TRUE,scientific=FALSE),expand=expansion(mult=c(.01,.03)))+
 labs(title='Conditional persistence probability: condensed BI validation',subtitle='psi = 0.5; reused stochastic estimates; local boundary-layer-independent leading approximation',x=expression(R[0]-1~'(log scale)'),y='Conditional persistence probability',caption='Colour identifies rho. Solid/22/44 curves identify theta = 0/0.5/1; circle/triangle/square points identify the same theta values.')+
 theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),strip.background=element_rect(fill='grey88'),legend.position='bottom',legend.box='vertical',legend.box.just='left',plot.title=element_text(size=15),plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0))+
 guides(colour=guide_legend(order=1),linetype=guide_legend(order=2,override.aes=list(linewidth=1.2)),shape=guide_legend(order=2))
dir.create('psi_validation/figures',FALSE,TRUE);cairo_pdf('psi_validation/figures/fig13_psi05_BI_comb.pdf',width=12,height=9,onefile=FALSE);print(p);dev.off();cat('old BI comb complete\n')
