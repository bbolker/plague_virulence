# Original fixed Fig. 1 grid and all its stochastic points; no simulation reruns.
source('validation/R/theory.R')
source('validation/R/kendall.R')
source('validation/R/ode_reference.R')
source('validation/R/validation_plot_data.R')
library(data.table)
library(ggplot2)
sim<-fread('validation/data/direct_entry_comparison_original_grid.csv',integer64='double')
stopifnot(nrow(sim)==280L)
C_explicit<-function(R0,theta)exp(unname(C_log_derivatives(R0,theta)['logC']))
rr<-sort(unique(c(1+exp(seq(log(.03),log(5),length.out=80)),sim$R0)))
keys<-CJ(R0=rr,theta=seq(0,1,.25))
cache<-'validation/data/direct_entry_comparison_dense_curves.csv'
if ('--plot-only' %in% commandArgs(trailingOnly=TRUE)) {
curves<-fread(cache)
} else {
rows<-list()
for(i in seq_len(nrow(keys))){
 r<-keys$R0[i];th<-keys$theta[i];rho<-.01
 const<-tryCatch(list(D=D_regularized(r,th),C=C_explicit(r,th),
                      ode=ode_reference(r,rho,th,1e6)),error=function(e)list(error=conditionMessage(e)))
 for(K in c(1e6,1e7,1e8,1e9)){
  z<-data.table(R0=r,theta=th,rho=rho,K=K,P_BI=stable_validation_probability(r,rho,th,K),
       P_first_K=NA_real_,P_second_K=NA_real_,entry_available=FALSE,status='OK')
  if(!is.null(const$error))z[,status:=const$error] else {
   xf<-x_final(r);delta<-1-r*xf;hf<-h_theta(xf,th);a<-hf/delta;b<-r*xf/delta
   cc<-delta*h_prime(xf,th)+r*hf;kappa<-a*cc/(2*delta^2)
   yb<-sqrt(rho*h_theta(1/r,th)/K);L<-log(const$C/yb)
   x1<-xf+rho*a*L+b*yb
   x2<-x1+rho^2*(kappa*L^2+const$D)+rho*yb/delta^2*(r*a*(L+1)-b*cc)+b*r*yb^2/(2*delta^2)
   available<-is.finite(const$ode$y_min)&&const$ode$y_min<yb
   z[,entry_available:=available]
   if(available)for(nm in c('first','second')){
    x<-if(nm=='first')x1 else x2
    if(is.finite(x)&&x>xf&&x<1/r){
     q<-tryCatch(kendall_quantities(x,r,rho,th,K,yline=yb)$P1,error=function(e)NA_real_)
     z[,(paste0('P_',nm,'_K')):=q]
    }
   }
  }
  rows[[length(rows)+1L]]<-z
 }
 if(i%%40==0 || i==nrow(keys)){
  fwrite(rbindlist(rows,fill=TRUE),cache)
  cat('curve keys',i,'/',nrow(keys),'\n');flush.console()
 }
}
curves<-rbindlist(rows)
}
stopifnot(all(curves$status=='OK'))
# NA gaps remain gaps; do not connect across inadmissible entry intervals.
long<-melt(curves,id.vars=c('R0','theta','K'),measure.vars=c('P_BI','P_first_K','P_second_K'),
 variable.name='method',value.name='prediction')
long[,method:=factor(method,levels=c('P_BI','P_first_K','P_second_K'),
 labels=c('BI','First-order entry + Kendall','Second-order entry + Kendall'))]
setorder(long,K,theta,method,R0)
long[,segment:=cumsum(is.na(prediction)),by=.(K,theta,method)]
long[,theta_label:=factor(theta,levels=seq(0,1,.25))]
sim[,theta_label:=factor(theta,levels=seq(0,1,.25))]
clip<-function(p)pmax(.001,pmin(.999,p))
long[,yp:=qlogis(clip(prediction))]
sim[,`:=`(yp=qlogis(clip(P_conditional)),lo=qlogis(clip(cond_low)),hi=qlogis(clip(cond_high)))]
breaks<-c(.001,.01,.05,.25,.5,.75,.95,.99,.999)
pal<-c('#0072B2','#56B4E9','#009E73','#E69F00','#D55E00')
if(capabilities('aqua')) quartzFonts(STIX=quartzFont(c('STIXGeneral-Regular','STIXGeneral-Bold','STIXGeneral-Italic','STIXGeneral-BoldItalic')))
base<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),
 panel.grid.major=element_line(colour='grey90',linewidth=.25),
 strip.background=element_rect(fill='grey94'),legend.position='bottom',
 plot.caption=element_text(hjust=0,size=8.5),legend.box='vertical',axis.title.x=element_text(family=if(capabilities('aqua'))'STIX' else 'serif'))
axes<-list(scale_x_log10(breaks=c(.03,.1,.3,1,3,5),labels=c('0.03','0.1','0.3','1','3','5')),
 scale_y_continuous(breaks=qlogis(breaks),labels=format(breaks,trim=TRUE),limits=qlogis(c(.001,.999))),
 labs(x=bquote(.(intToUtf8(0x211B))[0]-1~'(log scale)'),y='Conditional persistence probability'))
p<-ggplot()+geom_line(data=long[is.finite(yp)],aes(R0-1,yp,colour=theta_label,linetype=method,
 group=interaction(theta,method,segment)),linewidth=.72)+
 geom_errorbar(data=sim,aes(R0-1,ymin=lo,ymax=hi,colour=theta_label),width=0,linewidth=.3,alpha=.65)+
 geom_point(data=sim,aes(R0-1,yp,colour=theta_label),shape=16,size=1.65)+
 facet_wrap(~K,ncol=2,labeller=as_labeller(function(x)paste0('italic(K)==10^',round(log10(as.numeric(x)))),default=label_parsed))+
 scale_colour_manual(values=pal,name=expression(theta))+
 scale_linetype_manual(values=c('BI'='solid','First-order entry + Kendall'='dashed','Second-order entry + Kendall'='dotdash'),name=NULL)+
 axes+base+labs(caption='All 280 original simulation points are retained (rho = 0.01). Curves stop where the ODE has no threshold crossing or the entry estimate is inadmissible.\nPoints: simulations; bars: 95% Wilson intervals. All three predictions are conditional on establishment.')
# A second layout separates theta to make the method differences easy to read.
method_pal<-c('BI'='#0072B2','First-order entry + Kendall'='#888888','Second-order entry + Kendall'='#D55E00')
for(d in c('long','sim')){
 z<-get(d);z[,theta_facet:=factor(sprintf('theta==%g',theta),levels=sprintf('theta==%g',seq(0,1,.25)))];
 z[,K_facet:=factor(sprintf('K==10^%d',round(log10(K))),levels=sprintf('K==10^%d',6:9))];assign(d,z)
}
pdetail<-ggplot()+geom_line(data=long[is.finite(yp)],aes(R0-1,yp,colour=method,linetype=method,
 group=interaction(method,segment)),linewidth=.65)+
 geom_errorbar(data=sim,aes(R0-1,ymin=lo,ymax=hi),width=0,linewidth=.25,alpha=.5)+
 geom_point(data=sim,aes(R0-1,yp),colour='black',shape=16,size=1.1)+
 facet_grid(theta_facet~K_facet,labeller=label_parsed)+
 scale_colour_manual(values=method_pal,name=NULL)+
 scale_linetype_manual(values=c('BI'='solid','First-order entry + Kendall'='dashed','Second-order entry + Kendall'='dotdash'),name=NULL)+axes+base+
 labs(caption='Same full parameter grid as the four-panel layout. Black points: simulations; bars: 95% Wilson intervals.\nOne theta per row; one carrying capacity per column. No selection by simulation probability or prediction error.')
out<-'validation/figures/previews';dir.create(out,recursive=TRUE,showWarnings=FALSE)
save_plot<-function(p,name,w,h){
 if(capabilities('aqua')){
  quartz(file=file.path(out,paste0(name,'.pdf')),type='pdf',width=w,height=h);print(p);dev.off()
  png(file.path(out,paste0(name,'.png')),type='quartz',width=w*180,height=h*180,res=180);print(p);dev.off()
 }else{ggsave(file.path(out,paste0(name,'.pdf')),p,width=w,height=h);ggsave(file.path(out,paste0(name,'.png')),p,width=w,height=h,dpi=180)}
}
save_plot(p,'entry_curve_comparison_fig1_layout',9,7.7)
save_plot(pdetail,'entry_curve_comparison_by_theta',11.5,10)
cat('DONE: curve comparison previews\n')
if ('--paper-export' %in% commandArgs(trailingOnly=TRUE)) {
 paper_plot<-pdetail+labs(caption=NULL)+
  scale_x_log10(breaks=c(.03,.1,.3,1,3),labels=c('0.03','0.1','0.3','1','3'))+
  scale_y_continuous(breaks=qlogis(c(.001,.01,.05,.5,.95,.99,.999)),
   labels=c('0.001','0.010','0.050','0.500','0.950','0.990','0.999'),
   limits=qlogis(c(.001,.999)))+
  theme(axis.text=element_text(size=8.5),strip.text=element_text(size=9),legend.text=element_text(size=9))
 paper_file<-'validation/figures/paper/entry_approximation_comparison.pdf'
 if(capabilities('aqua')){
  quartz(file=paper_file,type='pdf',width=7.15,height=8);print(paper_plot);dev.off()
 }else ggsave(paper_file,paper_plot,width=7.15,height=8,device=cairo_pdf)
}
