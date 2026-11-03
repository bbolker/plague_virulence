source('validation/R/theory.R');source('validation/R/extrema.R');library(data.table);library(ggplot2)
m<-fread('validation/data/global_shape_audit_regions.csv');rt<-fread('validation/data/global_shape_audit_roots.csv');rt[,R0:=as.numeric(R0)];base<-readRDS('validation/data/global_shape_audit_base_stable.rds')
out<-'validation/figures/previews';theme0<-theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),legend.position='bottom')
quartzFonts(STIX=quartzFont(c('STIXGeneral-Regular','STIXGeneral-Bold','STIXGeneral-Italic','STIXGeneral-BoldItalic')))
saveplot<-function(p,name,w,h){quartz(file=file.path(out,name),type='pdf',width=w,height=h);print(p);dev.off()}
m[,label:=factor(ifelse(nroots==1,'Increase - decrease','Increase - decrease - increase - decrease'))]
m[,Klab:=paste0('K = 10^',round(log10(K)))]
p<-ggplot(m,aes(rho,theta,fill=label))+geom_tile()+facet_wrap(~Klab,nrow=1)+scale_fill_manual(values=c('Increase - decrease'='#4477AA','Increase - decrease - increase - decrease'='#EEAA33'),name='BI derivative sign sequence')+scale_x_log10(breaks=c(.002,.005,.01,.02,.04),labels=c('.002','.005','.01','.02','.04'))+labs(x=expression(rho),y=expression(theta))+theme0
saveplot(p,'global_BI_structure_audit.pdf',9,4.5)
configs<-data.table(theta=c(0,0,.5,1),wanted=c(.01,.02,.01,.01))
curves<-list();marks<-list();summary<-list()
for(i in 1:nrow(configs)){
 th<-configs$theta[i];rv<-unique(m$rho)[which.min(abs(unique(m$rho)-configs$wanted[i]))];k<-1e6
 roots<-rt[K==k & theta==th & rho==rv];z<-base[theta==th]
 ymax<-if(th==1)log(600) else max(log(100),roots$logR)+log(1e10)
 y<-sort(unique(c(log(1+exp(seq(log(1e-4),log(20),length.out=800))),seq(log(20),ymax,length.out=1600),roots$logR)))
 r<-exp(y);sel<-r<=600;lc<-A<-numeric(length(r));lc[sel]<-approx(log(z$R0),z$logC,xout=y[sel],rule=2)$y;A[sel]<-approx(log(z$R0),z$A,xout=y[sel],rule=2)$y
 if(any(!sel)){d<-1-th;stopifnot(d>0);yy<-y[!sel];n<-0:12;lc[!sel]<--.5772156649015329-log(d)-yy;A[!sel]<-exp(-d*yy)*vapply(yy,function(qq)sum(exp(-n*qq)/((n+d)*(n+d+1))),0.)}
 lb<-log(k)+lc+.5*(log(rv)+log(r-1)-th*log(r)-log(2*pi))-A/rv
 B<-exp(pmin(lb,700));lpcond<-ifelse(lb< -20,lb,log(-expm1(-B)))
 lp<-log1p(-1/r)+lpcond;logitp<-lp-log(-expm1(lp))
 # Failure flags use >=1, not a claim of accuracy where all flags are below 1.
 hs<-(1-1/r)*r^-th;ys<-rv*hs;ybl<-sqrt(ys/k);L<-lc-log(ybl)
 flags<-pmax(rv*(r-1)*r^-th,rv*r^(1-th)/(r-1),rv*L^2,1/(k*ys))
 name<-sprintf('%s: rho=%.4g, theta=%g, K=10^6',LETTERS[i],rv,th)
 dd<-data.table(R0=r,logitP=logitp,P=exp(lp),flag=flags>=1,case=name)
 dd[,segment:=rleid(flag)]
 full<-copy(dd);full[,view:='Full curve (logit probability)'];full[,vertical:=pmax(-40,logitP)]
 local<-copy(dd[R0<=20]);local[,view:='Detail: R0 <= 20 (probability)'];local[,vertical:=P]
 curves[[i]]<-rbind(full,local)
 roots[,`:=`(case=name,rho=rv)]
 for(j in 1:nrow(roots)){
  idx<-which.min(abs(y-roots$logR[j]));roots$P[j]<-dd$P[idx];roots$logitP[j]<-dd$logitP[idx];roots$flag[j]<-dd$flag[idx]
 }
 mf<-copy(roots);mf[,`:=`(view='Full curve (logit probability)',vertical=pmax(-40,logitP))];ml<-copy(roots[R0<=20]);ml[,`:=`(view='Detail: R0 <= 20 (probability)',vertical=P)]
 marks[[i]]<-rbindlist(list(mf,ml),use.names=TRUE,fill=TRUE);summary[[i]]<-roots
}
d<-rbindlist(curves);mks<-rbindlist(marks,use.names=TRUE,fill=TRUE);mks[,R0:=as.numeric(R0)];fwrite(rbindlist(summary),'validation/data/global_shape_representative_extrema.csv')
d[,view:=factor(view,levels=c('Detail: R0 <= 20 (probability)','Full curve (logit probability)'))];mks[,view:=factor(view,levels=levels(d$view))]
p<-ggplot(d,aes(log10(R0-1),vertical))+geom_line(aes(group=interaction(case,view,segment),alpha=flag),colour='#0072B2',linewidth=.7)+geom_point(data=mks,aes(shape=type),size=2,colour='#D55E00')+facet_wrap(case~view,ncol=2,scales='free')+scale_alpha_manual(values=c('FALSE'=1,'TRUE'=.23),labels=c('No flag triggered','Approximation warning'),name='Approximation warning')+scale_shape_manual(values=c(Maximum=17,Minimum=16),name='Stationary points')+labs(x=expression(log[10](R[0]-1)),y='Probability (detail); logit probability (full, lower limit -40)')+theme0
saveplot(p,'global_BI_representative_curves.pdf',11,10)
