source('validation/R/theory.R')
source('validation/R/extrema.R')
library(data.table); library(ggplot2)
out <- 'validation/figures/previews'
rgrid <- unique(c(1+exp(seq(log(.0001),log(9),length.out=550)),exp(seq(log(10),log(600),length.out=240))))
ths <- c(seq(0,.95,by=.025),.975,1)
rhos <- exp(seq(log(.002),log(.04),length.out=55))
ks <- c(1e6,1e9)
basefile <- 'validation/data/global_shape_audit_base_stable.rds'
if(file.exists(basefile)) base <- readRDS(basefile) else {
 base <- rbindlist(lapply(ths,function(th){
  cat('constituents theta',th,'\n')
  ans<-rbindlist(lapply(rgrid,function(r){
   cd<-C_log_derivatives(r,th);ad<-action_barrier_derivatives(r,th)
   data.table(R0=r,theta=th,logC=unname(cd['logC']),dC=unname(cd['d1']),A=unname(ad['value']),dA=unname(ad['d1']))
  }));ans
 }));saveRDS(base,basefile)
}
# Scaled derivative R0*G/B, computed without overflow or probability rounding.
signfun <- function(r,lc,dc,A,dA,rho,th,K){
 lb<-log(K)+lc+.5*(log(rho)+log(r-1)-th*log(r)-log(2*pi))-A/rho
 L<-dc+1/(2*(r-1))-th/(2*r)-dA/rho
 B<-exp(pmin(lb,700));lex<-ifelse(lb< -10, log1p(B/2+B^2/6),ifelse(B<700,log(expm1(pmin(B,700)))-lb,B-lb))
 r*L+exp(pmin(700,lex-log(r-1)))
}
tailq <- function(y,th){
 d<-1-th;r<-exp(y);n<-0:12
 list(r=r,lc=-.5772156649015329-log(d)-y,dc=-1/r,
 A=exp(-d*y)*sum(exp(-n*y)/((n+d)*(n+d+1))),
 dA=-exp(-(1+d)*y)*sum(exp(-n*y)/(n+1+d)))
}
# Check the asymptotic join, not an assumed root count.
joincheck<-rbindlist(lapply(ths[ths<1],function(th){
 z<-base[theta==th][.N];q<-tailq(log(600),th)
 data.table(theta=th,logC_difference=q$lc-z$logC,A_relative=q$A/z$A-1,dC_scaled_difference=600*(q$dc-z$dC))
}));fwrite(joincheck,'validation/data/global_shape_audit_join.csv')
roots<-list();maps<-list();j<-0
for(k in ks)for(th in ths){
 q<-base[theta==th]
 ys<-if(th<1)seq(log(600),max(50,12/(1-th)),length.out=1800)[-1] else numeric()
 tq<-if(length(ys))tailq(ys,th) else NULL
 # Vectorized tail power series.
 if(length(ys)) {d<-1-th;n<-0:12;tq$A<-exp(-d*ys)*vapply(ys,function(y)sum(exp(-n*y)/((n+d)*(n+d+1))),0.);tq$dA<--exp(-(1+d)*ys)*vapply(ys,function(y)sum(exp(-n*y)/(n+1+d)),0.)}
 for(rv in rhos){
  v<-signfun(q$R0,q$logC,q$dC,q$A,q$dA,rv,th,k);xx<-log(q$R0)
  if(length(ys)){v<-c(v,signfun(tq$r,tq$lc,tq$dc,tq$A,tq$dA,rv,th,k));xx<-c(xx,ys)}
  ii<-which(sign(v[-length(v)])!=sign(v[-1]))
  ff<-function(y){r<-exp(y);if(r<=600){cd<-C_log_derivatives(r,th);ad<-action_barrier_derivatives(r,th);signfun(r,cd['logC'],cd['d1'],ad['value'],ad['d1'],rv,th,k)}else{t<-tailq(y,th);signfun(r,t$lc,t$dc,t$A,t$dA,rv,th,k)}}
  rr<-if(length(ii))vapply(ii,function(i)uniroot(ff,xx[c(i,i+1)],tol=1e-9)$root,0.) else numeric()
  types<-if(length(ii))ifelse(v[ii]>0,'Maximum','Minimum') else character()
  j<-j+1;maps[[j]]<-data.table(rho=rv,theta=th,K=k,nroots=length(rr),pattern=paste(c(if(v[1]>0)'+' else '-',if(length(ii))ifelse(v[ii+1]>0,'+','-')),collapse=' '))
  if(length(rr))roots[[j]]<-data.table(rho=rv,theta=th,K=k,logR=rr,R0=exp(rr),type=types,tail=rr>log(600))
 }
 cat('classified',k,th,'\n')
}
m<-rbindlist(maps);rt<-rbindlist(roots);fwrite(m,'validation/data/global_shape_audit_regions.csv');fwrite(rt,'validation/data/global_shape_audit_roots.csv')
print(m[,.N,by=.(K,nroots,pattern)])
