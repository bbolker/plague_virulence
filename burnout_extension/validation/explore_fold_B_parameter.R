# Inverse parameterization avoids folds with exponentially small B being
# invisible on an R grid. The B=0 slice describes the limiting K=0 fold.
source('validation/R/fold_geometry.R');library(data.table)
x<-readRDS('validation/data/fold_constituents.rds');bs<-unique(c(0,exp(seq(-30,log(25),length.out=280))));out<-list();kk<-0
for(th in unique(x$theta)){
 a<-x[theta==th][order(R0)];a<-a[!duplicated(log(R0))];lr<-log(a$R0)
 for(B in bs){h<-fold_HD(B);gamma<-a$V*a$U1-a$U*a$V1;g<-a$V*h['D']+a$V1+gamma/h['H'];ii<-which(head(g,-1)*tail(g,-1)<0)
  for(i in ii){y<-lr[i]-g[i]*(lr[i+1]-lr[i])/(g[i+1]-g[i]);point<-function(name)approx(lr,a[[name]],xout=y)$y
   U<-point('U');V<-point('V');la<-(h['H']-U)/V;if(!is.finite(la)||la<=0)next
   r<-exp(y);rho<-1/la;lc<-point('logC');A<-point('A');logb<-lc+.5*(log(rho)+log(r-1)-th*y-log(2*pi))-A/rho
   kk<-kk+1;out[[kk]]<-data.table(theta=th,B=B,R0=r,rho=rho,logK=if(B>0)log(B)-logb else -Inf,logb=logb,lower_R=a$R0[i],upper_R=a$R0[i+1])
  }
 }
}
z<-rbindlist(out);fwrite(z,'validation/data/fold_B_parameter_scan.csv');print(z[,.N,by=.(theta,B)][,.(max_roots=max(N)),by=theta]);print(z[B==0]);cat('multiple R roots at fixed theta,B:',nrow(z[,.N,by=.(theta,B)][N>1]),'\n')
