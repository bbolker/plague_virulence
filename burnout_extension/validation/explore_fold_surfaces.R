source('validation/R/theory.R');source('validation/R/extrema.R');source('validation/R/large_shape_jets.R');source('validation/R/fold_geometry.R');library(data.table)
ths<-c(seq(0,.95,.025),.975,.99,.995,.999,.9999,.99999,.999999,1)
all<-list();jj<-0;constituents<-list();cc<-0
for(th in ths){
 de<-1-th
 rr<-unique(c(seq(1.1,20,length.out=320),exp(seq(log(20),log(600),length.out=100)),if(de>0)exp(seq(log(600),log(max(1000,80/de)),length.out=180))else numeric()))
 for(r in rr){
  if(r<=600){cd<-C_log_derivatives(r,th);ad<-action_barrier_derivatives(r,th);C1<-cd['d1'];C2<-cd['d2'];logC<-cd['logC'];A<-ad['value'];A1<-ad['d1'];A2<-ad['d2']
  }else{cd<-large_shape_derivatives4(r,th);C1<-cd['C1'];C2<-cd['C2'];logC<-cd['logC'];A<-large_A(r,th);A1<-large_action_d1(r,th);A2<-cd['A2']}
  q<-r*(r-1);qp<-2*r-1;c0<-C1+.5/(r-1)-th/(2*r);c1<-C2-.5/(r-1)^2+th/(2*r^2)
  U<-unname(-q*c0);V<-unname(q*A1);U1<-unname(-q*(qp*c0+q*c1));V1<-unname(q*(qp*A1+q*A2))
  cc<-cc+1;constituents[[cc]]<-data.table(R0=r,theta=th,U=U,V=V,U1=U1,V1=V1,logC=unname(logC),A=unname(A))
  pts<-fold_points(U,V,U1,V1,unname(logC),unname(A),r,th)
  if(!is.null(pts)){jj<-jj+1;all[[jj]]<-cbind(data.table(R0=r,theta=th),as.data.table(pts))}
 }
 cat('theta',th,'fold records',if(length(all))nrow(rbindlist(all))else 0,'\n');flush.console()
}
saveRDS(rbindlist(constituents),'validation/data/fold_constituents.rds')
x<-rbindlist(all);fwrite(x,'validation/data/fold_surface_scan.csv');print(x[,.(n=.N,min_logK=min(logK),max_logK=max(logK),minrho=min(rho),maxrho=max(rho)),by=theta]);cat('multiple fold roots per R:',max(x[,.N,by=.(R0,theta)]$N),'\n')
