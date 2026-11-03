source('validation/R/theory.R');source('validation/R/extrema.R');source('validation/R/fold_geometry.R');library(data.table)
th<-.55;K<-exp(3.75)
F<-function(par){r<-par[1];la<-par[2];rho<-1/la;cd<-C_log_derivatives(r,th);ad<-action_barrier_derivatives(r,th);q<-r*(r-1);qp<-2*r-1;c0<-cd['d1']+.5/(r-1)-th/(2*r);c1<-cd['d2']-.5/(r-1)^2+th/(2*r^2)
 U<-unname(-q*c0);V<-unname(q*ad['d1']);U1<-unname(-q*(qp*c0+q*c1));V1<-unname(q*(qp*ad['d1']+q*ad['d2']));s<-U+la*V;st<-U1+la*V1
 lb<-log(K)+unname(cd['logC'])+.5*(log(rho)+log(r-1)-th*log(r)-log(2*pi))-unname(ad['value'])*la;B<-exp(lb);hd<-fold_HD(B)
 c(unname(hd['H'])-s,-s*unname(hd['D'])-st)
}
x<-fread('validation/data/narrow_fold_geometry_turns.csv')[theta==th];starts<-list()
for(di in unique(x$direction)){z<-x[direction==di][order(rho)];ix<-which(head(z$logK,-1)-log(K)<0 & tail(z$logK,-1)-log(K)>0 | head(z$logK,-1)-log(K)>0 & tail(z$logK,-1)-log(K)<0)
 for(i in ix){frac<-(log(K)-z$logK[i])/(z$logK[i+1]-z$logK[i]);starts[[length(starts)+1]]<-c(z$R0[i]+frac*(z$R0[i+1]-z$R0[i]),1/(z$rho[i]+frac*(z$rho[i+1]-z$rho[i])))}
}
results<-lapply(starts,function(p){for(it in 1:10){v<-F(p);hs<-c(1e-4,1e-5);J<-sapply(1:2,function(j){pp<-pm<-p;pp[j]<-pp[j]+hs[j];pm[j]<-pm[j]-hs[j];(F(pp)-F(pm))/(2*hs[j])});step<-solve(J,v);p<-p-step;if(max(abs(step))<1e-10)break};v<-F(p);data.table(theta=th,K=K,R0=p[1],lambda=p[2],rho=1/p[2],E=v[1],Et=v[2])})
z<-rbindlist(results)[order(rho)];print(z);fwrite(z,'validation/data/reentrant_fold_centers.csv')
