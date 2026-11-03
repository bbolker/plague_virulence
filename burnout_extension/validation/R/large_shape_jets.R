source('validation/R/shape_jets.R')
large_shape_derivatives4<-function(r,theta,nquad=192){
 R<-jvar(r);one<-jconst(1);de<-1-theta
 if(de*r>80){lc<-jsub(jconst(-.5772156649015329-log(de)),jlog(R))}else{
  quad<-statmod::gauss.quad(nquad,kind='legendre');v<-(quad$nodes+1)/2;w<-quad$weights/2
  A<-jscale(R,de);aa<-jdiv(jexpm1(jscale(A,-v)),jconst(v));K<-lapply(aa,function(q)sum(w*q))
  tt<-55*v;ww<-55*w;u<-exp(-tt);om<--expm1(-tt);S<-jsub(R,jconst(tt));num<-jscale(jsub(one,jscale(R,u)),om);den<-jsub(jscale(R,om),jconst(tt));co<-jmul(jexp(jscale(S,-de)),jsub(jdiv(num,den),jinv(S)));co<-lapply(co,function(q)sum(ww*q));lc<-jadd(K,co)
 }
 f<-jexp(jscale(R,-1));A2<-jsub(jdiv(jpow(R,theta-2),jsub(R,one)),jdiv(jmul(jexp(jscale(R,-de)),jadd(jconst(de),jscale(f,theta))),jsub(one,jmul(R,f))))
 c(logC=lc[[1]],C1=lc[[2]],C2=2*lc[[3]],C3=6*lc[[4]],C4=24*lc[[5]],A2=A2[[1]],A3=A2[[2]],A4=2*A2[[3]])
}
large_action_d1<-function(r,th){de<-1-th;n<-0:18;exp(-de*r)-r^(-de-1)*sum(r^-n/(n+de+1))}

large_A<-function(r,th){d<-1-th;n<-1:18
 if(d==0)return(r-log(r)-1+sum(r^-n/(n*(n+1))))
 t<-d*r
 lead<-if(t<40)exp(-t)*expm1(d*(r-log(r))-log1p(d))/d else exp(-d*log(r))/(d*(1+d))
 lead+exp(-d*log(r))*sum(r^-n/((n+d)*(n+d+1)))
}
