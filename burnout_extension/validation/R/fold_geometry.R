# Research-only fold enumeration, retaining the full numerical amplitude.
# A fold is P_R=P_RR=0. At fixed R,theta the exact convexity reduction
# gives at most two roots in B, except for the V=V1=0 exceptional case.
fold_HD<-function(B){
 if(B==0)return(c(H=1,D=0,Fs=1))
 if(B<1e-3){n<-0:12;H<-sum(B^n/factorial(n+1));D<-sum(n*B^n/factorial(n+1));Dp<-sum((1:12)^2*B^(0:11)/factorial(2:13))
 }else{eb<-exp(B);H<-expm1(B)/B;D<-eb-H;Dp<-eb-D/B}
 c(H=H,D=D,Fs=D+H*B*Dp/D)
}
fold_H_inverse<-function(s){
 if(s<=1)return(0)
 hi<-max(2,log(s)+log1p(log(s))+2)
 while(fold_HD(hi)['H']<s)hi<-2*hi
 uniroot(function(B)fold_HD(B)['H']-s,c(0,hi),tol=1e-12)$root
}
fold_points<-function(U,V,U1,V1,logC,A,R0,theta){
 if(!all(is.finite(c(U,V,U1,V1,logC,A))))return(NULL)
 if(V==0)return(NULL) # exceptional slice handled separately in the analysis
 lower<-if(V>0)fold_H_inverse(max(1,U))else 0
 if(V<0 && U<=1)return(NULL)
 upper<-if(V<0)fold_H_inverse(U)else Inf
 gamma<-V*U1-U*V1
 g<-function(B){h<-fold_HD(B);unname(V*h['D']+V1+gamma/h['H'])}
 # G_s=V F'(s)+V1 has at most one zero, because F''>0.
 extremum<-NA_real_;target<--V1/V
 if(target>1){hi<-max(2,log(target)+3)
  while(fold_HD(hi)['Fs']<target)hi<-2*hi
  extremum<-uniroot(function(B)fold_HD(B)['Fs']-target,c(0,hi),tol=1e-12)$root
 }
 if(!is.finite(upper)){
  upper<-max(2,lower+2,if(is.finite(extremum))extremum+2 else 0)
  while(sign(g(upper))!=sign(V) && upper<600)upper<-min(600,upper*1.5)
 }
 cuts<-sort(unique(c(lower,upper,if(is.finite(extremum)&&extremum>lower&&extremum<upper)extremum)))
 roots<-numeric()
 for(j in seq_len(length(cuts)-1)){
  aa<-cuts[j];bb<-cuts[j+1];ga<-g(aa);gb<-g(bb)
  if(ga*gb<0)roots<-c(roots,uniroot(g,c(aa,bb),tol=1e-12)$root)
 }
 if(!length(roots))return(NULL)
 out<-lapply(roots,function(B){s<-unname(fold_HD(B)['H']);lambda<-(s-U)/V
  if(lambda<=0||B<=0)return(NULL)
  rho<-1/lambda;logb<-logC+.5*(log(rho)+log(R0-1)-theta*log(R0)-log(2*pi))-A/rho
  data.frame(B=B,s=s,rho=rho,logK=log(B)-logb,lambda=lambda,fold_residual=g(B))
 })
 do.call(rbind,out)
}
