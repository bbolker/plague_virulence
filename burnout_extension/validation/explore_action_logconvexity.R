library(data.table)
logPhi<-function(t,de){f<-exp(-t);r<-t/(-expm1(-t));lo<-log(r);# integrate in v from log r to t, with stable normalization
 ans<-integrate(function(v)exp(de*(t-v)-v)/(-expm1(-v)),lo,t,rel.tol=1e-10,abs.tol=1e-12,subdivisions=800)$value
 log(ans)
}
z<-CJ(t=exp(seq(log(.001),log(200),length.out=220)),delta=seq(0,1,.025));z[,second:=mapply(function(t,d){h<-t*1e-3;(logPhi(t+h,d)-2*logPhi(t,d)+logPhi(t-h,d))/h^2},t,delta)];fwrite(z,'validation/data/action_logconvexity_exploration.csv');print(z[order(second)][1:8]);cat('negative count',nrow(z[second< -1e-6]),'\n')
