source('validation/R/theory.R');source('validation/R/extrema.R');library(data.table)
st<-function(t,th,rv){r<-1/(-expm1(t));cd<-C_log_derivatives(r,th);ad<-action_barrier_derivatives(r,th);q<-r*(r-1);qp<-2*r-1;c0<-cd['d1']+1/(2*(r-1))-th/(2*r);c1<-cd['d2']-1/(2*(r-1)^2)+th/(2*r^2);q*(-qp*c0-q*c1+(qp*ad['d1']+q*ad['d2'])/rv)}
z<-fread('validation/data/critical_geometry_UV.csv');z<-z[R0>1.01&R0<100];out<-list();j<-0
for(rv in c(.001,.003,.01,.03,.1,.3,1,3,10)){
 d<-z[U+V/rv>1 & U1+V1/rv<0];if(!nrow(d))next;d<-d[seq(1,.N,by=4)]
 for(i in 1:nrow(d)){t<-log1p(-1/d$R0[i]);th<-d$theta[i];h<-min(abs(t)*1e-3,1e-4);v<-c(st(t-h,th,rv),st(t,th,rv),st(t+h,th,rv));if(any(v>=0))next
 lc<-(log(-v[1])-2*log(-v[2])+log(-v[3]))/h^2;j<-j+1;out[[j]]<-data.table(R0=d$R0[i],theta=th,rho=rv,logcurvature=as.numeric(lc))
 };cat('rho',rv,'done\n')
}
a<-rbindlist(out);fwrite(a,'validation/data/negative_slope_logconcavity_exploration.csv');print(a[order(-logcurvature)][1:10]);cat('positive',nrow(a[logcurvature>1e-3]),'of',nrow(a),'\n')
