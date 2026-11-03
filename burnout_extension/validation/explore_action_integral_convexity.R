source('validation/R/theory.R');source('validation/R/extrema.R');source('validation/R/shape_jets.R');library(data.table)
rr<-c(seq(1.03,6,length.out=100),exp(seq(log(6),log(300),length.out=100)))
out<-list();j<-0
for(th in seq(0,1,.05))for(r in rr){
 f<-x_final(r)
 a<-action_barrier_derivatives(r,th);aa<-shape_derivatives4(r,th)
 de<-1-th; J<-1-a['d1']/f^de
 if(J<1e-5)next
 fp<-f*(f-1)/(1-r*f); z1<-fp/f
 z2<-fp/(1-r*f)+(f-1)*(f+r*fp)/(1-r*f)^2
 JR<--f^(-de)*(aa['A2']-de*z1*a['d1'])
 JRR<--f^(-de)*(aa['A3']-2*de*z1*aa['A2']+(de^2*z1^2-de*z2)*a['d1'])
 # S=-log f, R_S= -1/z1; R_SS= -z2/z1^3
 logJSS<-(JRR/J-(JR/J)^2)/z1^2-JR/J*z2/z1^3
 j<-j+1;out[[j]]<-data.table(R0=r,theta=th,J=as.numeric(J),logJSS=as.numeric(logJSS))
}
x<-rbindlist(out);print(x[order(logJSS)][1:12]);fwrite(x,'validation/data/action_integral_logconvexity.csv')
