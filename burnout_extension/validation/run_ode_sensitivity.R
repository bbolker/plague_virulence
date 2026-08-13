source('validation/R/ode_reference.R');library(data.table)
pts<-data.frame(R0=c(1.2,3,5),rho=c(.005,.01,.05),theta=c(0,.5,1),K=1e4)
out<-list();j<-0
for(i in 1:nrow(pts)) for(y0 in c(1e-8,1e-10,1e-12)) {
 j<-j+1;o<-do.call(ode_reference,c(as.list(pts[i,]),list(y0=y0,init='manifold',dt=.02)))
 out[[j]]<-data.frame(pts[i,],init='manifold',y0=y0,status=o$status,x_in=o$x_in)
}
for(i in 1:nrow(pts)) {
 j<-j+1;o<-do.call(ode_reference,c(as.list(pts[i,]),list(init='finiteK',dt=.02)))
 out[[j]]<-data.frame(pts[i,],init='finiteK',y0=1/pts$K[i],status=o$status,x_in=o$x_in)
}
fwrite(rbindlist(out),'validation/data/ode_sensitivity.csv')
