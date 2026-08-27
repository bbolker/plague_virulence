source('validation/R/theory.R')

ode_reference <- function(R0,rho,theta,K,y0=1e-10,init=c('manifold','finiteK'),
                          dt=.05,tmax=NULL,trajectory=FALSE,yline=NULL) {
  init<-match.arg(init); xs<-1/R0; ystar<-rho*h_theta(xs,theta)
  ybl<-if(is.null(yline)) boundary_layer_height(ystar,K,'sqrt') else yline
  if(init=='manifold') x0<-1-R0/(R0-1+rho)*y0 else {x0<-1-1/K;y0<-1/K}
  rhs<-function(t,z,p) list(c(rho*h_theta(z[1],theta)-R0*z[1]*z[2],
                                (R0*z[1]-1)*z[2]))
  if(is.null(tmax)) tmax<-max(200,20/rho,2*log(1/y0)/(R0-1))
  times<-seq(0,tmax,by=dt)
  z<-deSolve::ode(c(x=x0,y=y0),times,rhs,NULL,method='lsoda',rtol=2e-10,
                  atol=c(2e-12,2e-14),maxsteps=2e6)
  x<-z[,'x']; y<-z[,'y']; n<-nrow(z)
  peak<-which(x[-n]>=xs & x[-1]<xs)[1]
  if(is.na(peak)) return(list(status='NO_PEAK',x_in=NA_real_,y_min=NA_real_,ybl=ybl))
  ii<-seq.int(peak+1,n-1)
  entry<-ii[y[ii]>=ybl & y[ii+1]<ybl]
  trough<-ii[x[ii]<xs & x[ii+1]>=xs]
  ie<-if(length(entry)) entry[1] else Inf; it<-if(length(trough)) trough[1] else Inf
  if(is.infinite(it)) return(list(status='NO_TROUGH',x_in=NA_real_,y_min=min(y[ii]),ybl=ybl))
  ymin<-min(y[seq.int(peak,it+1)])
  if(it<ie) ans<-list(status='NO_ENTRY',x_in=NA_real_,y_min=ymin,ybl=ybl,
                     y_min_ratio=ymin/ybl)
  else {
    w<-(ybl-y[ie])/(y[ie+1]-y[ie]); xin<-x[ie]+w*(x[ie+1]-x[ie])
    ans<-list(status='OK',x_in=xin,y_min=ymin,ybl=ybl,y_min_ratio=ymin/ybl)
  }
  if(trajectory) ans$trajectory<-as.data.frame(z)
  ans
}
