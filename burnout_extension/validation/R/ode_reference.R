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

# First post-epidemic trough: the first local minimum of y after the epidemic
# peak.  Since y'=(R0*x-1)y and y remains positive in the deterministic model,
# it is the first upward crossing of the infective nullcline x=1/R0 after the
# downward crossing at the peak.  We bracket on the solver output, then reintegrate
# only across that bracket and root-find x(t)-1/R0; y_min is therefore not a
# coarse-grid minimum.
first_deterministic_trough <- function(R0,rho,theta,K=Inf,y0=1e-10,
    init=c('manifold','finiteK'),dt=.1,tmax=NULL,trajectory=FALSE) {
  init<-match.arg(init); xs<-1/R0
  if(init=='finiteK' && (!is.finite(K)||K<=1))
    stop('finiteK initialization requires finite K > 1')
  if(init=='manifold') x0<-1-R0/(R0-1+rho)*y0 else {x0<-1-1/K;y0<-1/K}
  # Log prevalence preserves positivity when the trough is far below an
  # absolute tolerance that could be placed directly on y.
  rhs<-function(t,z,p) list(c(rho*h_theta(z[1],theta)-R0*z[1]*exp(z[2]),
                              R0*z[1]-1))
  if(is.null(tmax)) tmax<-max(200,20/rho,2*log(1/y0)/(R0-1))
  times<-seq(0,tmax,by=dt)
  z<-deSolve::ode(c(x=x0,log_y=log(y0)),times,rhs,NULL,method='lsoda',rtol=2e-10,
                  atol=c(2e-12,2e-11),maxsteps=2e6)
  x<-z[,'x']; n<-nrow(z)
  peak<-which(x[-n]>=xs & x[-1]<xs)[1]
  if(is.na(peak)) return(list(status='NO_PEAK',x_min=NA_real_,
    y_min=NA_real_,t_min=NA_real_,nullcline_residual=NA_real_))
  candidates<-seq.int(peak+1L,n-1L)
  hit<-candidates[x[candidates]<xs & x[candidates+1L]>=xs]
  if(!length(hit)) return(list(status='NO_TROUGH',x_min=NA_real_,
    y_min=exp(min(z[seq.int(peak+1L,n),'log_y'])),t_min=NA_real_,
    nullcline_residual=NA_real_))
  i<-hit[1L]; tlo<-z[i,'time']; thi<-z[i+1L,'time']; state0<-z[i,c('x','log_y')]
  state_at<-function(tt) {
    if(tt==tlo) return(state0)
    drop(deSolve::ode(state0,c(tlo,tt),rhs,NULL,method='lsoda',rtol=5e-12,
      atol=c(5e-14,5e-13),maxsteps=1e5)[2,c('x','log_y')])
  }
  tmin<-uniroot(function(tt) state_at(tt)[1]-xs,c(tlo,thi),
                 tol=max(1e-12,.Machine$double.eps*max(1,thi)))$root
  at_min<-state_at(tmin)
  ans<-list(status='OK',x_min=unname(at_min[1]),y_min=exp(unname(at_min[2])),
            log_y_min=unname(at_min[2]),
            t_min=tmin,nullcline_residual=R0*unname(at_min[1])-1)
  if(trajectory) {
    ans$trajectory<-as.data.frame(z)
    ans$trajectory$y<-exp(ans$trajectory$log_y)
  }
  ans
}
