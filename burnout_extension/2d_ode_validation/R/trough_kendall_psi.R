# Full two-dimensional deterministic trajectory and finite-horizon Kendall
# validation for burnout_finite_prevalence_trough_theory.tex.

h_theta_2d <- function(x,theta)(1-x)*x^theta

log1pexp_2d <- function(x){
  ifelse(x>35,x,ifelse(x< -35,exp(x),log1p(exp(x))))
}

prob_from_logB_2d <- function(logB){
  if(!is.finite(logB))return(if(is.na(logB))NA_real_ else if(logB>0)1 else 0)
  if(logB>log(.Machine$double.xmax))return(1)
  -expm1(-exp(logB))
}

empty_2d_result <- function(status){
  vals<-c(t_minus=NA_real_,t_t=NA_real_,t_plus=NA_real_,x_minus=NA_real_,
    y_minus=NA_real_,x_t=NA_real_,y_t=NA_real_,log_y_t=NA_real_,
    x_plus=NA_real_,y_plus=NA_real_,alpha_minus=NA_real_,alpha_t=NA_real_,
    alpha_plus=NA_real_,tau_G=NA_real_,
    L_minus=NA_real_,L_plus=NA_real_,Delta_minus=NA_real_,Delta_plus=NA_real_,
    kendall_scaled_integral=NA_real_,log_kendall_integral=NA_real_,
    single_lineage_u=NA_real_,log_B_trough=NA_real_,P_trough=NA_real_,
    P_kendall=NA_real_)
  c(list(status=status),as.list(vals))
}

trough_kendall_2d <- function(R0,rho,theta,psi,K,tmax=NULL){
  stopifnot(R0>1,rho>0,theta>=0,psi>=0,K>1)
  if(!requireNamespace('deSolve',quietly=TRUE))stop('deSolve is required')
  if(is.null(tmax))tmax<-max(500,30/rho,4*log(K)/(R0-1))
  # Log x is the stable default across the full grid. For psi=1, theta=0 the
  # high-R0 node regime has no finite recovery crossing; late near-equilibrium
  # root noise is classified after the solve rather than switching coordinates.
  direct_x<-FALSE
  log_x_of<-function(qx)if(direct_x)log(pmax(qx,.Machine$double.xmin))else qx
  first_atol<-if(direct_x)1e-100 else 1e-10
  log_sum_exp2<-function(a,b){m<-max(a,b);m+log(exp(a-m)+exp(b-m))}
  log_growth<-function(qx,log_y){
    log_x<-log_x_of(qx)
    log_n<-log_sum_exp2(log_x,log_y)
    log(R0)+log_x-psi*log_n
  }
  growth<-function(log_x,log_y)exp(log_growth(log_x,log_y))
  rhs<-function(time,z,parms){
    qx<-z[1];log_x<-log_x_of(qx);log_y<-z[2];x<-exp(log_x)
    log_n<-log_sum_exp2(log_x,log_y);g<-growth(qx,log_y)
    recruitment_per_x<-rho*(1-x)*exp((theta-1)*log_x)
    infection_per_x<-exp(log(R0)+log_y-psi*log_n)
    dlogx<-recruitment_per_x-infection_per_x
    list(c(if(direct_x)x*dlogx else dlogx,g-1))
  }
  # log(g) has exactly the same zeros as g-1, but remains well scaled when
  # both x and y are exponentially small (notably psi=1, theta=0).
  rootfun<-function(time,z,parms)log_growth(z[1],z[2])
  identity_event<-function(time,z,parms)z
  root_once<-function(state,time0){
    solve_root<-function(rtol,atol)tryCatch(
      deSolve::ode(state,c(time0,tmax),rhs,NULL,method='lsodar',
        rootfunc=rootfun,events=list(func=identity_event,root=TRUE,terminalroot=1),
        rtol=rtol,atol=atol,hmax=if(direct_x).1 else 5,maxsteps=2e6),
      error=function(e)NULL)
    z<-solve_root(1e-9,c(first_atol,1e-9))
    if(is.null(z))z<-solve_root(1e-8,c(if(direct_x)1e-90 else 1e-8,1e-8))
    if(is.null(z))return(NULL)
    roots<-attr(z,'troot')
    if(!length(roots)||!is.finite(roots[1]))return(NULL)
    root_time<-unname(roots[1])
    # Some deSolve versions report troot correctly but retain a later output
    # row even with terminalroot set. Reintegrate exactly to the reported root
    # rather than assuming the matrix's last row is the event state.
    zr<-deSolve::ode(state,c(time0,root_time),rhs,NULL,method='lsoda',
      rtol=1e-10,atol=c(if(direct_x)1e-110 else 1e-11,1e-10),
      hmax=if(direct_x).05 else 1,maxsteps=2e6)
    list(time=root_time,state=unname(zr[nrow(zr),c('q_x','log_y')]))
  }
  advance<-function(hit){
    # Move far enough away from an event that DLSODAR does not re-detect the
    # same root as an initial-time zero at extremely deep, sharp valleys.
    eps<-1e-2
    z<-deSolve::ode(setNames(hit$state,c('q_x','log_y')),
      c(hit$time,hit$time+eps),rhs,NULL,method='lsoda',
      rtol=1e-9,atol=c(first_atol,1e-9),maxsteps=1e5)
    list(time=hit$time+eps,state=unname(z[nrow(z),c('q_x','log_y')]))
  }

  initial<-c(q_x=if(direct_x)1-1/K else log1p(-1/K),log_y=log(1/K))
  if(rootfun(0,initial,NULL)<=0)return(empty_2d_result('NO_INITIAL_GROWTH'))
  peak<-root_once(initial,0)
  if(is.null(peak))return(empty_2d_result('NO_FIRST_PEAK'))
  after_peak<-advance(peak)
  trough<-root_once(setNames(after_peak$state,c('q_x','log_y')),after_peak$time)
  if(is.null(trough))return(empty_2d_result('NO_TROUGH'))
  after_trough<-advance(trough)
  next_peak<-root_once(setNames(after_trough$state,c('q_x','log_y')),after_trough$time)
  if(is.null(next_peak))return(empty_2d_result('NO_RECOVERY_PEAK'))

  unpack<-function(hit){
    lx<-log_x_of(hit$state[1])
    c(q_x=hit$state[1],log_x=lx,x=exp(lx),log_y=hit$state[2],y=exp(hit$state[2]))
  }
  zm<-unpack(peak);zt<-unpack(trough);zp<-unpack(next_peak)
  crossing_rate<-function(z){
    log_n<-log_sum_exp2(z['log_x'],z['log_y'])
    w_x<-exp(z['log_x']-log_n)
    dlogx<-rho*(1-z['x'])*exp((theta-1)*z['log_x'])-
      exp(log(R0)+z['log_y']-psi*log_n)
    unname((1-psi*w_x)*dlogx)
  }
  alpha_minus<-crossing_rate(zm);alpha_t<-crossing_rate(zt)
  alpha_plus<-crossing_rate(zp)
  if(!(alpha_minus<0&&alpha_t>0&&alpha_plus<0)){
    bad<-empty_2d_result('BAD_CROSSING_DIRECTION')
    bad$alpha_minus<-alpha_minus;bad$alpha_t<-alpha_t;bad$alpha_plus<-alpha_plus
    bad$t_minus<-peak$time;bad$t_t<-trough$time;bad$t_plus<-next_peak$time
    return(bad)
  }

  # Evaluate the Kendall integral in trough-scaled form.  Its integrand is at
  # most one on the physical first-peak to recovery-peak interval, so even an
  # exponentially deep trough remains numerically stable.
  rhs_quad<-function(time,z,parms){
    qx<-z[1];log_x<-log_x_of(qx);log_y<-z[2];x<-exp(log_x)
    log_n<-log_sum_exp2(log_x,log_y);g<-growth(qx,log_y)
    recruitment_per_x<-rho*(1-x)*exp((theta-1)*log_x)
    infection_per_x<-exp(log(R0)+log_y-psi*log_n)
    dlogx<-recruitment_per_x-infection_per_x
    list(c(if(direct_x)x*dlogx else dlogx,g-1,
      exp(-(log_y-zt['log_y']))))
  }
  q0<-c(q_x=zm['q_x'],log_y=zm['log_y'],J=0)
  qz<-deSolve::ode(q0,c(peak$time,next_peak$time),rhs_quad,NULL,
    method='lsoda',rtol=2e-10,atol=c(if(direct_x)1e-110 else 2e-12,2e-11,2e-11),
    hmax=max(.05,min(1,(1/sqrt(alpha_t))/5)),maxsteps=2e6)
  Jscaled<-unname(qz[nrow(qz),'J'])
  if(!is.finite(Jscaled)||Jscaled<=0)return(empty_2d_result('BAD_KENDALL_INTEGRAL'))
  logI<-unname(zm['log_y']-zt['log_y']+log(Jscaled))
  logu<- -log1pexp_2d(logI)
  u<-if(logu<log(.Machine$double.xmin))0 else exp(logu)
  logm<-log(K)+unname(zm['log_y'])
  if(logu< -35){
    Pk<-prob_from_logB_2d(logm+logu)
  }else{
    m<-exp(logm)
    Pk<- -expm1(m*log1p(-u))
  }
  logBtr<-log(K)+unname(zt['log_y'])+.5*log(alpha_t/(2*pi))
  Ptr<-prob_from_logB_2d(logBtr)
  sqrta<-sqrt(alpha_t)
  ans<-list(status='OK',t_minus=peak$time,t_t=trough$time,t_plus=next_peak$time,
    x_minus=unname(zm['x']),y_minus=unname(zm['y']),x_t=unname(zt['x']),
    y_t=unname(zt['y']),log_y_t=unname(zt['log_y']),x_plus=unname(zp['x']),
    y_plus=unname(zp['y']),alpha_minus=alpha_minus,alpha_t=alpha_t,
    alpha_plus=alpha_plus,tau_G=1/sqrta,
    L_minus=sqrta*(trough$time-peak$time),
    L_plus=sqrta*(next_peak$time-trough$time),
    Delta_minus=unname(zm['log_y']-zt['log_y']),
    Delta_plus=unname(zp['log_y']-zt['log_y']),
    kendall_scaled_integral=Jscaled,log_kendall_integral=logI,
    single_lineage_u=u,log_B_trough=logBtr,P_trough=Ptr,P_kendall=Pk)
  ans
}

trough_kendall_row_2d <- function(R0,rho,theta,psi,K){
  z<-trough_kendall_2d(R0,rho,theta,psi,K)
  as.data.frame(c(list(rho=rho,theta=theta,K=K,R0=R0,psi=psi),z),
    stringsAsFactors=FALSE)
}
