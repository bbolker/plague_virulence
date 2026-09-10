source('psi_validation/R/theory_psi.R')

# Second-order regular-outer initialization for the exact nonlinear tail.
outer_initial_psi <- local({
  cache<-new.env(parent=emptyenv())
  function(R0,theta,psi,xbar_fraction=.45){
    key<-sprintf('%.12g|%.12g|%.12g|%.6g',R0,theta,psi,xbar_fraction)
    if(exists(key,cache,inherits=FALSE))return(get(key,cache))
    xf<-x_final_psi(R0,psi);xs<-x_star_psi(R0,psi)
    xb<-xf+xbar_fraction*(xs-xf);yb<-F_psi(xb,R0,psi)
    n0<-function(x)n0_psi(x,R0,psi);F<-function(x)F_psi(x,R0,psi)
    Fp<-function(x)n0(x)^psi/(R0*x)-1
    Fpp<-function(x)psi*n0(x)^(2*psi-1)/(R0^2*x^2)-n0(x)^psi/(R0*x^2)
    rhs<-function(x,z,p){
      nn<-n0(x);ff<-F(x);hh<-h_theta(x,theta);Y1<-nn^psi*z[1]
      w1p<-hh*Fp(x)/(R0*x*ff)
      w2p<-hh/(R0^2*x^2*ff^2)*(psi*nn^(psi-1)*ff-(nn^psi-R0*x))*Y1+
        hh^2*nn^psi*(nn^psi-R0*x)/(R0^3*x^3*ff^2)
      list(c(w1p,w2p))
    }
    xhi<-1-1e-8
    z<-deSolve::ode(c(W1=0,W2=0),c(xhi,xb),rhs,NULL,method='lsoda',
      rtol=2e-11,atol=c(2e-13,2e-13),maxsteps=1e6)
    W1<-z[2,'W1'];W2<-z[2,'W2'];nn<-n0(xb)
    Y1<-nn^psi*W1;Y2<-nn^psi*W2+psi/(2*nn)*Y1^2
    W1p<-h_theta(xb,theta)*Fp(xb)/(R0*xb*F(xb))
    n0p<-nn^psi/(R0*xb)
    Y1p<-psi*nn^(psi-1)*n0p*W1+nn^psi*W1p
    X1<- -Y1/Fp(xb)
    X2<- -(0.5*Fpp(xb)*X1^2+Y1p*X1+Y2)/Fp(xb)
    ans<-c(xbar=xb,ybar=yb,X1=unname(X1),X2=unname(X2))
    assign(key,ans,cache);ans
  }
})

# Locate the first recovery trough on the continuous deterministic background.
# This is the finite-prevalence Kendall saddle from Section "Finite-prevalence
# Kendall crossing" of burnout_psi_nonlinear_tail_theory.tex.  Integrating in
# time avoids the inverse-flow singularity at g = 1.
finite_prevalence_trough_psi <- function(R0,rho,theta,psi,
                                        xbar_fraction=.45,
                                        time_max=1e6){
  oi<-outer_initial_psi(R0,theta,psi,xbar_fraction)
  x0<-unname(oi['xbar']+rho*oi['X1']+rho^2*oi['X2'])
  y0<-unname(oi['ybar'])
  empty<-c(t_t=NA_real_,x_t=NA_real_,y_t=NA_real_,log_y_t=NA_real_,
    s_t=NA_real_,alpha_t=NA_real_,gaussian_width=NA_real_)
  if(!is.finite(x0)||!is.finite(y0)||x0<=0||y0<=0)return(empty)
  growth<-function(x,y)R0*x*(x+y)^(-psi)
  if(growth(x0,y0)>=1)return(empty)
  rhs<-function(time,z,parms){
    x<-z[1];y<-exp(z[2]);g<-growth(x,y)
    list(c(rho*h_theta(x,theta)-g*y,g-1))
  }
  rootfun<-function(time,z,parms)growth(z[1],exp(z[2]))-1
  identity_event<-function(time,z,parms)z
  z<-deSolve::ode(c(x=x0,log_y=log(y0)),c(0,time_max),rhs,NULL,method='lsodar',
    rootfunc=rootfun,
    events=list(func=identity_event,root=TRUE,terminalroot=1),
    rtol=2e-10,atol=c(2e-12,2e-11),maxsteps=1e6)
  roots<-attr(z,'troot')
  if(!length(roots)||!is.finite(roots[1]))return(empty)
  last<-z[nrow(z),];xt<-unname(last['x']);logyt<-unname(last['log_y']);yt<-exp(logyt)
  alpha<-(rho*h_theta(xt,theta)-yt)*(1/xt-psi/(xt+yt))
  if(!is.finite(logyt)||!is.finite(alpha)||alpha<=0)return(empty)
  c(t_t=unname(roots[1]),x_t=xt,y_t=yt,log_y_t=logyt,s_t=yt/xt,
    alpha_t=alpha,gaussian_width=1/sqrt(alpha))
}

finite_prevalence_kendall_psi <- function(R0,rho,theta,psi,K_values,
                                          xbar_fraction=.45){
  tr<-finite_prevalence_trough_psi(R0,rho,theta,psi,xbar_fraction)
  logB<-log(K_values)+unname(tr['log_y_t'])+
    .5*log(unname(tr['alpha_t'])/(2*pi))
  B<-exp(pmin(logB,log(.Machine$double.xmax)))
  B[logB>log(.Machine$double.xmax)]<-Inf
  probability<-ifelse(is.infinite(B),1,-expm1(-B))
  data.frame(K=K_values,trough=probability,log_B_trough=logB,
    t_t=unname(tr['t_t']),x_t=unname(tr['x_t']),y_t=unname(tr['y_t']),
    log_y_t=unname(tr['log_y_t']),s_t=unname(tr['s_t']),alpha_t=unname(tr['alpha_t']),
    gaussian_width=unname(tr['gaussian_width']))
}

# Evaluate all requested K values with a single nonlinear-tail integration.
nonlinear_tail_quantities <- function(R0,rho,theta,psi,K_values,
                                      xbar_fraction=.45){
  oi<-outer_initial_psi(R0,theta,psi,xbar_fraction)
  xf<-x_final_psi(R0,psi);xs<-x_star_psi(R0,psi);eta<-1-psi
  hs<-h_theta(xs,theta);ygeo<-sqrt(rho*hs/K_values)
  x<-unname(oi['xbar']+rho*oi['X1']+rho^2*oi['X2']);dl<-.005
  empty_result<-function()data.frame(K=K_values,x_h=NA_real_,y_h=NA_real_,
    C_eff=NA_real_,s=NA_real_,d=NA_real_,epsilon_A=NA_real_,g=NA_real_,
    overlap_score=NA_real_,leading=NA_real_,next_order=NA_real_,
    c_L=laplace_correction_psi(R0,theta,psi),trough=NA_real_,
    log_B_trough=NA_real_,t_t=NA_real_,x_t=NA_real_,y_t=NA_real_,
    log_y_t=NA_real_,s_t=NA_real_,alpha_t=NA_real_,gaussian_width=NA_real_)
  if(!is.finite(x)||x<=xf||x>=xs)return(empty_result())
  ell_max<-max(5,log(oi['ybar']/min(ygeo))+2);nmax<-ceiling(ell_max/dl)
  path<-matrix(NA_real_,nmax+1L,2L,dimnames=list(NULL,c('ell','x')));path[1,]<-c(0,x);n<-1L
  f<-function(l,xx){y<-oi['ybar']*exp(-l);g<-R0*xx*(xx+y)^(-psi)
    (rho*h_theta(xx,theta)-g*y)/(1-g)}
  for(j in seq_len(nmax)){
    l<-(j-1)*dl;y<-oi['ybar']*exp(-l);gg<-R0*x*(x+y)^(-psi)
    if(!is.finite(gg)||gg>=.98||x<=0)break
    k1<-f(l,x);k2<-f(l+dl/2,x+dl*k1/2);k3<-f(l+dl/2,x+dl*k2/2);k4<-f(l+dl,x+dl*k3)
    xn<-x+dl*(k1+2*k2+2*k3+k4)/6
    if(!is.finite(xn)||xn<=0)break
    x<-xn;n<-n+1L;path[n,]<-c(j*dl,x)
  }
  path<-path[seq_len(n),,drop=FALSE];yy<-oi['ybar']*exp(-path[,'ell']);xx<-path[,'x']
  keep<-is.finite(xx)&xx>xf&xx<1
  if(!any(keep))return(empty_result())
  yy<-yy[keep];xx<-xx[keep]
  gg<-R0*xx*(xx+yy)^(-psi);ghh<-R0*xx^eta;ss<-yy/xx
  dd<-gg*yy/(rho*h_theta(xx,theta));ee<-(psi*ghh*ss+abs(1-ghh)*dd)/(1-gg)
  # Use the traditional geometric handoff when it is reached before the safe
  # pre-trough cutoff; otherwise use the deepest available descending-tail
  # point.  This continuous rule avoids turning a discrete diagnostic optimum
  # into artificial kinks in parameter curves.
  choose<-function(target)which.min(abs(log(yy/target)))
  pick<-vapply(ygeo,choose,1L);xh<-xx[pick];yh<-yy[pick];g<-gg[pick]
  s<-ss[pick];d<-dd[pick];epsA<-ee[pick]
  score<-pmax(s/.1,d/.1,epsA/.1,10/(K_values*yh));gh<-R0*xh^eta
  logC<-log(yh)+vapply(xh,function(x)action_diff_psi(xf,x,R0,theta,psi),0.)/rho
  Ceff<-exp(logC);dA<-action_diff_psi(xf,xs,R0,theta,psi)
  cL<-laplace_correction_psi(R0,theta,psi)
  logB<-log(K_values)+logC+.5*log(eta*rho*hs/(2*pi*xs))-dA/rho
  B<-exp(pmin(logB,log(.Machine$double.xmax)));B[logB>log(.Machine$double.xmax)]<-Inf
  lead<-ifelse(is.infinite(B),1,-expm1(-B))
  Bnext<-B*exp(-rho*cL);nextp<-ifelse(is.infinite(Bnext),1,-expm1(-Bnext))
  ans<-data.frame(K=K_values,x_h=xh,y_h=yh,C_eff=Ceff,s=s,d=d,epsilon_A=epsA,g=g,
    overlap_score=score,
    leading=lead,next_order=nextp,c_L=cL)
  tr<-finite_prevalence_kendall_psi(R0,rho,theta,psi,K_values,xbar_fraction)
  merge(ans,tr,by='K',sort=FALSE)
}
