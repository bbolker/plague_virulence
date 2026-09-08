# Boundary-independent burnout theory for population-size-dependent incidence.
source('validation/R/theory.R')

n0_psi <- function(x,R0,psi) {
  eta <- 1-psi
  z <- 1+eta/R0*log(x)
  ifelse(z>0,z^(1/eta),NA_real_)
}
F_psi <- function(x,R0,psi) n0_psi(x,R0,psi)-x

x_star_psi <- function(R0,psi) R0^(-1/(1-psi))
x_final_psi <- function(R0,psi) {
  eta <- 1-psi; xs <- x_star_psi(R0,psi)
  lo <- exp(-R0/eta)*(1+1e-10)
  uniroot(function(x) x^eta-1-eta/R0*log(x),c(lo,xs),tol=2e-13)$root
}

action_diff_psi <- function(x0,x1,R0,theta,psi) {
  if(x0==x1) return(0)
  eta <- 1-psi
  integrate(function(z) (1-R0*z^eta)/h_theta(z,theta),x0,x1,
    rel.tol=2e-10,abs.tol=2e-12,subdivisions=1000L,stop.on.error=TRUE)$value
}

laplace_correction_psi <- function(R0,theta,psi) {
  eta <- 1-psi; xs <- x_star_psi(R0,psi)
  hs <- h_theta(xs,theta); hp <- h_prime(xs,theta); hpp <- h_second(xs,theta)
  ((2*eta^2-eta-1)*hs^2+(eta-1)*xs*hs*hp+2*xs^2*hp^2-
     3*xs^2*hs*hpp)/(24*eta*hs*xs)
}

# Compute C and D from the inverse-flow coefficient hierarchy.  Several small
# endpoints are used and the known O(y log y), O(y log^2 y) remainders are
# removed by a local fit.  xbar is deliberately well inside (xf,xstar).
CD_psi <- local({
  cache <- new.env(parent=emptyenv())
  function(R0,theta,psi,xbar_fraction=.45) {
    stopifnot(requireNamespace('deSolve',quietly=TRUE))
    key <- sprintf('%.12g|%.12g|%.12g|%.6g',R0,theta,psi,xbar_fraction)
    if(exists(key,cache,inherits=FALSE)) return(get(key,cache))
    if(abs(psi)<1e-14) {
      ans <- c(C=C_explicit(R0,theta),D=D_regularized(R0,theta))
      assign(key,ans,cache); return(ans)
    }
    xf <- x_final_psi(R0,psi); xs <- x_star_psi(R0,psi)
    xb <- xf+xbar_fraction*(xs-xf); yb <- F_psi(xb,R0,psi)
    n0 <- function(x) n0_psi(x,R0,psi)
    F <- function(x) F_psi(x,R0,psi)
    Fp <- function(x) n0(x)^psi/(R0*x)-1
    Fpp <- function(x) psi*n0(x)^(2*psi-1)/(R0^2*x^2)-n0(x)^psi/(R0*x^2)
    outer_rhs <- function(x,z,p) {
      nn<-n0(x); ff<-F(x); hh<-h_theta(x,theta); y1<-nn^psi*z[1]
      w1p<-hh*Fp(x)/(R0*x*ff)
      w2p<-hh/(R0^2*x^2*ff^2)*
        (psi*nn^(psi-1)*ff-(nn^psi-R0*x))*y1+
        hh^2*nn^psi*(nn^psi-R0*x)/(R0^3*x^3*ff^2)
      list(c(w1p,w2p))
    }
    xhi <- 1-1e-8
    oz <- deSolve::ode(c(W1=0,W2=0),c(xhi,xb),outer_rhs,NULL,method='lsoda',
      rtol=2e-11,atol=c(2e-13,2e-13),maxsteps=1e6)
    W1 <- oz[2,'W1']; W2 <- oz[2,'W2']; nn <- n0(xb)
    Y1 <- nn^psi*W1; Y2 <- nn^psi*W2+psi/(2*nn)*Y1^2
    y1p <- h_theta(xb,theta)*Fp(xb)/(R0*xb*F(xb))
    Y1p <- psi*nn^(psi-1)*(nn^psi/(R0*xb))*W1+nn^psi*y1p
    X1b <- -Y1/Fp(xb)
    X2b <- -(0.5*Fpp(xb)*X1b^2+Y1p*X1b+Y2)/Fp(xb)
    inv_rhs <- function(y,z,p) {
      x<-z[1]; x1<-z[2]; x2<-z[3]
      g<-R0*x*(x+y)^(-psi)
      gx<-g*(1/x-psi/(x+y))
      gxx<-gx*(1/x-psi/(x+y))+g*(-1/x^2+psi/(x+y)^2)
      phix<-gx/(g-1)^2
      phixx<-(gxx*(g-1)-2*gx^2)/(g-1)^3
      psix<-(h_prime(x,theta)*(g-1)-h_theta(x,theta)*gx)/(g-1)^2
      list(c(-g/(g-1),phix*x1+h_theta(x,theta)/((g-1)*y),
        phix*x2+.5*phixx*x1^2+psix*x1/y))
    }
    eps <- yb*c(2e-4,1e-4,5e-5,2e-5,1e-5)
    iz <- deSolve::ode(c(X0=unname(xb),X1=unname(X1b),X2=unname(X2b)),c(yb,sort(eps,decreasing=TRUE)),
      inv_rhs,NULL,method='lsoda',rtol=2e-11,atol=c(2e-13,2e-11,2e-9),maxsteps=1e6)
    iz <- iz[-1,,drop=FALSE]; yy<-iz[,'time']; q<-R0*xf^(1-psi);delta<-1-q
    a<-h_theta(xf,theta)/delta
    z1<-iz[,'X1']+a*log(yy)
    C<-exp(coef(lm(z1~I(yy*log(yy))+yy))[1]/a)
    ccoef<-delta*h_prime(xf,theta)+(1-psi)*q/xf*h_theta(xf,theta)
    kappa<-h_theta(xf,theta)*ccoef/(2*delta^3)
    z2<-iz[,'X2']-kappa*log(C/yy)^2
    D<-coef(lm(z2~I(yy*log(yy)^2)+I(yy*log(yy))+yy))[1]
    ans<-c(C=unname(C),D=unname(D));assign(key,ans,cache);ans
  }
})

bi_quantities_psi <- function(R0,rho,theta,psi,K,I0=1,next_order=FALSE) {
  eta<-1-psi;xf<-x_final_psi(R0,psi);xs<-x_star_psi(R0,psi)
  cd<-CD_psi(R0,theta,psi);q<-R0*xf^eta;delta<-1-q
  a<-h_theta(xf,theta)/delta
  dA<-action_diff_psi(xf,xs,R0,theta,psi)
  cL<-laplace_correction_psi(R0,theta,psi)
  shift<-if(next_order) rho*(cd['D']/a-cL) else 0
  logB<-log(K)+log(cd['C'])+.5*log(eta*rho*h_theta(xs,theta)/(2*pi*xs))-
    dA/rho+shift
  B<-if(logB>log(.Machine$double.xmax)) Inf else exp(logB)
  pc<-if(is.infinite(B)) 1 else -expm1(-B)
  est<- -expm1(-I0*log(R0))
  c(P_conditional=unname(pc),P_unconditional=unname(est*pc),log_B=unname(logB),
    C=unname(cd['C']),D=unname(cd['D']),Delta_A_f=dA,c_L=cL,a=a)
}
