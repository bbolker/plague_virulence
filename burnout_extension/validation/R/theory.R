# General-recruitment burnout asymptotics. Base R + deSolve only.
h_theta <- function(x, theta) (1-x)*x^theta
h_prime <- function(x, theta) x^(theta-1)*(theta-(theta+1)*x)
F0 <- function(x, R0) 1-x+log(x)/R0
Fp <- function(x, R0) 1/(R0*x)-1
Fpp <- function(x, R0) -1/(R0*x^2)

x_final <- function(R0) {
  uniroot(function(x) F0(x,R0), c(1e-15, 1/R0-1e-12),
          tol=1e-13)$root
}

action_derivative <- function(x,R0,theta) (1-R0*x)/h_theta(x,theta)
action_diff <- function(x0,x1,R0,theta) {
  if (x0==x1) return(0)
  integrate(function(z) action_derivative(z,R0,theta),x0,x1,
            rel.tol=2e-11,abs.tol=2e-12,subdivisions=500L,stop.on.error=TRUE)$value
}

C_explicit <- local({
  cache <- new.env(parent=emptyenv())
  function(R0,theta) {
    key <- sprintf('%.14g|%.14g',R0,theta)
    if (exists(key,cache,inherits=FALSE)) return(get(key,cache))
    xf <- x_final(R0); xs <- 1/R0; at <- h_theta(xf,theta)/xf
    raw <- function(u) (xs-u)*h_theta(u,theta)/(u^2*F0(u,R0))-at/(u-xf)
    # Estimate removable endpoint value by Richardson extrapolation.
    ee <- max(1e-7*xf,1e-9)
    endpoint <- 2*raw(xf+ee/2)-raw(xf+ee)
    J <- endpoint*ee/4 + integrate(raw,xf+ee/4,1,rel.tol=2e-10,
          abs.tol=2e-11,subdivisions=1000L)$value
    ans <- (1-xf)*(xs-xf)/xf*exp(J/at)
    assign(key,ans,cache); ans
  }
})

outer_coefficients_at <- function(R0,theta,xvals) {
    stopifnot(requireNamespace('deSolve',quietly=TRUE))
    xhi <- 1-1e-10
    rhs <- function(x,z,p) {
      ff<-F0(x,R0); hh<-h_theta(x,theta)
      y1p <- -hh*(R0*x-1)/(R0^2*x^2*ff)
      y2p <- hh*(R0*x-1)/(R0^2*x^2*ff^2)*(z[1]-hh/(R0*x))
      list(c(y1p,y2p))
    }
    times <- c(xhi,sort(unique(xvals),decreasing=TRUE))
    zz <- deSolve::ode(c(Y1=0,Y2=0),times,rhs,NULL,method='lsoda',
                       rtol=3e-11,atol=c(2e-13,2e-13),maxsteps=1e6)
    ans <- zz[-1,,drop=FALSE]
    ans[match(xvals,ans[,1]),,drop=FALSE]
}

inverse_coefficients_many <- function(ys,R0,theta) {
  xf <- x_final(R0)
  x0 <- vapply(ys,function(y) uniroot(function(x) F0(x,R0)-y,
                                      c(xf,1/R0),tol=1e-13)$root,0.)
  zz <- outer_coefficients_at(R0,theta,x0); y1<-zz[,'Y1']; y2<-zz[,'Y2']
  y1p <- -h_theta(x0,theta)*(R0*x0-1)/(R0^2*x0^2*F0(x0,R0))
  x1 <- -y1/Fp(x0,R0)
  x2 <- -(0.5*Fpp(x0,R0)*x1^2+y1p*x1+y2)/Fp(x0,R0)
  cbind(X0=x0,X1=x1,X2=x2)
}

inverse_coefficients <- function(y,R0,theta) {
  drop(inverse_coefficients_many(y,R0,theta))
}

C_inverse <- function(R0,theta,ys=c(2e-6,1e-6,5e-7,2e-7)) {
  xf<-x_final(R0); a<-h_theta(xf,theta)/(1-R0*xf)
  inv<-inverse_coefficients_many(ys,R0,theta)
  v<-inv[,'X1']/a+log(ys)
  exp(coef(lm(v~I(ys*log(ys))))[1])
}

D_extract <- function(R0,theta,ys=exp(seq(log(2e-3),log(3e-5),length.out=24)),
                      basis=c('ylog2','ylog','sqrt'),n=NULL) {
  basis<-match.arg(basis); xf<-x_final(R0)
  # Near threshold the zero-order epidemic peak can lie below the default window.
  ymax <- F0(1/R0,R0)
  if(max(ys)>=.8*ymax) {
    hi <- .5*ymax; lo <- hi/70
    ys <- exp(seq(log(hi),log(lo),length.out=length(ys)))
  }
  delta<-1-R0*xf; hf<-h_theta(xf,theta); a<-hf/delta
  cc<-delta*h_prime(xf,theta)+R0*hf; kappa<-hf*cc/(2*delta^3)
  C<-C_explicit(R0,theta)
  inv<-inverse_coefficients_many(ys,R0,theta)
  val<-inv[,'X2']-kappa*log(C/ys)^2
  L<-log(ys)
  X<-switch(basis,ylog2=cbind(1,ys*L^2,ys*L,ys),
            ylog=cbind(1,ys*L,ys),sqrt=cbind(1,sqrt(ys),ys*L^2,ys*L))
  fit<-lm.fit(X,val); list(D=unname(fit$coefficients[1]),
    rmse=sqrt(mean(fit$residuals^2)),ys=ys,raw=val,basis=basis,n=n)
}

matching <- function(R0,rho,theta,K,D=NULL) {
  xf<-x_final(R0); xs<-1/R0; delta<-1-R0*xf; hf<-h_theta(xf,theta)
  a<-hf/delta; b<-R0*xf/delta; cc<-delta*h_prime(xf,theta)+R0*hf
  C<-C_explicit(R0,theta); if(is.null(D)) D<-D_extract(R0,theta)$D
  ys<-rho*h_theta(xs,theta); ybl<-sqrt(ys/K); L<-log(C/ybl)
  Q<-(R0*a-b*cc)/(a*delta^2)
  M1<-rho*L+(b/a)*ybl
  M2<-M1+rho^2*D/a+rho*ybl*Q*(L+1)+(b*Q/(2*a))*ybl^2
  total<-action_diff(xf,xs,R0,theta)
  root<-function(M) if(!is.finite(M)||M<=0||M>=total) NA_real_ else
    uniroot(function(x) action_diff(xf,x,R0,theta)-M,c(xf+1e-12,xs-1e-12),tol=2e-11)$root
  list(xf=xf,ystar=ys,ybl=ybl,C=C,D=D,M1=M1,M2=M2,x_first=root(M1),
       x_second=root(M2),Lambda_first=(total-M1)/rho,Lambda_second=(total-M2)/rho)
}
