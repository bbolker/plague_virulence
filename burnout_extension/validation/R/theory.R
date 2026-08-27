# General-recruitment burnout asymptotics. Base R + deSolve only.
h_theta <- function(x, theta) (1-x)*x^theta
h_prime <- function(x, theta) x^(theta-1)*(theta-(theta+1)*x)

# Boundary-layer matching heights.  The logarithmic compromises are
# y_BL = K^(-1/3) * y_star^(2/3) and K^(-1/4) * y_star^(3/4), with no
# extra theta correction.
boundary_layer_height <- function(ystar,K,
    choice=c('sqrt','compromise','compromise_3_4','ystar')) {
  choice <- match.arg(choice)
  switch(choice,sqrt=sqrt(ystar/K),
    compromise=K^(-1/3)*ystar^(2/3),
    compromise_3_4=K^(-1/4)*ystar^(3/4),ystar=ystar)
}

boundary_layer_labels <- c(
  sqrt='sqrt(y*/K)',compromise='2/3 compromise',
  compromise_3_4='3/4 compromise',ystar='y*')
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
    # Integrate the removable endpoint from a local polynomial fit; evaluating
    # the subtraction at machine-scale offsets needlessly loses digits.
    ee <- max(min(1e-4*(1-xf),.03*xf),1e-10)
    ss <- ee*c(1,1.4,2,2.8,4)
    fit <- lm.fit(cbind(1,ss,ss^2,ss^3),raw(xf+ss))$coefficients
    endpoint <- fit[1]*ee+fit[2]*ee^2/2+fit[3]*ee^3/3+fit[4]*ee^4/4
    J <- endpoint + integrate(raw,xf+ee,1,rel.tol=2e-11,
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

D_regularized <- local({
  cache <- new.env(parent=emptyenv())
  function(R0,theta) {
    stopifnot(requireNamespace('deSolve',quietly=TRUE))
    key <- sprintf('%.14g|%.14g',R0,theta)
    if (exists(key,cache,inherits=FALSE)) return(get(key,cache))

    xf <- x_final(R0); span <- 1-xf; delta <- 1-R0*xf
    hf <- h_theta(xf,theta); hfp <- h_prime(xf,theta)
    f1 <- delta/(R0*xf); f2 <- -1/(R0*xf^2)
    alpha <- hf/(R0*xf); a <- hf/delta
    C <- C_explicit(R0,theta); beta <- alpha*log(f1/C)
    p0 <- (2*delta*xf*hfp-(1+2*delta)*hf)/(2*R0*delta*xf^2)
    cc <- delta*hfp+R0*hf; kappa <- hf*cc/(2*delta^3)
    Am1 <- (hf/xf-cc)/delta^2
    upper_eps <- 1e-6
    xhi <- 1-upper_eps
    upper_slope <- (R0-1.5)/(R0-1)-theta
    y1hi <- upper_eps/R0+upper_slope*upper_eps^2/(2*R0)
    beta_s <- min(xf,span)*seq(.002,.016,length.out=10)
    y1_rhs <- function(x,z,p) list(-h_theta(x,theta)*(R0*x-1)/
      (R0^2*x^2*F0(x,R0)))
    pilot <- deSolve::ode(c(Y1=y1hi),c(xhi,sort(xf+beta_s,decreasing=TRUE)),
      y1_rhs,NULL,method='lsoda',rtol=5e-13,atol=2e-14,maxsteps=1e6)[-1,,drop=FALSE]
    py <- pilot[match(xf+beta_s,pilot[,1]),'Y1']-alpha*log(beta_s)
    tt <- beta_s/span
    beta_ode <- lm.fit(cbind(1,tt,tt^2,tt^3,tt^4,tt^5),py)$coefficients[1]
    y1_shift <- beta-beta_ode
    B0 <- beta-alpha
    B1 <- -hf/(2*R0*delta*xf^2)

    qsing <- function(s) -a*alpha*log(s)/s^2-a*B0/s^2+
      Am1*alpha*log(s)/s+(Am1*B0-a*B1)/s
    P2 <- function(s) (a*alpha*log(s)+a*beta)/s+
      .5*Am1*alpha*log(s)^2+(Am1*B0-a*B1)*log(s)

    # Stop before the endpoint, where direct subtraction loses digits.  The
    # omitted piece is integrated from its local log-power expansion below.
    eps <- max(min(4e-4*span,.03*xf),2e-8)
    sample_s <- eps*c(1,1.35,1.8,2.5,3.4,4.6,6.2,8.4)
    rhs <- function(x,z,p) {
      ff <- F0(x,R0); hh <- h_theta(x,theta)
      y1p <- -hh*(R0*x-1)/(R0^2*x^2*ff)
      q2 <- hh*(R0*x-1)/(R0^2*x^2*ff^2)*(z[1]+y1_shift-hh/(R0*x))
      list(c(y1p,q2-qsing(x-xf)))
    }
    targets <- sort(unique(c(xf+eps,xf+sample_s)),decreasing=TRUE)
    # Y1(1-e)=e/R0+O(e^2); the resulting O(e^2) initialization error has an
    # O(e) effect on the regularized integral.
    zz <- deSolve::ode(c(Y1=y1hi,Ireg=0),c(xhi,targets),rhs,NULL,method='lsoda',
      rtol=2e-12,atol=c(2e-14,2e-12),maxsteps=1e6)
    zz <- zz[-1,,drop=FALSE]
    ord <- match(xf+sample_s,zz[,1])
    y1 <- zz[ord,'Y1']+y1_shift; xx <- xf+sample_s
    q2 <- h_theta(xx,theta)*(R0*xx-1)/(R0^2*xx^2*F0(xx,R0)^2)*
      (y1-h_theta(xx,theta)/(R0*xx))
    qr <- q2-qsing(sample_s)
    X <- cbind(log(sample_s),1,sample_s*log(sample_s),sample_s,
               sample_s^2*log(sample_s),sample_s^2)
    cf <- lm.fit(X,qr)$coefficients
    endpoint <- cf[1]*eps*(log(eps)-1)+cf[2]*eps+
      cf[3]*eps^2*(log(eps)/2-1/4)+cf[4]*eps^2/2+
      cf[5]*eps^3*(log(eps)/3-1/9)+cf[6]*eps^3/3
    I_from_eps <- -zz[match(xf+eps,zz[,1]),'Ireg']
    # The missing interval (xhi,1) is O(1-xhi) and below the requested accuracy.
    K2 <- -P2(span)-I_from_eps-endpoint
    N20 <- K2+alpha*beta*f2/f1^2-(alpha+beta)*p0/f1+
      f2*beta^2/(2*f1^2)
    ans <- -N20/f1-kappa*beta^2/alpha^2
    assign(key,ans,cache); ans
  }
})

matching <- function(R0,rho,theta,K,D=NULL,yline=NULL) {
  xf<-x_final(R0); xs<-1/R0; delta<-1-R0*xf; hf<-h_theta(xf,theta)
  a<-hf/delta; b<-R0*xf/delta; cc<-delta*h_prime(xf,theta)+R0*hf
  C<-C_explicit(R0,theta); if(is.null(D)) D<-D_regularized(R0,theta)
  ys<-rho*h_theta(xs,theta)
  if(is.null(yline)) yline<-boundary_layer_height(ys,K,'sqrt')
  ybl<-yline; L<-log(C/ybl)
  Q<-(R0*a-b*cc)/(a*delta^2)
  M1<-rho*L+(b/a)*ybl
  M2<-M1+rho^2*D/a+rho*ybl*Q*(L+1)+(b*Q/(2*a))*ybl^2
  total<-action_diff(xf,xs,R0,theta)
  root<-function(M) if(!is.finite(M)||M<=0||M>=total) NA_real_ else
    uniroot(function(x) action_diff(xf,x,R0,theta)-M,c(xf+1e-12,xs-1e-12),tol=2e-11)$root
  list(xf=xf,ystar=ys,ybl=ybl,C=C,D=D,M1=M1,M2=M2,x_first=root(M1),
       x_second=root(M2),Lambda_first=(total-M1)/rho,Lambda_second=(total-M2)/rho)
}
