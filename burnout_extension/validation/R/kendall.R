source('validation/R/theory.R')

kendall_quantities <- function(xin,R0,rho,theta,K,yline=NULL) {
  xs<-1/R0
  if(is.null(yline)) yline<-sqrt(rho*h_theta(xs,theta)/K)
  ybl<-yline
  Lambda<-action_diff(xin,xs,R0,theta)/rho
  scaled<-function(X) {
    vapply(X,function(xx) if(xx>=1) 0 else
      exp(action_diff(xs,xx,R0,theta)/rho)/(rho*h_theta(xx,theta)),0.)
  }
  J<-integrate(scaled,xin,1,rel.tol=2e-9,abs.tol=0,subdivisions=1200L,
               stop.on.error=TRUE)$value
  logI<-Lambda+log(J)
  B<-K*ybl*log1p(exp(-logI)); P1<- -expm1(-B)
  logIL<-0.5*log(2*pi/(R0*rho*h_theta(xs,theta)))+Lambda
  BL<-K*ybl*log1p(exp(-logIL)); P1L<- -expm1(-BL)
  list(Lambda=Lambda,logI=logI,B=B,P1=P1,B_L=BL,P1_L=P1L)
}
