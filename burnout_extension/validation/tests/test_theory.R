source('validation/R/theory.R')
for(th in c(0,.5,1)) for(x in c(.13,.37,.71)) {
  e<-1e-6; got<-action_diff(x-e,x+e,3,th)/(2*e)
  stopifnot(abs(got-action_derivative(x,3,th))<2e-7)
}
expected<-c('2'=0.07968023,'3'=0.14061515,'5'=0.16570580)
for(nm in names(expected)) stopifnot(abs(C_explicit(as.numeric(nm),.5)-expected[nm])<3e-6)
expected_D<-c('2'=-0.04784,'3'=-0.07927,'5'=-0.02572)
for(nm in names(expected_D))
  stopifnot(abs(D_regularized(as.numeric(nm),.5)-expected_D[nm])<8e-5)
z<-matching(3,.01,.5,1e4,D=-.07962)
stopifnot(abs(action_diff(z$x_second,1/3,3,.5)/.01-z$Lambda_second)<2e-8)
bi<-bi_quantities(3,.01,.5,1e4)
stopifnot(is.finite(bi$log_y_min),bi$P_conditional>=0,bi$P_conditional<=1,
          abs(bi$P_unconditional-(1-1/3)*bi$P_conditional)<1e-14)
bi4<-bi_quantities(3,.01,.5,1e4,I0=4)
stopifnot(abs(bi4$P_unconditional-(1-3^-4)*bi4$P_conditional)<1e-14)
for(R0 in c(1.05,1.2,2,5))
  stopifnot(abs(laplace_correction(R0,0)-1/(12*(R0-1)))<2e-13)
bn<-bi_next_quantities(3,.01,.5,1e4)
stopifnot(all(is.finite(bn$log_B)),all(bn$P_conditional>=0),
          all(bn$P_conditional<=1),
          abs((bn$log_B['D_only']-bn$log_B['leading'])-.01*bn$D/bn$a)<1e-14,
          abs((bn$log_B['Laplace_only']-bn$log_B['leading'])+.01*bn$c_L)<1e-14,
          abs(bn$log_B['combined']-bn$log_B['D_only']+.01*bn$c_L)<1e-14)
bn0<-bi_next_quantities(3,1e-9,.5,1e4,D=bn$D)
stopifnot(abs(bn0$log_B['combined']-bn0$log_B['leading'])<1e-8)
ys<-.01*h_theta(1/3,.5); hs<-boundary_layer_height(ys,1e4,c('sqrt'))
hc<-boundary_layer_height(ys,1e4,'compromise'); hy<-boundary_layer_height(ys,1e4,'ystar')
hc34<-boundary_layer_height(ys,1e4,'compromise_3_4')
stopifnot(abs(log(hc)-(2*log(hs)+log(hy))/3)<1e-12,hs<hc,hc<hy)
stopifnot(abs(log(hc34)-(log(hs)+log(hy))/2)<1e-12,hc<hc34,hc34<hy)
stopifnot(abs(boundary_layer_height(ys,1e4,'log_inverse',2)-2/log(1e4))<1e-14)
stopifnot(abs(boundary_layer_height(ys,1e4,'ystar_log_inverse',2)-2*ys/log(1e4))<1e-14)
source('validation/R/ode_reference.R')
tr<-first_deterministic_trough(3,.01,.5,1e4,dt=.2)
stopifnot(tr$status=='OK',abs(tr$x_min-1/3)<2e-10,abs(tr$nullcline_residual)<1e-9,
          is.finite(tr$y_min),tr$y_min>0,is.finite(tr$t_min))
cat('theory unit tests passed\n')
