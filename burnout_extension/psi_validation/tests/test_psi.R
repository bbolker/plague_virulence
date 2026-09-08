source('psi_validation/R/theory_psi.R')
source('psi_validation/R/stochastic_psi.R')
for(theta in c(0,.5,1))for(R0 in c(1.2,2,4)){
  z<-CD_psi(R0,theta,0)
  stopifnot(abs(z['C']/C_explicit(R0,theta)-1)<1e-12,
            abs(z['D']-D_regularized(R0,theta))<1e-12)
}
z<-bi_quantities_psi(1.2,.02,.5,.5,10000)
stopifnot(is.finite(z['P_conditional']),z['P_conditional']>0,z['P_conditional']<1)
cat('psi theory regression tests passed\n')
