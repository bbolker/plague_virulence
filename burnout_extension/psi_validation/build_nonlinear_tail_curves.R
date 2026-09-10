source('psi_validation/R/theory_nonlinear_tail.R')
library(data.table)
psi<-0.5;K_values<-c(1000,3000,10000,30000)
grid<-expand.grid(rho=c(.01,.02,.05,.10),theta=c(0,.5,1),
  R0=1+exp(seq(log(.05),log(5),length.out=101)),KEEP.OUT.ATTRS=FALSE)
ans<-vector('list',nrow(grid))
for(i in seq_len(nrow(grid))){
  gr<-grid[i,];rho0<-gr$rho;theta0<-gr$theta;R00<-gr$R0
  z<-try(nonlinear_tail_quantities(R00,rho0,theta0,psi,K_values),silent=TRUE)
  if(inherits(z,'try-error'))stop(sprintf('curve failed at row %d: rho=%g theta=%g R0=%g: %s',
    i,rho0,theta0,R00,as.character(z)))
  setDT(z)
  z[,`:=`(rho=rho0,theta=theta0,R0=R00,psi=psi)]
  ans[[i]]<-z
}
curves<-rbindlist(ans)
setcolorder(curves,c('rho','theta','K','R0','psi','leading','next_order','trough',
  'x_h','y_h','C_eff','s','d','epsilon_A','g','overlap_score','c_L',
  'log_B_trough','t_t','x_t','y_t','log_y_t','s_t','alpha_t','gaussian_width'))
fwrite(curves,'psi_validation/data/psi05_nonlinear_tail_curves.csv')
cat('nonlinear-tail curves complete\n')
