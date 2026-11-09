# Compare the manuscript's direct entry formulas on the original validation grid.
# Reuses saved stochastic counts; does not run simulations or edit paper figures.
source('validation/R/theory.R')
source('validation/R/kendall.R')
source('validation/R/ode_reference.R')
source('validation/R/validation_plot_data.R')
library(data.table)
sim <- load_validation_plot_data()$sim[R0 <= 6]
stopifnot(nrow(sim)==280L)
outfile <- 'validation/data/direct_entry_comparison_original_grid.csv'
# Stable evaluation of the same C_theta quadrature used by the current figures.
C_explicit <- function(R0,theta) exp(unname(C_log_derivatives(R0,theta)['logC']))
constants <- new.env(parent=emptyenv())
rows <- list()
for(i in seq_len(nrow(sim))) {
 z<-sim[i];r<-z$R0;rho<-z$rho;th<-z$theta;K<-z$K
 key<-paste(r,th)
 w<-character()
 ans<-tryCatch(withCallingHandlers({
  if(!exists(key,constants,inherits=FALSE))assign(key,D_regularized(r,th),constants)
  D<-get(key,constants);C<-C_explicit(r,th);xf<-x_final(r);delta<-1-r*xf
  hf<-h_theta(xf,th);a<-hf/delta;b<-r*xf/delta
  cc<-delta*h_prime(xf,th)+r*hf;kappa<-a*cc/(2*delta^2)
  ystar<-rho*h_theta(1/r,th);yb<-sqrt(ystar/K);L<-log(C/yb)
  x1<-xf+rho*a*L+b*yb
  x2<-x1+rho^2*(kappa*L^2+D)+rho*yb/delta^2*(r*a*(L+1)-b*cc)+b*r*yb^2/(2*delta^2)
  ode<-ode_reference(r,rho,th,K,yline=yb)
  valid<-function(x)is.finite(x)&&x>xf&&x<1/r
  result<-list(status='OK',ODE_status=ode$status,C=C,D=D,ybl=yb,rhoL2=rho*L^2,
   xf=xf,x_first=x1,x_second=x2,x_ode=ode$x_in,
   first_admissible=valid(x1),second_admissible=valid(x2),
   P_BI=stable_validation_probability(r,rho,th,K),
   P_first_K=NA_real_,P_second_K=NA_real_,P_first_L=NA_real_,P_second_L=NA_real_,P_ODE_K=NA_real_,P_ODE_L=NA_real_)
  for(pair in list(list('first',x1),list('second',x2),list('ODE',ode$x_in))) {
   nm<-pair[[1]];x<-pair[[2]]
   if(valid(x)){
    q<-tryCatch(kendall_quantities(x,r,rho,th,K,yline=yb),error=function(e)e)
    if(inherits(q,'error'))result[[paste0(nm,'_integration_error')]]<-conditionMessage(q)
    else{result[[paste0('P_',nm,'_K')]]<-q$P1;result[[paste0('P_',nm,'_L')]]<-q$P1_L}
   }
  }
  result
 },warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)list(status=conditionMessage(e)))
 ans$warnings<-paste(unique(w),collapse=' | ')
 rows[[i]]<-cbind(z,as.data.table(ans))
 if(i%%20==0 || i==nrow(sim)){
  fwrite(rbindlist(rows,fill=TRUE),outfile)
  cat(i,'/',nrow(sim),'completed\n');flush.console()
 }
}
d<-rbindlist(rows,fill=TRUE)
# Reproduce the previously inspected direct-formula calculations before summarizing.
old<-fread('validation/data/fig1_right_corner_theory_comparison.csv')
check<-merge(d,old,by=c('R0','rho','theta','K'),suffixes=c('_new','_old'))
stopifnot(nrow(check)==3L)
for(nm in c('x_first','x_second','P_BI','P_first_K','P_second_K'))
 stopifnot(max(abs(check[[paste0(nm,'_new')]]-check[[paste0(nm,'_old')]]))<2e-5)
common<-d[ODE_status=='OK' & is.finite(P_first_K)&is.finite(P_second_K)&is.finite(P_BI)]
methods<-c('P_BI','P_first_K','P_second_K','P_first_L','P_second_L','P_ODE_K','P_ODE_L')
summary<-rbindlist(lapply(methods,function(nm){e<-abs(common[[nm]]-common$P_conditional);data.table(method=nm,n=sum(is.finite(e)),MAE=mean(e,na.rm=TRUE),median=median(e,na.rm=TRUE),max=max(e,na.rm=TRUE))}))
fwrite(summary,'validation/data/direct_entry_comparison_summary.csv')
cat('Common valid grid:',nrow(common),'of',nrow(d),'\n');print(summary)
cat('ODE status\n');print(d[,.N,by=ODE_status])
cat('Numerical exceptions\n');print(d[status!='OK',.(R0,theta,K,status)])
cat('Completed:',outfile,'\n')
