library(data.table)
simdt<-fread('psi_validation/data/psi05_stochastic_results.csv')
curves<-fread('psi_validation/data/psi05_nonlinear_tail_curves.csv')
setkey(curves,rho,theta,K,R0);setkey(simdt,rho,theta,K,R0)
# Evaluate summaries at the stochastic R0 grid using linear interpolation of
# the dense curves on log(R0-1).
at<-curves[,.(R0_sim=sort(unique(simdt[rho==.BY$rho&theta==.BY$theta&K==.BY$K,R0]))),by=.(rho,theta,K)]
at[,`:=`(leading=approx(log(curves[rho==.BY$rho&theta==.BY$theta&K==.BY$K,R0]-1),
  curves[rho==.BY$rho&theta==.BY$theta&K==.BY$K,leading],log(R0_sim-1),rule=1)$y,
  next_order=approx(log(curves[rho==.BY$rho&theta==.BY$theta&K==.BY$K,R0]-1),
  curves[rho==.BY$rho&theta==.BY$theta&K==.BY$K,next_order],log(R0_sim-1),rule=1)$y),by=.(rho,theta,K)]
setnames(at,'R0_sim','R0');z<-merge(simdt,at,by=c('rho','theta','K','R0'))
out<-z[,.(points=.N,defined_leading=sum(is.finite(leading)),
  MAE_conditional=mean(abs(P_conditional-leading),na.rm=TRUE),
  MAE_next_conditional=mean(abs(P_conditional-next_order),na.rm=TRUE))]
diag<-curves[,.(fraction_nominal_overlap=mean(overlap_score<=1,na.rm=TRUE),
  median_overlap_score=median(overlap_score,na.rm=TRUE),max_epsilon_A=max(epsilon_A,na.rm=TRUE))]
fwrite(cbind(out,diag),'psi_validation/data/psi05_nonlinear_tail_summary.csv');print(cbind(out,diag))
