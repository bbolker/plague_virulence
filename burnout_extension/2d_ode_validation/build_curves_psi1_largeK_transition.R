source('2d_ode_validation/R/trough_kendall_psi.R')
stopifnot(requireNamespace('data.table',quietly=TRUE))
library(data.table)

grid<-fread('2d_ode_validation/data/psi1_largeK_transition_pilot.csv')[,
  .(rho,theta,K,R0,psi,q,trough_scale)]
ans<-vector('list',nrow(grid))
for(i in seq_len(nrow(grid))){
  g<-grid[i]
  z<-tryCatch(trough_kendall_row_2d(g$R0,g$rho,g$theta,g$psi,g$K),
    error=function(e)data.frame(rho=g$rho,theta=g$theta,K=g$K,R0=g$R0,
      psi=g$psi,status=paste0('ERROR: ',conditionMessage(e))))
  ans[[i]]<-cbind(data.frame(q=g$q,trough_scale=g$trough_scale),z)
}
out<-rbindlist(ans,fill=TRUE)
fwrite(out,'2d_ode_validation/data/psi1_largeK_transition_2d_ode.csv')
print(out[,.(points=.N,OK=sum(status=='OK'),
  statuses=paste(names(table(status)),table(status),collapse='; ')),by=.(rho,K)])
