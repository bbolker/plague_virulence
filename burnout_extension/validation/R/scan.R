source('validation/R/theory.R')
source('validation/R/ode_reference.R')
source('validation/R/kendall.R')

run_point <- function(R0,rho,theta,K,D=NULL) {
  ode<-ode_reference(R0,rho,theta,K); m<-matching(R0,rho,theta,K,D)
  base<-data.frame(theta=theta,R0=R0,rho=rho,K=K,status=ode$status,
    y_min_ratio=ode$y_min_ratio %||% NA_real_,x_in_ODE=ode$x_in,
    x_in_first=m$x_first,x_in_second=m$x_second,C=m$C,D=m$D)
  if(ode$status!='OK'||anyNA(c(m$x_first,m$x_second))) return(base)
  q0<-kendall_quantities(ode$x_in,R0,rho,theta,K)
  q1<-kendall_quantities(m$x_first,R0,rho,theta,K)
  q2<-kendall_quantities(m$x_second,R0,rho,theta,K)
  transform(base,
    x_relerr_first=(m$x_first-ode$x_in)/ode$x_in,
    x_relerr_second=(m$x_second-ode$x_in)/ode$x_in,
    Lambda_ODE=q0$Lambda,Lambda_first=q1$Lambda,Lambda_second=q2$Lambda,
    DeltaLambda_first=q1$Lambda-q0$Lambda,DeltaLambda_second=q2$Lambda-q0$Lambda,
    order_gap_Lambda=abs(q2$Lambda-q1$Lambda),B_ref=q0$B,B_K_first=q1$B,
    B_K_second=q2$B,B_L_ODE=q0$B_L,B_L_first=q1$B_L,B_L_second=q2$B_L,
    P1_ref=q0$P1,P1_K_first=q1$P1,P1_K_second=q2$P1,
    P1_L_ODE=q0$P1_L,P1_L_first=q1$P1_L,P1_L_second=q2$P1_L)
}
`%||%`<-function(x,y) if(is.null(x)) y else x

run_grid <- function(grid,file) {
  keys<-unique(grid[c('R0','theta')]); const<-vector('list',nrow(keys))
  for(i in seq_len(nrow(keys))) {
    d<-D_regularized(keys$R0[i],keys$theta[i]); const[[i]]<-d
    cat('constants',i,'/',nrow(keys),'R0=',keys$R0[i],'theta=',keys$theta[i],'D=',d,'\n')
  }
  ans<-vector('list',nrow(grid))
  for(i in seq_len(nrow(grid))) {
    k<-match(paste(grid$R0[i],grid$theta[i]),paste(keys$R0,keys$theta))
    ans[[i]]<-tryCatch(do.call(run_point,c(as.list(grid[i,]),list(D=const[[k]]))),
                       error=function(e) data.frame(grid[i,],status=paste0('ERROR:',conditionMessage(e))))
    if(i%%10==0) cat('points',i,'/',nrow(grid),'\n')
  }
  out<-data.table::rbindlist(ans,fill=TRUE); dir.create(dirname(file),FALSE,TRUE)
  data.table::fwrite(out,file); out
}
