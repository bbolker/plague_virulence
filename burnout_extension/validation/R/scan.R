source('validation/R/theory.R')
source('validation/R/ode_reference.R')
source('validation/R/kendall.R')

run_point <- function(R0,rho,theta,K,D=NULL,boundary_choice='sqrt') {
  ystar<-rho*h_theta(1/R0,theta)
  yline<-boundary_layer_height(ystar,K,boundary_choice)
  ode<-ode_reference(R0,rho,theta,K,yline=yline)
  m<-matching(R0,rho,theta,K,D,yline=yline)
  base<-data.frame(theta=theta,R0=R0,rho=rho,K=K,
    boundary_choice=boundary_choice,boundary=unname(boundary_layer_labels[boundary_choice]),status=ode$status,
    y_min_ratio=ode$y_min_ratio %||% NA_real_,x_in_ODE=ode$x_in,
    x_in_first=m$x_first,x_in_second=m$x_second,C=m$C,D=m$D)
  if(ode$status!='OK'||anyNA(c(m$x_first,m$x_second))) return(base)
  q0<-kendall_quantities(ode$x_in,R0,rho,theta,K,yline)
  q1<-kendall_quantities(m$x_first,R0,rho,theta,K,yline)
  q2<-kendall_quantities(m$x_second,R0,rho,theta,K,yline)
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

run_grid <- function(grid,file,
    boundary_choices=c('sqrt','compromise','compromise_3_4','ystar')) {
  keys<-unique(grid[c('R0','theta')]); const<-vector('list',nrow(keys))
  for(i in seq_len(nrow(keys))) {
    d<-D_regularized(keys$R0[i],keys$theta[i]); const[[i]]<-d
    cat('constants',i,'/',nrow(keys),'R0=',keys$R0[i],'theta=',keys$theta[i],'D=',d,'\n')
  }
  jobs<-expand.grid(grid_row=seq_len(nrow(grid)),boundary_choice=boundary_choices,
    KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
  worker <- function(i,jobs,grid,keys,const) {
    j<-jobs$grid_row[i]
    k<-match(paste(grid$R0[j],grid$theta[j]),paste(keys$R0,keys$theta))
    tryCatch(do.call(run_point,c(as.list(grid[j,]),
      list(D=const[[k]],boundary_choice=jobs$boundary_choice[i]))),
      error=function(e) data.frame(grid[j,],boundary_choice=jobs$boundary_choice[i],
        boundary=unname(boundary_layer_labels[jobs$boundary_choice[i]]),
        status=paste0('ERROR:',conditionMessage(e))))
  }
  cores <- min(8L,parallel::detectCores(logical=FALSE),nrow(jobs))
  cl <- parallel::makeCluster(cores); on.exit(parallel::stopCluster(cl),add=TRUE)
  parallel::clusterEvalQ(cl,source('validation/R/scan.R'))
  ans <- parallel::parLapplyLB(cl,seq_len(nrow(jobs)),worker,
    jobs=jobs,grid=grid,keys=keys,const=const)
  out<-data.table::rbindlist(ans,fill=TRUE); dir.create(dirname(file),FALSE,TRUE)
  data.table::fwrite(out,file); out
}
