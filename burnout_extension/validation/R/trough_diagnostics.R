source('validation/R/theory.R')
source('validation/R/ode_reference.R')

trough_phaseplane_approx <- function(R0,rho,theta,order=c('first','second'),D=NULL) {
  order<-match.arg(order); xf<-x_final(R0); xs<-1/R0
  delta<-1-R0*xf; hf<-h_theta(xf,theta); a<-hf/delta; b<-R0*xf/delta
  cc<-delta*h_prime(xf,theta)+R0*hf; C<-C_explicit(R0,theta)
  if(order=='second' && is.null(D)) D<-D_regularized(R0,theta)
  Q<-(R0*a-b*cc)/(a*delta^2); target<-action_diff(xf,xs,R0,theta)
  balance<-function(logy) {
    y<-exp(logy); L<-log(C/y); M<-rho*L+(b/a)*y
    if(order=='second') M<-M+rho^2*D/a+rho*y*Q*(L+1)+(b*Q/(2*a))*y^2
    M-target
  }
  # The low-y root is the trough-side solution of the corner balance.  Search
  # in log(y) so exponentially shallow troughs remain numerically resolvable.
  upper<-log(min(C*.999,1)); mesh<-seq(-740,upper,length.out=3000)
  value<-vapply(mesh,balance,0.); cross<-which(value[-length(value)]*value[-1]<=0)
  if(!length(cross)) return(list(status='NO_ASYMPTOTIC_ROOT',y_min=NA_real_))
  j<-cross[1L]
  root<-uniroot(balance,c(mesh[j],mesh[j+1L]),tol=1e-11)$root
  list(status='OK',y_min=exp(root))
}

boundary_entry_from_trajectory <- function(trajectory,y_boundary,peak_index) {
  y<-trajectory$y; x<-trajectory$x; n<-nrow(trajectory)
  ii<-seq.int(peak_index,n-1L)
  entry<-ii[y[ii]>=y_boundary & y[ii+1L]<y_boundary]
  if(!length(entry)) return(NA_real_)
  i<-entry[1L]; w<-(y_boundary-y[i])/(y[i+1L]-y[i])
  x[i]+w*(x[i+1L]-x[i])
}

critical_K_estimate <- function(y_min,y_star,boundary_choice,log_constant=1) {
  if(!is.finite(y_min)||y_min<=0) return(NA_real_)
  ans<-switch(boundary_choice,
    sqrt=y_star/y_min^2,
    compromise=y_star^2/y_min^3,
    compromise_3_4=y_star^3/y_min^4,
    ystar=if(y_min<=y_star) Inf else NA_real_,
    log_inverse=exp(pmin(709,log_constant/y_min)),
    ystar_log_inverse=exp(pmin(709,log_constant*y_star/y_min)),NA_real_)
  ans
}

# Diagnose boundary entry only; this does not alter the burnout probability.
run_trough_diagnostics <- function(grid,file,log_constant=1,
    boundary_choices=c('sqrt','compromise','compromise_3_4','ystar',
                       'log_inverse','ystar_log_inverse'),dt=.1) {
  required<-c('R0','rho','theta','K')
  if(!all(required %in% names(grid))) stop('grid must contain: ',paste(required,collapse=', '))
  dyn<-unique(grid[c('R0','rho','theta')]); trajectories<-vector('list',nrow(dyn))
  for(i in seq_len(nrow(dyn))) {
    p<-dyn[i,]; tr<-do.call(first_deterministic_trough,
      c(as.list(p),list(K=Inf,init='manifold',dt=dt,trajectory=TRUE)))
    D<-tryCatch(D_regularized(p$R0,p$theta),error=function(e) NA_real_)
    ap1<-trough_phaseplane_approx(p$R0,p$rho,p$theta,'first')
    ap2<-if(is.finite(D)) trough_phaseplane_approx(p$R0,p$rho,p$theta,'second',D) else
      list(status='NO_CONSTANT',y_min=NA_real_)
    trajectories[[i]]<-list(trough=tr,D=D,approx1=ap1,approx2=ap2)
    cat('trajectory',i,'/',nrow(dyn),'R0=',p$R0,'rho=',p$rho,'theta=',p$theta,
        'status=',tr$status,'\n')
  }
  rows<-vector('list',nrow(grid)*length(boundary_choices)); out_i<-0L
  for(j in seq_len(nrow(grid))) {
    k<-match(paste(grid$R0[j],grid$rho[j],grid$theta[j]),
             paste(dyn$R0,dyn$rho,dyn$theta)); q<-trajectories[[k]]; tr<-q$trough
    z<-tr$trajectory; xs<-1/grid$R0[j]
    peak<-which(z$x[-nrow(z)]>=xs & z$x[-1]<xs)[1]
    y_peak<-if(is.na(peak)) NA_real_ else max(z$y[seq_len(peak+1L)])
    ystar<-grid$rho[j]*h_theta(xs,grid$theta[j])
    for(choice in boundary_choices) {
      out_i<-out_i+1L
      yb<-boundary_layer_height(ystar,grid$K[j],choice,log_constant)
      ratio<-tr$y_min/yb
      xin<-if(is.na(peak)) NA_real_ else boundary_entry_from_trajectory(z,yb,peak+1L)
      rows[[out_i]]<-data.frame(R0=grid$R0[j],theta=grid$theta[j],rho=grid$rho[j],
        K=grid$K[j],y_star=ystar,boundary_type=choice,
        boundary_label=unname(boundary_layer_labels[choice]),log_constant=log_constant,
        y_boundary=yb,x_trough=tr$x_min,y_trough=tr$y_min,t_trough=tr$t_min,
        log_y_trough=tr$log_y_min,
        K_y_trough=grid$K[j]*tr$y_min,y_trough_over_boundary=ratio,
        log10_trough_over_boundary=log10(ratio),
        trough_below_boundary=is.finite(ratio)&&ratio<=1,
        hits_boundary=is.finite(xin),x_boundary_entry=xin,y_peak=y_peak,
        trough_status=tr$status,nullcline_residual=tr$nullcline_residual,
        y_trough_phaseplane_first=q$approx1$y_min,
        phaseplane_first_status=q$approx1$status,
        y_trough_phaseplane_second=q$approx2$y_min,
        phaseplane_second_status=q$approx2$status,
        phaseplane_first_rel_error=(q$approx1$y_min-tr$y_min)/tr$y_min,
        phaseplane_second_rel_error=(q$approx2$y_min-tr$y_min)/tr$y_min,
        K_critical_estimate=critical_K_estimate(tr$y_min,ystar,choice,log_constant))
    }
  }
  out<-data.table::rbindlist(rows,fill=TRUE)
  dir.create(dirname(file),FALSE,TRUE);data.table::fwrite(out,file);out
}
