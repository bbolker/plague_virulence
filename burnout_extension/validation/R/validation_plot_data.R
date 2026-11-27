# Shared data/predictions for conditional and unconditional validation figures.
# High-R counts are included only after the entire fixed-design scan completes.
source('validation/R/extrema.R')

# Recompute 95% Clopper--Pearson intervals from retained simulation counts.
exact_binomial_intervals <- function(sim, conf.level=.95) {
  bounds <- function(k,n) {
    stopifnot(all(is.na(k) | is.na(n) | (k>=0 & k<=n)))
    lo <- hi <- rep(NA_real_,length(k))
    ok <- is.finite(k) & is.finite(n) & n>0
    alpha <- 1-conf.level
    lo[ok & k==0] <- 0
    hi[ok & k==n] <- 1
    j <- ok & k>0
    lo[j] <- stats::qbeta(alpha/2,k[j],n[j]-k[j]+1)
    j <- ok & k<n
    hi[j] <- stats::qbeta(1-alpha/2,k[j]+1,n[j]-k[j])
    list(low=lo,high=hi)
  }
  cond <- bounds(sim$persistent,sim$established)
  uncond <- bounds(sim$persistent,sim$attempts)
  sim[, `:=`(cond_low=cond$low,cond_high=cond$high,
              uncond_low=uncond$low,uncond_high=uncond$high)]
  sim
}

load_validation_plot_data <- function() {
  original <- 'validation/data/fig1_display_results.csv'
  if(!file.exists(original)) original <- 'validation/data/fig11_scan_results.csv'
  sim <- data.table::fread(original)
  extension <- 'validation/data/fig1_high_R_results.csv'
  corrected_file <- 'validation/data/fig1_high_R_guarded_results.csv'
  corrected <- FALSE
  if(file.exists(corrected_file)) {
    candidate <- data.table::fread(corrected_file,integer64="double")
    if(nrow(candidate)==180L && all(candidate$attempts==10000L) &&
       all(candidate$unresolved==0L)) {
      extension <- corrected_file
      corrected <- TRUE
    }
  }
  extended <- FALSE
  if(file.exists(extension)) {
    z <- data.table::fread(extension,integer64="double")
    if(nrow(z)==180L && all(z$attempts==10000L) && all(z$unresolved==0L)) {
      stopifnot(all(z$R0>max(sim$R0)),all(is.finite(z$P_conditional)),
                all(z$established>0),all(z$persistent<=z$established))
      z[,sample_source:=if(corrected)'boundary_guarded_high_R_fixed10000' else
                          'high_R_fixed10000_scan']
      sim <- data.table::rbindlist(list(sim,z),use.names=TRUE,fill=TRUE)
      extended <- TRUE
    }
  }
  sim <- exact_binomial_intervals(sim)
  list(sim=sim,extended=extended,corrected=corrected,
       source=extension,Rmax=if(extended)128 else 6)
}

stable_validation_probability <- function(R0,rho,theta,K,conditional=TRUE) {
  logC <- unname(C_log_derivatives(R0,theta)['logC'])
  A <- unname(action_barrier_derivatives(R0,theta)['value'])
  logB <- log(K)+logC+.5*(log(rho)+log(R0-1)-theta*log(R0)-log(2*pi))-A/rho
  p <- -expm1(-exp(pmin(logB,700)))
  if(conditional) p else (1-1/R0)*p
}
