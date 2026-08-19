source('validation/R/stochastic.R')

stopifnot(requireNamespace('data.table',quietly=TRUE))

checkpoint_file <- 'validation/data/stochastic_scan_checkpoint.rds'
results_file <- 'validation/data/stochastic_scan_results.csv'
scan_version <- 'exact_ctmc_v2_rho002_batch500_Dregularized'
batch_size <- 500L
base_attempts <- 3000L
expanded_attempts <- 10000L
ci_width_threshold <- .05
seed_base <- 1900000000L

R0_minus_1 <- c(.05,.075,.10,.15,.20,.30,.50,.75,1,1.5,2,3,5)
grid <- expand.grid(
  rho=c(.01,.02,.05,.10),theta=c(0,.5,1),K=c(1000,3000,10000,30000),
  R0=1+R0_minus_1,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
grid <- grid[order(grid$rho,grid$theta,grid$K,grid$R0),]
grid$point_id <- seq_len(nrow(grid))

new_state <- function() {
  counts <- transform(grid,attempts=0L,established=0L,persistent=0L,
                      unresolved=0L,target_attempts=base_attempts,
                      completed_batches=0L,elapsed_seconds=0)
  list(version=scan_version,created=as.character(Sys.time()),grid=grid,
       counts=counts)
}

checkpoint_candidates <- c(checkpoint_file,paste0(checkpoint_file,'.tmp'))
valid_checkpoint <- Filter(function(f) file.exists(f) &&
  !inherits(try(readRDS(f),silent=TRUE),'try-error'),checkpoint_candidates)
if(length(valid_checkpoint)) {
  newest <- valid_checkpoint[which.max(file.info(valid_checkpoint)$mtime)]
  state <- readRDS(newest)
  if(!identical(state$version,scan_version) || !identical(state$grid,grid)) {
    key <- function(z) do.call(paste,c(z[c('rho','theta','K','R0')],sep='|'))
    old_key <- key(state$grid); new_key <- key(grid)
    if(!all(old_key %in% new_key))
      stop('Existing checkpoint does not match this scan specification.')
    old_counts <- state$counts; migrated <- new_state()
    at <- match(old_key,new_key)
    copy_columns <- setdiff(intersect(names(old_counts),names(migrated$counts)),
                            c('rho','theta','K','R0','point_id'))
    for(nm in copy_columns) migrated$counts[at,nm] <- old_counts[[nm]]
    migrated$created <- state$created
    migrated$migrated <- as.character(Sys.time())
    state <- migrated
    cat('Migrated',length(at),'completed parameter points into expanded rho grid\n')
  }
  cat('Resuming checkpoint with',sum(state$counts$attempts),'completed trajectories\n')
} else state <- new_state()

atomic_save <- function(state) {
  tmp <- paste0(checkpoint_file,'.tmp')
  saveRDS(state,tmp,compress=FALSE)
  if(!file.copy(tmp,checkpoint_file,overwrite=TRUE))
    stop('Could not replace checkpoint')
  z <- data.table::as.data.table(state$counts)
  z[,P_unconditional:=persistent/attempts]
  z[,`:=`(uncond_low=NA_real_,uncond_high=NA_real_,
          P_conditional=NA_real_,cond_low=NA_real_,cond_high=NA_real_)]
  for(j in which(z$attempts>0L)) {
    u <- wilson(z$persistent[j],z$attempts[j])
    z[j,`:=`(uncond_low=u[1],uncond_high=u[2])]
    if(z$established[j]>0L) {
      cc <- wilson(z$persistent[j],z$established[j])
      z[j,`:=`(P_conditional=persistent/established,
               cond_low=cc[1],cond_high=cc[2])]
    }
  }
  data.table::fwrite(z,results_file)
  file.remove(tmp)
}

cl <- parallel::makeCluster(8L)
on.exit(parallel::stopCluster(cl),add=TRUE)
parallel::clusterExport(cl,'one_ctmc',envir=environment())

for(i in seq_len(nrow(state$counts))) {
  repeat {
    row <- state$counts[i,]
    if(row$attempts>=row$target_attempts) {
      if(row$attempts==base_attempts) {
        wu <- diff(wilson(row$persistent,row$attempts))
        wc <- if(row$established>0L)
          diff(wilson(row$persistent,row$established)) else Inf
        if(max(wu,wc)>ci_width_threshold) {
          state$counts$target_attempts[i] <- expanded_attempts
          atomic_save(state)
          cat(sprintf('point %d/%d expands to %d (CI widths %.4f, %.4f)\n',
            i,nrow(state$counts),expanded_attempts,wu,wc))
          next
        }
      }
      break
    }
    n_now <- min(batch_size,row$target_attempts-row$attempts)
    batch_id <- row$completed_batches+1L
    seed <- seed_base+i*100L+batch_id
    t0 <- proc.time()[['elapsed']]
    got <- simulate_counts_batch(cl,row$R0,row$rho,row$theta,row$K,
                                 n_now,seed)
    elapsed <- proc.time()[['elapsed']]-t0
    state$counts$attempts[i] <- row$attempts+got['attempts']
    state$counts$established[i] <- row$established+got['established']
    state$counts$persistent[i] <- row$persistent+got['persistent']
    state$counts$unresolved[i] <- row$unresolved+got['unresolved']
    state$counts$completed_batches[i] <- batch_id
    state$counts$elapsed_seconds[i] <- row$elapsed_seconds+elapsed
    state$updated <- as.character(Sys.time())
    atomic_save(state)
    cat(sprintf('point %d/%d batch %d: %d/%d attempts (%.1f s)\n',
      i,nrow(state$counts),batch_id,state$counts$attempts[i],
      state$counts$target_attempts[i],elapsed))
    if(got['unresolved']>0L)
      stop('Unresolved trajectories at point ',i,', batch ',batch_id)
  }
}

state$completed <- as.character(Sys.time())
atomic_save(state)
cat('stochastic scan complete\n')

# Retain a small reproducible set of full event histories for diagnostics; the
# complete scan itself is stored as aggregate counts because full histories for
# millions of trajectories would be unnecessarily large.
diagnostic_file <- 'validation/data/stochastic_scan_diagnostic_trajectories.rds'
if(!file.exists(diagnostic_file)) {
  diagnostic_points <- expand.grid(R0=2,rho=c(.01,.02,.05,.10),theta=c(0,.5,1),
    K=10000,KEEP.OUT.ATTRS=FALSE)
  diagnostics <- vector('list',nrow(diagnostic_points))
  for(i in seq_len(nrow(diagnostic_points))) {
    set.seed(2026082900+i)
    traces <- vector('list',5L);summaries <- vector('list',5L)
    for(j in seq_len(5L)) {
      z <- do.call(one_ctmc_trace,as.list(diagnostic_points[i,]))
      traces[[j]] <- z$trajectory
      z$summary$trajectory_id <- j;summaries[[j]] <- z$summary
    }
    diagnostics[[i]] <- list(parameters=diagnostic_points[i,],
      summaries=data.table::rbindlist(summaries),trajectories=traces)
  }
  saveRDS(diagnostics,diagnostic_file,compress='xz')
}
