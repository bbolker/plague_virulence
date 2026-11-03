source('validation/R/stochastic.R')
stopifnot(requireNamespace('data.table',quietly=TRUE),
          requireNamespace('adaptivetau',quietly=TRUE))

checkpoint_file <- 'validation/data/fig11_scan_checkpoint.rds'
results_file <- 'validation/data/fig11_scan_results.csv'
scan_version <- 'theta_common_grid_rho001_adaptivetau_eps001_v7_fixed10000'
classification_version <- 'first_trough_turning_chronological_v1'
batch_size <- 500L
fixed_attempts <- 10000L
seed_base <- 2100000000L

R0_minus_1 <- c(.03,.05,.075,.10,.15,.20,.30,.50,.75,1,1.5,2,3,5)
grid <- expand.grid(rho=.01,theta=c(0,.25,.5,.75,1),K=c(1e6,1e7,1e8,1e9),
  R0=1+R0_minus_1,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
dir.create(dirname(checkpoint_file),recursive=TRUE,showWarnings=FALSE)
row.names(grid) <- NULL
grid <- grid[order(grid$theta,grid$K,grid$R0),]
grid$point_id <- seq_len(nrow(grid))

new_state <- function() list(version=scan_version,classification_version=classification_version,
  sampling_policy='fixed_attempts_per_point',fixed_attempts=fixed_attempts,
  seed_base=seed_base,batch_size=batch_size,workers=8L,
  created=as.character(Sys.time()),
  grid=grid,counts=transform(grid,method='adaptive tau',tau_epsilon=.01,
    attempts=0L,established=0L,persistent=0L,unresolved=0L,
    target_attempts=fixed_attempts,completed_batches=0L,elapsed_seconds=0))

if(file.exists(checkpoint_file)) {
  state <- readRDS(checkpoint_file)
  if(!identical(state$classification_version,classification_version))
    stop('Checkpoint uses an older classification rule. Archive it and start a fresh scan; old counts cannot be reclassified without trajectories.')
  if(!identical(state$grid,grid))
    stop('Checkpoint grid differs from this scan. Preserve the original point ordering to continue the same random-seed sequences.')
  if(!identical(state$version,scan_version)) {
    old <- state$counts; migrated <- new_state()
    key <- function(z) do.call(paste,c(z[c('rho','theta','K','R0')],sep='|'))
    at <- match(key(old),key(migrated$counts))
    keep <- which(!is.na(at))
    copy_columns <- setdiff(intersect(names(old),names(migrated$counts)),
      c('rho','theta','K','R0','point_id','target_attempts'))
    for(nm in copy_columns) migrated$counts[at[keep],nm] <- old[[nm]][keep]
    migrated$created <- state$created
    migrated$migrated <- as.character(Sys.time())
    state <- migrated
    atomic_migration <- paste0(checkpoint_file,'.migration.tmp')
    saveRDS(state,atomic_migration,compress=FALSE)
    if(!file.copy(atomic_migration,checkpoint_file,overwrite=TRUE))
      stop('Could not save migrated checkpoint.')
    file.remove(atomic_migration)
    cat('Migrated',length(keep),'overlapping points into common grid\n')
  }
  cat('Resuming with',sum(state$counts$attempts),'saved trajectories\n')
} else state <- new_state()

stopifnot(all(state$counts$attempts<=fixed_attempts),
          all(state$counts$attempts==batch_size*state$counts$completed_batches))
state$counts$target_attempts[] <- fixed_attempts
if(any(state$counts$attempts<fixed_attempts)) state$completed <- NULL

atomic_save <- function(state) {
  tmp <- paste0(checkpoint_file,'.tmp')
  saveRDS(state,tmp,compress=FALSE)
  if(!file.copy(tmp,checkpoint_file,overwrite=TRUE)) stop('Checkpoint replace failed')
  z <- data.table::as.data.table(state$counts)
  z[,`:=`(P_unconditional=persistent/attempts,
          P_conditional=persistent/established,
          uncond_low=NA_real_,uncond_high=NA_real_,
          cond_low=NA_real_,cond_high=NA_real_)]
  for(j in which(z$attempts>0L)) {
    u <- wilson(z$persistent[j],z$attempts[j])
    z[j,`:=`(uncond_low=u[1],uncond_high=u[2])]
    if(z$established[j]>0L) {
      cc <- wilson(z$persistent[j],z$established[j])
      z[j,`:=`(cond_low=cc[1],cond_high=cc[2])]
    }
  }
  data.table::fwrite(z,results_file)
  file.remove(tmp)
}

atomic_save(state)

cl <- parallel::makeCluster(8L)
on.exit(parallel::stopCluster(cl),add=TRUE)
parallel::clusterExport(cl,c('one_adaptive_tau','classify_first_trough'),envir=environment())

for(i in seq_len(nrow(state$counts))) repeat {
  row <- state$counts[i,]
  if(row$attempts>=row$target_attempts) break
  n_now <- min(batch_size,row$target_attempts-row$attempts)
  batch_id <- row$completed_batches+1L
  seed <- seed_base+i*100L+batch_id
  t0 <- proc.time()[['elapsed']]
  got <- simulate_counts_batch_tau(cl,row$R0,row$rho,row$theta,row$K,
    n_now,seed,tau_epsilon=.01)
  state$counts$attempts[i] <- row$attempts+got['attempts']
  state$counts$established[i] <- row$established+got['established']
  state$counts$persistent[i] <- row$persistent+got['persistent']
  state$counts$unresolved[i] <- row$unresolved+got['unresolved']
  state$counts$completed_batches[i] <- batch_id
  state$counts$elapsed_seconds[i] <- row$elapsed_seconds+
    proc.time()[['elapsed']]-t0
  state$updated <- as.character(Sys.time())
  atomic_save(state)
  cat(sprintf('point %d/%d: %d/%d attempts\n',i,nrow(state$counts),
    state$counts$attempts[i],state$counts$target_attempts[i]))
  if(got['unresolved']>0L) stop('Unresolved tau trajectories at point ',i)
}

stopifnot(all(state$counts$attempts==fixed_attempts),
          all(state$counts$unresolved==0L))
state$completed <- as.character(Sys.time())
atomic_save(state)
cat('Figure 1 fixed-size stochastic scan complete\n')
