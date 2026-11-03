# Fixed-design convergence audit. No plot-data changes in this script.
source('validation/R/stochastic.R')
outdir<-'validation/data/tau_bias_audit_20261104'
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
work<-tempfile('exact_audit_');dir.create(work)
src<-file.path(work,'exact.c')
stopifnot(file.copy('validation/src/exact_first_trough.c',src))
sdk<-'/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk'
log<-system2(file.path(R.home('bin'),'R'),c('CMD','SHLIB',shQuote(src)),
 env=if(dir.exists(sdk))paste0('SDKROOT=',sdk) else character(),stdout=TRUE,stderr=TRUE)
dll<-file.path(work,'exact.so');if(!file.exists(dll))stop(paste(log,collapse='\n'))
cl<-parallel::makeCluster(8L)
parallel::clusterCall(cl,function(p)dyn.load(p),dll)

one_tau_audit<-function(R0,rho,theta,K,epsilon,tf=100){
 trans<-list(c(S=-1,I=1),c(I=-1),c(S=1))
 rates<-function(z,p,t){S<-z['S'];I<-z['I'];c(p$R0*S*I/p$K,I,
   if(S<p$K)p$K*p$rho*(1-S/p$K)*(S/p$K)^p$theta else 0)}
 z<-adaptivetau::ssa.adaptivetau(c(S=K-1,I=1),trans,rates,
    list(R0=R0,rho=rho,theta=theta,K=K),tf=tf,
    tl.params=list(epsilon=epsilon,extraChecks=TRUE))
 a<-classify_first_trough(z[,'S'],z[,'I'],R0,rho,theta,K)
 c(established=as.integer(a['established']),persistent=as.integer(isTRUE(a['persist'])),
   unresolved=as.integer(is.na(a['persist'])),ever_zero=as.integer(any(z[,'S']==0)))
}
parallel::clusterExport(cl,c('one_tau_audit','classify_first_trough'))
design<-data.frame(method=c('exact',rep('tau',5)),epsilon=c(NA,.03,.01,.003,.001,.0003),
 attempts=c(10000,rep(4000,5)),seed=202611400L+1:6,R0=128,rho=.01,theta=.25,K=1e6,tf=100)
write.csv(design,file.path(outdir,'convergence_design.csv'),row.names=FALSE)
results<-list()
for(i in seq_len(nrow(design))){
 row<-design[i,];totals<-c(attempts=0,established=0,persistent=0,unresolved=0,ever_zero=0,events=0)
 start<-proc.time()[3]
 for(batch in seq_len(row$attempts/500)){
  parallel::clusterSetRNGStream(cl,iseed=row$seed+batch*100L)
  got<-parallel::parLapply(cl,c(63L,63L,63L,63L,62L,62L,62L,62L),function(n,row){
   if(row$method=='exact'){
    a<-.C('exact_counts',as.double(c(row$R0,row$rho,row$theta,row$K)),as.integer(n),out=double(6))$out
    return(c(attempts=n,established=a[1],persistent=a[2],unresolved=a[3],ever_zero=a[4],events=a[6]))
   }
   a<-replicate(n,one_tau_audit(row$R0,row$rho,row$theta,row$K,row$epsilon,row$tf))
   z<-rowSums(a);c(attempts=n,established=unname(z[1]),persistent=unname(z[2]),
     unresolved=unname(z[3]),ever_zero=unname(z[4]),events=NA)
  },row)
  totals<-totals+colSums(do.call(rbind,got))
  print(c(method=row$method,epsilon=row$epsilon,attempts=totals['attempts']),quote=FALSE)
  saveRDS(list(design=design,finished=results,current=list(row=row,totals=totals)),file.path(outdir,'convergence_checkpoint.rds'))
 }
 ci<-wilson(totals['persistent'],totals['established'])
 results[[i]]<-cbind(row,data.frame(established=totals['established'],persistent=totals['persistent'],
  unresolved=totals['unresolved'],ever_zero=totals['ever_zero'],
  P=totals['persistent']/totals['established'],CI_low=max(0,ci[1]),CI_high=min(1,ci[2]),
  seconds=proc.time()[3]-start))
 write.csv(do.call(rbind,results),file.path(outdir,'convergence_results.csv'),row.names=FALSE)
 print(results[[i]],row.names=FALSE)
 stopifnot(totals['unresolved']==0)
}
parallel::stopCluster(cl)
