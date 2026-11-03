# Small, fixed-design high-R validation. No manuscript or figure changes.
# Exact embedded jump chain: waiting times are unnecessary for these outcomes.
source('validation/R/theory.R')
source('validation/R/extrema.R')
source('validation/R/stochastic.R')
large_K <- '--large-K' %in% commandArgs(trailingOnly=TRUE)
outdir <- if(large_K) 'validation/data/high_R_large_K_20261102' else
 'validation/data/high_R_spotcheck_20261102'
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
work <- tempfile('burnout_exact_');dir.create(work)
cfile <- file.path(work,'exact.c')
writeLines('
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>
void exact_counts(double *p,int *n,double *out) {
  double R=p[0],rho=p[1],theta=p[2],K=p[3],xs=1/R;
  GetRNGstate();
  for(int k=0;k<*n;k++) {
    double S=K-1,I=1;
    int peak=0,trough=0,positive=0,turned=0,zero=0,done=0;
    for(long ev=0;ev<20000000;ev++) {
      if(ev%100000==0) R_CheckUserInterrupt();
      double old=S,inf=R*S*I/K,rem=I;
      double rec=S<K?K*rho*(1-S/K)*pow(S/K,theta):0;
      double total=inf+rem+rec;
      if(total<=0) break;
      double u=unif_rand()*total;
      if(u<inf){S--;I++;}else if(u<inf+rem){I--;}else{S++;}
      if(S==0)zero=1;
      if(!peak && old/K>xs && S/K<=xs && I>0)peak=1;
      if(peak && !trough && old/K<xs && S/K>=xs && I>0)trough=1;
      if(trough){
        double g=K*rho*(1-S/K)*pow(S/K,theta)-R*S*I/K;
        if(g>0)positive=1;
        if(positive && g<0)turned=1;
      }
      out[5]++;
      if(I==0){done=1;break;}
      if(trough && turned && old/K>xs && S/K<=xs && I>0){out[1]++;done=1;break;}
    }
    out[0]+=peak;out[2]+=!done;out[3]+=zero;out[4]+=trough;
  }
  PutRNGstate();
}',cfile)
rbin <- file.path(R.home('bin'),'R')
sdk <- '/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk'
compile_env <- if(dir.exists(sdk)) paste0('SDKROOT=',sdk) else character()
log <- system2(rbin,c('CMD','SHLIB',shQuote(cfile)),env=compile_env,stdout=TRUE,stderr=TRUE)
if(!file.exists(file.path(work,'exact.so')))stop(paste(log,collapse='\n'))
dyn.load(file.path(work,'exact.so'))
counts <- function(R0,rho,theta,K,n,seed){
 set.seed(seed)
 .C('exact_counts',as.double(c(R0,rho,theta,K)),as.integer(n),out=double(6))$out
}
# Exact reproducibility check against the existing R event-by-event simulator.
for(z in list(c(32,.01,0,1000),c(64,.01,.25,10000))){
 a<-counts(z[1],z[2],z[3],z[4],25,10203)
 set.seed(10203)
 b<-replicate(25,one_ctmc(z[1],z[2],z[3],z[4]))
 stopifnot(a[1]==sum(b['established',]),a[2]==sum(b['persist',]),a[3]==sum(is.na(b['persist',])))
}
grid<-if(large_K) data.frame(R0=rep(c(128,32,64,128),each=2),rho=.01,
 theta=rep(c(0,.25,.25,.25),each=2),K=rep(c(1e5,1e6),4)) else
 data.frame(R0=c(16,32,128,32,64,128),rho=.01,
 theta=c(0,0,0,.25,.25,.25),K=c(5000,1000,1000,60000,10000,4000))
grid$n_attempts<-2000L
grid$seed<-(if(large_K) 202611040L else 202611020L)+seq_len(nrow(grid))
write.csv(grid,file.path(outdir,'design.csv'),row.names=FALSE)
results<-list()
for(i in seq_len(nrow(grid))){
 z<-grid[i,];start<-proc.time()[3]
 a<-counts(z$R0,z$rho,z$theta,z$K,z$n_attempts,z$seed)
 elapsed<-proc.time()[3]-start
 stopifnot(a[3]==0)
 theory<-extremum_quantities(z$R0,z$rho,z$theta,z$K)
 pc<- -expm1(-exp(unname(theory['logB'])))
 ci<-wilson(a[2],a[1]);ui<-wilson(a[2],z$n_attempts)
 results[[i]]<-cbind(z,data.frame(established=a[1],persistent=a[2],unresolved=a[3],
  ever_zero_susceptibles=a[4],crossed_trough=a[5],events=a[6],seconds=elapsed,
  BI_conditional=pc,sim_conditional=a[2]/a[1],CI_low=ci[1],CI_high=ci[2],
  BI_unconditional=(1-1/z$R0)*pc,sim_unconditional=a[2]/z$n_attempts,
  uncond_CI_low=ui[1],uncond_CI_high=ui[2]))
 write.csv(do.call(rbind,results),file.path(outdir,'results.csv'),row.names=FALSE)
 print(results[[i]],row.names=FALSE);flush.console()
}
writeLines(c('Exact Gillespie embedded jump chain; no tau approximation.',
 'C outcome counts checked for exact same-seed agreement with existing R one_ctmc().',
 paste('All',nrow(grid),'parameter combinations and 2000 attempts each fixed before simulation.'),
 paste('Seeds',min(grid$seed),'through',max(grid$seed),'; no sample selection or seed selection.'),
 'Conditional denominator = established epidemics; unconditional denominator = all attempts.',
 'Intervals = 95% Wilson intervals; unresolved paths trigger an error.',
 capture.output(sessionInfo())),file.path(outdir,'provenance.txt'))
