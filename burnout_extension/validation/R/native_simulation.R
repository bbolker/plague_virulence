compile_native_simulation <- function(source_file) {
  work<-tempfile('burnout_native_');dir.create(work)
  src<-file.path(work,basename(source_file));stopifnot(file.copy(source_file,src))
  sdk<-'/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk'
  log<-system2(file.path(R.home('bin'),'R'),c('CMD','SHLIB',shQuote(src)),
    env=if(dir.exists(sdk))paste0('SDKROOT=',sdk) else character(),stdout=TRUE,stderr=TRUE)
  dll<-sub('\\.c$','.so',src)
  if(!file.exists(dll))stop(paste(log,collapse='\n'))
  dll
}

native_batch <- function(cl,row,n,seed) {
  parallel::clusterSetRNGStream(cl,iseed=seed)
  chunks<-rep(n%/%length(cl),length(cl))
  if(n%%length(cl))chunks[seq_len(n%%length(cl))]<-chunks[seq_len(n%%length(cl))]+1L
  a<-parallel::parLapply(cl,chunks,function(n,row){
    if(row$method=='exact') {
      z<-.C('exact_counts',as.double(c(row$R0,row$rho,row$theta,row$K)),as.integer(n),out=double(6))$out
      z<-c(z,0,0,0)
    }else{
      z<-.C('guarded_counts',as.double(c(row$R0,row$rho,row$theta,row$K,
        row$epsilon,row$event_budget,row$small_count,row$tf)),as.integer(n),out=double(9))$out
    }
    c(attempts=n,established=z[1],persistent=z[2],unresolved=z[3],
      ever_zero=z[4],trough=z[5],events=z[6],tau_steps=z[7],exact_steps=z[8],rejected=z[9])
  },row)
  colSums(do.call(rbind,a))
}
