source('validation/R/trough_diagnostics.R')
full <- '--full' %in% commandArgs(trailingOnly=TRUE)
if(full) {
  base<-data.table::fread('validation/data/K_R0_scan.csv',select=c('R0','rho','theta','K'))
  grid<-unique(as.data.frame(base))
  file<-'validation/data/first_trough_diagnostics.csv'
} else {
  grid<-expand.grid(R0=c(1.2,2,3,5),rho=c(.005,.02,.05),theta=c(0,.5,1),
                    K=c(1e3,1e4,1e5))
  file<-'validation/data/first_trough_diagnostics_small.csv'
}
run_trough_diagnostics(grid,file)
cat('wrote',file,'\n')
