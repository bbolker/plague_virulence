# Run only once the high-R scan is complete; preserve the pre-extension PDFs.
source('validation/R/theory.R')
stopifnot(requireNamespace('data.table',quietly=TRUE))
source('validation/R/validation_plot_data.R')
z<-load_validation_plot_data()
stopifnot(z$extended,nrow(z$sim)==460L)
backup<-'validation/figures/before_high_R_extension_20261102'
dir.create(backup,recursive=TRUE,showWarnings=FALSE)
for(name in c('main_persistence_validation.pdf','unconditional_persistence_validation.pdf')){
 src<-file.path('validation/figures/paper',name);dst<-file.path(backup,name)
 if(!file.exists(dst))stopifnot(file.copy(src,dst,overwrite=FALSE))
}
data.table::fwrite(z$sim,'validation/data/fig1_extended_display_results.csv')
rbin<-file.path(R.home('bin'),'Rscript')
for(script in c('validation/make_paper_figures.R','validation/make_unconditional_validation_figure.R')){
 args<-c(script,if(grepl('make_paper_figures',script))'--fig1-only')
 status<-system2(rbin,args)
 stopifnot(status==0)
}
writeLines(c('Both validation PDFs updated after all 180 high-R parameter combinations completed.',
 'Original 280 display rows, including the pooled replacement point, preserved.',
 '180 new combinations, 10000 attempts each, epsilon=.01, 8 workers.',
 'Conditional and unconditional probabilities use the same epidemic trajectories.',
 'Stable log-coordinate C/action calculations used for both theoretical curves.',
 paste('Completed:',Sys.time())), 'validation/data/fig1_high_R_figures_completed.txt')
cat('Both extended validation PDFs complete.\n')
