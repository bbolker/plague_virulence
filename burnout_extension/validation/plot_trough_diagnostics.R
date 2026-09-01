library(data.table);library(ggplot2)
theme_set(theme_bw(base_size=9)+theme(panel.grid.minor=element_blank(),legend.position='bottom'))

plot_trough_diagnostics <- function(file='validation/data/first_trough_diagnostics.csv',
    outdir='validation/figures/trough_diagnostics') {
  d<-fread(file);dir.create(outdir,FALSE,TRUE)
  d[,boundary_label:=factor(boundary_label,levels=unique(boundary_label))]
  available_R0<-sort(unique(d$R0));target_R0<-c(1.2,2,3,5)
  representative_R0<-unique(vapply(target_R0,function(v)
    available_R0[which.min(abs(available_R0-v))],0.))
  representative<-d[R0 %in% representative_R0 & theta %in% c(0,.5,1)]
  p1<-ggplot(representative,aes(K,y_boundary,colour=boundary_label))+geom_line()+
    geom_line(aes(y=y_trough,group=interaction(R0,rho,theta)),colour='black',linewidth=.8)+
    scale_x_log10()+scale_y_log10()+facet_grid(theta+rho~R0,labeller=label_both)+
    labs(title='First deterministic trough and candidate boundaries',
      subtitle='Black: numerical trough (K-independent for manifold initialization)',
      x='K',y='infective fraction',colour='Boundary')
  ggsave(file.path(outdir,'A_trough_depth_vs_K.pdf'),p1,width=12,height=9,device=cairo_pdf)
  p2<-ggplot(representative,aes(K,y_trough_over_boundary,colour=boundary_label))+geom_line()+
    geom_hline(yintercept=1,linetype=2)+scale_x_log10()+scale_y_log10()+
    facet_grid(theta+rho~R0,labeller=label_both)+labs(
      title='First-trough depth relative to candidate boundary',x='K',
      y=expression(y[min]/y[boundary]),colour='Boundary')
  ggsave(file.path(outdir,'B_trough_boundary_ratio.pdf'),p2,width=12,height=9,device=cairo_pdf)
  p3<-ggplot(d,aes(log10(R0-1),factor(K),fill=hits_boundary))+geom_tile()+
    facet_grid(boundary_label+theta~rho,labeller=label_both)+
    scale_fill_manual(values=c('TRUE'='#2166AC','FALSE'='#B2182B'))+
    labs(title='Deterministic downward crossing of each boundary',x=expression(log[10](R[0]-1)),
      y='K',fill='Crosses')
  ggsave(file.path(outdir,'C_boundary_entry_heatmap.pdf'),p3,width=12,height=14,device=cairo_pdf)
  approx<-unique(d[,.(R0,rho,theta,phaseplane_first_rel_error,phaseplane_second_rel_error)])
  a<-melt(approx,id.vars=c('R0','rho','theta'),measure.vars=c(
    'phaseplane_first_rel_error','phaseplane_second_rel_error'),
    variable.name='order',value.name='relative_error')
  p4<-ggplot(a,aes(R0,relative_error,colour=order))+geom_hline(yintercept=0,linetype=2)+
    geom_line()+geom_point(size=1)+facet_grid(theta~rho,labeller=label_both)+
    scale_x_log10()+labs(title='Phase-plane trough approximation error',
      x=expression(R[0]),y='(approx - numeric) / numeric',colour='Approximation')
  ggsave(file.path(outdir,'D_phaseplane_trough_relative_error.pdf'),p4,width=10,height=7,device=cairo_pdf)
  invisible(c(p1=p1,p2=p2,p3=p3,p4=p4))
}

args<-commandArgs(trailingOnly=TRUE)
if(sys.nframe()==0L) plot_trough_diagnostics(if(length(args)) args[1] else
  'validation/data/first_trough_diagnostics.csv')
