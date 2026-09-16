library(data.table)
library(ggplot2)
library(grid)

peak<-fread('2d_ode_validation/data/psi1_peak_handoff_isolation.csv')
trough<-fread('2d_ode_validation/data/psi1_handoff_isolation.csv')
scan<-fread('2d_ode_validation/data/psi1_largeK_tau_results.csv')
orig<-scan[rho==.01&theta==.5&K==1e6&R0==4]

pdat<-rbindlist(list(
  data.table(label='Original CTMC\nconditioned on establishment',kind='Full CTMC',
    estimate=orig$P_conditional,low=orig$cond_low,high=orig$cond_high),
  data.table(label='Full CTMC from\ndeterministic first peak',kind='Full CTMC',
    estimate=peak[1]$estimate,low=peak[1]$low,high=peak[1]$high),
  data.table(label='Kendall from\ndeterministic first peak',kind='Kendall',
    estimate=peak[2]$estimate,low=NA_real_,high=NA_real_)))
tdat<-rbindlist(list(
  data.table(label='Full CTMC from\ndeterministic trough',kind='Full CTMC',
    estimate=trough[1]$estimate,low=trough[1]$low,high=trough[1]$high),
  data.table(label='Frozen-background\nbranching simulation',kind='Frozen branching',
    estimate=trough[2]$estimate,low=trough[2]$low,high=trough[2]$high),
  data.table(label='Kendall on\nfrozen background',kind='Kendall',
    estimate=trough[3]$estimate,low=NA_real_,high=NA_real_)))

cols<-c('Full CTMC'='#2C7FB8','Frozen branching'='#E69F00','Kendall'='#009E73')
base_theme<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position='none',
  plot.title=element_text(size=14),plot.subtitle=element_text(size=10),axis.title.y=element_blank())

p1<-ggplot(pdat,aes(estimate,factor(label,levels=rev(label)),colour=kind))+
  geom_errorbarh(aes(xmin=low,xmax=high),height=.12,linewidth=.7,na.rm=TRUE)+
  geom_point(size=3)+geom_label(aes(label=sprintf('%.3f',estimate)),nudge_y=.13,
    size=3.3,label.size=0,fill='white',show.legend=FALSE)+
  scale_colour_manual(values=cols)+scale_x_continuous(limits=c(0,1.04),breaks=seq(0,1,.2))+
  labs(title='A. Start at the first peak',
    subtitle='The full CTMC still matches the original 31% result, not Kendall',
    x='Probability alive at the deterministic recovery-peak horizon')+base_theme

p2<-ggplot(tdat,aes(estimate,factor(label,levels=rev(label)),colour=kind))+
  geom_errorbarh(aes(xmin=low,xmax=high),height=.12,linewidth=.7,na.rm=TRUE)+
  geom_point(size=3)+geom_label(aes(label=sprintf('%.4f',estimate)),nudge_y=.13,
    size=3.3,label.size=0,fill='white',show.legend=FALSE)+
  scale_colour_manual(values=cols)+
  scale_x_continuous(limits=c(.9845,1.0005),breaks=c(.985,.99,.995,1),labels=c('.985','.990','.995','1.000'))+
  labs(title='B. Start at the deterministic trough (zoomed)',
    subtitle='Frozen branching agrees with Kendall; full CTMC is only 1.1 points lower',
    x='Probability alive at the deterministic recovery-peak horizon')+base_theme

draw_all<-function(){
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2,2,widths=unit(c(1,1),'null'),
    heights=unit(c(.12,.88),'null'))))
  grid.text('Where does the psi = 1 persistence error arise?',
    vp=viewport(layout.pos.row=1,layout.pos.col=1:2),gp=gpar(fontsize=19,fontface='bold'))
  print(p1,vp=viewport(layout.pos.row=2,layout.pos.col=1))
  print(p2,vp=viewport(layout.pos.row=2,layout.pos.col=2))
  upViewport()
}

dir.create('2d_ode_validation/figures',FALSE,TRUE)
cairo_pdf('2d_ode_validation/figures/fig17_psi1_error_isolation.pdf',12,6.2,onefile=FALSE)
draw_all();dev.off()
png('2d_ode_validation/figures/fig17_psi1_error_isolation.png',width=1800,height=930,res=150)
draw_all();dev.off()
