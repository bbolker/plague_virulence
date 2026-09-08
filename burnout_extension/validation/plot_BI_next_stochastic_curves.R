source('validation/R/theory.R')
stopifnot(requireNamespace('data.table',quietly=TRUE),
          requireNamespace('ggplot2',quietly=TRUE))
library(data.table)
library(ggplot2)

sim <- fread('validation/data/stochastic_scan_results.csv')
analytic <- unique(fread('validation/data/stochastic_scan_analytic.csv')[,
  .(rho,theta,K,R0)])

add_bi_curves <- function(z) {
  keys <- unique(z[,.(R0,theta)])
  keys[,D:=vapply(seq_len(.N),function(i) D_regularized(R0[i],theta[i]),0.)]
  z <- merge(z,keys,by=c('R0','theta'))
  vals <- t(vapply(seq_len(nrow(z)),function(i) {
    q <- bi_next_quantities(z$R0[i],z$rho[i],z$theta[i],z$K[i],D=z$D[i])
    q$P_conditional[c('leading','combined')]
  },numeric(2)))
  z[,`:=`(leading=vals[,1],combined=vals[,2])]
  z
}
analytic <- add_bi_curves(analytic)

draw_page <- function(rho_value,theta_value,curves=analytic,points=sim,
                      large_K=FALSE) {
  a <- melt(curves[rho==rho_value & theta==theta_value],
    id.vars=c('rho','theta','K','R0'),measure.vars=c('leading','combined'),
    variable.name='method',value.name='value')
  a[,method:=factor(method,c('leading','combined'),
    c('Leading BI','Combined next-order BI'))]
  s <- copy(points[rho==rho_value & theta==theta_value])
  ggplot()+
    geom_line(data=a,aes(R0-1,value,colour=method),linewidth=.9)+
    geom_errorbar(data=s,aes(R0-1,ymin=cond_low,ymax=cond_high),width=0,
      colour='grey25',linewidth=.42)+
    geom_point(data=s,aes(R0-1,P_conditional),colour='#2C7FB8',size=2)+
    facet_wrap(~K,ncol=2,labeller=labeller(K=function(x) paste0('K = ',
      if(large_K) format(as.numeric(x),scientific=TRUE) else x)))+
    scale_x_log10(breaks=c(.03,.05,.1,.2,.5,1,2,5),
      labels=c('0.03','0.05','0.10','0.20','0.50','1','2','5'))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),
      expand=expansion(mult=c(.01,.03)))+
    scale_colour_manual(values=c('Leading BI'='#009E73',
                                  'Combined next-order BI'='#D55E00'))+
    labs(title='Conditional persistence probability',
      subtitle=sprintf('Boundary-independent comparison; rho = %.2f, theta = %g; I0 = 1',
                       rho_value,theta_value),
      x=expression(R[0]-1~'(log scale)'),
      y='Persistence probability, conditional on not fizzling',colour=NULL,
      caption=paste0('Blue points and grey bars are the existing stochastic estimates and 95% Wilson intervals. ',
        'Green: leading BI; orange: combined next-order BI. No simulations were rerun.'))+
    theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),
      strip.background=element_rect(fill='grey88'),plot.title=element_text(size=15),
      plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0),
      legend.position='bottom')
}

dir.create('validation/figures',FALSE,TRUE)
outfile <- 'validation/figures/fig15_BI_next_order_stochastic_curves.pdf'
cairo_pdf(outfile,width=10.5,height=8.4,onefile=TRUE)
first_page <- TRUE
for(rho_value in c(.01,.02,.05,.10)) for(theta_value in c(0,.5,1)) {
  if(!first_page) grid::grid.newpage()
  print(draw_page(rho_value,theta_value),newpage=FALSE)
  first_page <- FALSE
}
dev.off()
cat('Fig. 12-style next-order stochastic comparison complete:',outfile,'\n')
