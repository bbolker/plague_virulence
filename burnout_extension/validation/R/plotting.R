library(ggplot2)
library(data.table)
theme_set(theme_bw(base_size=9)+theme(panel.grid.minor=element_blank(),legend.position='bottom'))

metric_specs <- list(
  xin_first=list(label='First-order entry x_in',estimate='x_in_first',truth='x_in_ODE',
    formula="E[abs] == abs(x['in']^{(1)}-x['in']^{plain(ODE)}) ~~ ',' ~~ E[rel] == frac(abs(x['in']^{(1)}-x['in']^{plain(ODE)}),abs(x['in']^{plain(ODE)}))"),
  xin_second=list(label='Second-order entry x_in',estimate='x_in_second',truth='x_in_ODE',
    formula="E[abs] == abs(x['in']^{(2)}-x['in']^{plain(ODE)}) ~~ ',' ~~ E[rel] == frac(abs(x['in']^{(2)}-x['in']^{plain(ODE)}),abs(x['in']^{plain(ODE)}))"),
  Lambda_first=list(label='First-order action Lambda',estimate='Lambda_first',truth='Lambda_ODE',
    formula="E[abs] == abs(Lambda^{(1)}-Lambda^{plain(ODE)}) ~~ ',' ~~ E[rel] == frac(abs(Lambda^{(1)}-Lambda^{plain(ODE)}),abs(Lambda^{plain(ODE)}))"),
  Lambda_second=list(label='Second-order action Lambda',estimate='Lambda_second',truth='Lambda_ODE',
    formula="E[abs] == abs(Lambda^{(2)}-Lambda^{plain(ODE)}) ~~ ',' ~~ E[rel] == frac(abs(Lambda^{(2)}-Lambda^{plain(ODE)}),abs(Lambda^{plain(ODE)}))"),
  B_K_first=list(label='Exact-Kendall B with first-order entry',estimate='B_K_first',truth='B_ref',
    definition="B == -log(P[burn]) ~~ ',' ~~ B == K*y[BL]*log(1+1/cal(I))",
    formula="E[abs] == abs(B[K]^{(1)}-B[ref]) ~~ ',' ~~ E[rel] == frac(abs(B[K]^{(1)}-B[ref]),abs(B[ref]))"),
  B_K_second=list(label='Exact-Kendall B with second-order entry',estimate='B_K_second',truth='B_ref',
    definition="B == -log(P[burn]) ~~ ',' ~~ B == K*y[BL]*log(1+1/cal(I))",
    formula="E[abs] == abs(B[K]^{(2)}-B[ref]) ~~ ',' ~~ E[rel] == frac(abs(B[K]^{(2)}-B[ref]),abs(B[ref]))"),
  B_L_ODE=list(label='Laplace-only B (ODE entry)',estimate='B_L_ODE',truth='B_ref',
    definition="B == -log(P[burn]) ~~ ',' ~~ B == K*y[BL]*log(1+1/cal(I))",
    formula="E[abs] == abs(B['L,ODE']-B[ref]) ~~ ',' ~~ E[rel] == frac(abs(B['L,ODE']-B[ref]),abs(B[ref]))")
)

paper_names <- c(
  xin_first='fig01_xin_first_error_maps.pdf',
  xin_second='fig02_xin_second_error_maps.pdf',
  Lambda_first='fig03_Lambda_first_error_maps.pdf',
  Lambda_second='fig04_Lambda_second_error_maps.pdf',
  B_K_first='fig05_B_exactKendall_first_error_maps.pdf',
  B_K_second='fig06_B_exactKendall_second_error_maps.pdf',
  B_L_ODE='fig07_B_Laplace_only_error_maps.pdf'
)

prepare_metric <- function(d,spec) {
  z<-copy(d)
  z[,boundary:=factor(boundary,levels=c(
    'sqrt(y*/K)','2/3 compromise','3/4 compromise','y*'))]
  z[,estimate:=get(spec$estimate)];z[,truth:=get(spec$truth)]
  z[,absolute_error:=abs(estimate-truth)]
  z[,relative_error:=absolute_error/abs(truth)]
  z[,reason:=fcase(status=='NO_ENTRY','NO_ENTRY',status!='OK',status,
                   is.na(estimate),'ANALYTIC_NO_ROOT',is.na(truth),'NUMERICAL_MISSING',default='VALID')]
  z[,u:=log10(R0-1)];z[,v:=log10(K)]
  z
}

draw_metric_page <- function(z,spec,rho_value) {
  zz<-z[abs(rho-rho_value)<1e-12]
  tick_values<-c(.05,.1,.2,.5,1,2,5,9)
  xb<-log10(tick_values)
  xl<-format(tick_values,trim=TRUE,scientific=FALSE)
  one<-function(column,label,show_x) {
    q<-zz[,.(boundary,theta,u,K_f=factor(K,levels=sort(unique(zz$K))),reason,value=get(column))]
    finite_values<-q[is.finite(value)&value>0,value]
    lower_power<-floor(log10(min(finite_values)))
    upper_power<-ceiling(log10(max(finite_values)))
    colour_limits<-10^c(lower_power,upper_power)
    exponent_span<-upper_power-lower_power
    if(exponent_span<=10) {
      exponent_step<-1
    } else {
      divisors<-which(exponent_span%%seq_len(exponent_span)==0)
      exponent_step<-divisors[which(exponent_span/divisors<=10)[1]]
    }
    exponent_breaks<-seq(lower_power,upper_power,by=exponent_step)
    colour_breaks<-10^exponent_breaks
    ggplot(q,aes(u,K_f,fill=value))+geom_tile()+
      facet_grid(boundary~theta,labeller=labeller(
        theta=function(x) paste0('theta = ',x),boundary=label_value))+
      scale_fill_viridis_c(trans='log10',option='viridis',na.value='grey82',
        limits=colour_limits,breaks=colour_breaks,
        labels=scales::label_scientific(digits=2),oob=scales::squish)+
      scale_x_continuous(breaks=xb,labels=if(show_x) xl else NULL,expand=c(0,0))+
      scale_y_discrete(labels=c('1e3','3e3','1e4','3e4','1e5'),expand=c(0,0))+
      labs(x=if(show_x) expression(R[0]-1) else NULL,
           y='K',fill=label)+
      theme_bw(base_size=9)+theme(panel.grid=element_blank(),panel.spacing=unit(.38,'cm'),
            panel.border=element_rect(colour='grey55',linewidth=.45),
            strip.background=element_blank(),strip.text=element_text(size=9),
            legend.position='right',legend.key.height=unit(1.05,'cm'),
            axis.text.x=element_text(size=7),aspect.ratio=.54,
            plot.margin=margin(2,3,2,3))
  }
  title_panel<-cowplot::ggdraw()+cowplot::draw_label(
    sprintf('%s    rho = %.2f',spec$label,rho_value),fontface='bold',size=14,x=.5,hjust=.5)
  formula_panel<-ggplot()+theme_void()+
    annotate('text',x=.5,y=if(is.null(spec$definition)) .5 else .05,
             label=spec$formula,parse=TRUE,size=if(is.null(spec$definition)) 4.6 else 4.3,colour='black')+
    coord_cartesian(xlim=c(0,1),ylim=c(-.5,1.5),clip='off')+
    theme(plot.margin=margin(4,4,4,4))
  if(!is.null(spec$definition)) formula_panel<-formula_panel+
    annotate('text',x=.5,y=.95,label=spec$definition,parse=TRUE,size=4.3,colour='black')
  top<-cowplot::plot_grid(title_panel,formula_panel,ncol=1,
                         rel_heights=if(is.null(spec$definition)) c(.42,.58) else c(.3,.7))
  header_height<-if(is.null(spec$definition)) .20 else .29
  cowplot::plot_grid(top,one('absolute_error','Absolute error',FALSE),
                     one('relative_error','Relative error',TRUE),ncol=1,
                     rel_heights=c(header_height,1,1.08))
}

write_metric_pdf <- function(d,name,spec,outdir='validation/figures') {
  z<-prepare_metric(d,spec);file<-file.path(outdir,unname(paper_names[name]))
  grDevices::cairo_pdf(file,width=11.6,height=21,onefile=TRUE)
  for(rv in c(.02,.04,.06,.08,.10)) print(draw_metric_page(z,spec,rv))
  dev.off();file
}

write_stochastic_pdf <- function(outdir='validation/figures') {
  s<-fread('validation/data/stochastic_results.csv')
  s[,panel:=sprintf('R0=%g, rho=%g, theta=%g, K=%g',R0,rho,theta,K)]
  a<-melt(s,id.vars=c('panel','P1_hat','CI_low','CI_high'),
          measure.vars=c('P1_uncond_K_ODE','P1_uncond_K_first','P1_uncond_K_second','P1_uncond_L_second'),
          variable.name='method',value.name='value')
  a[,method:=factor(method,
    levels=c('P1_uncond_K_ODE','P1_uncond_K_first','P1_uncond_K_second','P1_uncond_L_second'),
    labels=c('Kendall + ODE entry','Kendall + first-order entry',
             'Kendall + second-order entry','Laplace + second-order entry'))]
  sim<-unique(s[,.(panel,value=P1_hat,CI_low,CI_high)])
  p<-ggplot(a,aes(method,value,colour=method))+geom_point(size=2.2)+
    geom_point(data=sim,aes(x='Gillespie',y=value),inherit.aes=FALSE,size=2.4)+
    geom_errorbar(data=sim,aes(x='Gillespie',ymin=CI_low,ymax=CI_high),inherit.aes=FALSE,width=.15)+
    facet_wrap(~panel,ncol=3,scales='free_y')+labs(
      title=expression(paste('Unconditional first-trough survival: Gillespie versus analytical ',P[1])),
      subtitle=expression(paste(P[1]==(1-1/R[0])*(1-exp(-B)),
        '; 10,000 trajectories per point; persistence requires formation of the second infectious peak')),
      x=NULL,y=expression(P[1]==Pr('persist beyond the first trough')))+
    theme(axis.text.x=element_text(angle=35,hjust=1),legend.position='none')
  ggsave(file.path(outdir,'fig08_P1_unconditional_stochastic_vs_analytic.pdf'),p,width=10,height=4.2,device=cairo_pdf)
}

make_current_plots <- function() {
  dir.create('validation/figures',FALSE,TRUE)
  d<-fread('validation/data/K_R0_scan.csv')
  files<-mapply(function(n,s) write_metric_pdf(d,n,s),names(metric_specs),metric_specs,SIMPLIFY=FALSE)
  write_stochastic_pdf();invisible(files)
}
