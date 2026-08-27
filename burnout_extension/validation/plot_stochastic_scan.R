source('validation/R/scan.R')
stopifnot(requireNamespace('data.table',quietly=TRUE),
          requireNamespace('ggplot2',quietly=TRUE))
library(data.table)
library(ggplot2)

simulation_file <- 'validation/data/stochastic_scan_results.csv'
analytic_file <- 'validation/data/stochastic_scan_analytic.csv'
figure_dir <- 'validation/figures'
dir.create(figure_dir,FALSE,TRUE)

rho_values <- c(.01,.02,.05,.10)
theta_values <- c(0,.5,1)
K_values <- c(1000,3000,10000,30000)
R0_values <- 1+exp(seq(log(.05),log(5),length.out=101))
analytic_grid <- expand.grid(rho=rho_values,theta=theta_values,K=K_values,
  R0=R0_values,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
analytic_grid <- analytic_grid[order(analytic_grid$rho,analytic_grid$theta,
                                     analytic_grid$K,analytic_grid$R0),]
analytic_grid$row_id <- seq_len(nrow(analytic_grid))

required_columns <- c('K_first','K_second','L_first','L_second',
  'K_first_compromise','K_second_compromise','L_first_compromise','L_second_compromise',
  'K_first_compromise_3_4','K_second_compromise_3_4',
  'L_first_compromise_3_4','L_second_compromise_3_4',
  'K_first_star','K_second_star','L_first_star','L_second_star')
if(file.exists(analytic_file)) {
  analytic <- fread(analytic_file)
  if(!all(required_columns %in% names(analytic))) analytic <- data.table()
  else {
    if('row_id.1' %in% names(analytic)) analytic[,row_id.1:=NULL]
    key <- function(z) do.call(paste,c(z[,c('rho','theta','K','R0')],sep='|'))
    analytic[,row_id:=match(key(analytic),key(analytic_grid))]
    analytic <- analytic[is.finite(row_id)]
    setorder(analytic,row_id)
    analytic <- unique(analytic,by='row_id')
  }
} else {
  analytic <- data.table()
}
done <- if(nrow(analytic)) analytic$row_id else integer()
missing <- setdiff(analytic_grid$row_id,done)
if(length(missing)) {
  cl <- parallel::makeCluster(8L)
  parallel::clusterEvalQ(cl,source('validation/R/scan.R'))
  blocks <- split(missing,ceiling(seq_along(missing)/80L))
  for(ids in blocks) {
    got <- parallel::parLapply(cl,ids,function(i,analytic_grid) {
      g <- analytic_grid[i,c('rho','theta','K','R0')]
      p <- try({
        R0<-g$R0;rho<-g$rho;theta<-g$theta;K<-g$K
        ode<-ode_reference(R0,rho,theta,K)
        D<-D_regularized(R0,theta)
        m<-matching(R0,rho,theta,K,D=D)
        ys<-rho*h_theta(1/R0,theta)
        yc<-boundary_layer_height(ys,K,'compromise')
        mc<-matching(R0,rho,theta,K,D=D,yline=yc)
        yc34<-boundary_layer_height(ys,K,'compromise_3_4')
        mc34<-matching(R0,rho,theta,K,D=D,yline=yc34)
        ms<-matching(R0,rho,theta,K,D=D,yline=ys)
        probs<-function(mm,yline=NULL) {
          q1<-if(is.finite(mm$x_first)) kendall_quantities(mm$x_first,R0,rho,theta,K,yline) else NULL
          q2<-if(is.finite(mm$x_second)) kendall_quantities(mm$x_second,R0,rho,theta,K,yline) else NULL
          c(K_first=q1$P1 %||% NA_real_,K_second=q2$P1 %||% NA_real_,
            L_first=q1$P1_L %||% NA_real_,L_second=q2$P1_L %||% NA_real_)
        }
        keys<-c('K_first','K_second','L_first','L_second')
        c(probs(m),setNames(probs(mc,yc),paste0(keys,'_compromise')),
          setNames(probs(mc34,yc34),paste0(keys,'_compromise_3_4')),
          setNames(probs(ms,ys),paste0(keys,'_star')),
          status=ode$status)
      },silent=TRUE)
      if(inherits(p,'try-error')) data.frame(row_id=i,g,
        status=paste0('ERROR:',as.character(p)),K_first=NA_real_,
        K_second=NA_real_,L_first=NA_real_,L_second=NA_real_,
        K_first_compromise=NA_real_,K_second_compromise=NA_real_,
        L_first_compromise=NA_real_,L_second_compromise=NA_real_,
        K_first_compromise_3_4=NA_real_,K_second_compromise_3_4=NA_real_,
        L_first_compromise_3_4=NA_real_,L_second_compromise_3_4=NA_real_,
        K_first_star=NA_real_,K_second_star=NA_real_,
        L_first_star=NA_real_,L_second_star=NA_real_)
      else data.frame(row_id=i,g,status=unname(p['status']),
        as.list(as.numeric(p[c(keys,paste0(keys,'_compromise'),
          paste0(keys,'_compromise_3_4'),paste0(keys,'_star'))])) |>
          setNames(c(keys,paste0(keys,'_compromise'),
            paste0(keys,'_compromise_3_4'),paste0(keys,'_star'))))
    },analytic_grid=analytic_grid)
    analytic <- rbind(analytic,rbindlist(got,fill=TRUE),fill=TRUE)
    setorder(analytic,row_id);fwrite(analytic,analytic_file)
    cat('analytic curves',max(ids),'/',nrow(analytic_grid),'\n')
  }
  parallel::stopCluster(cl)
}

sim <- fread(simulation_file)
stopifnot(all(sim$attempts>=3000L),all(sim$unresolved==0L))

method_labels <- c(
  K_first='Exact Kendall + first-order entry',
  K_second='Exact Kendall + second-order entry',
  L_first='Laplace + first-order entry',
  L_second='Laplace + second-order entry')

boundary_math_label <- function(x) {
  labels <- c(
    'sqrt(y*/K)'='y[BL] == sqrt(y[star]/K)',
    '2/3 compromise'='y[BL] == K^{-1/3} * y[star]^{2/3}',
    '3/4 compromise'='y[BL] == K^{-1/4} * y[star]^{3/4}',
    'y*'='y[BL] == y[star]',
    'Early-establishment reference'="'Early-establishment reference'")
  parse(text=unname(labels[x]))
}

make_long <- function(probability=c('unconditional','conditional')) {
  probability <- match.arg(probability)
  a <- copy(analytic)
  if(probability=='unconditional') {
    factor_est <- 1-1/a$R0
    for(nm in required_columns) set(a,j=nm,value=a[[nm]]*factor_est)
  }
  standard <- names(method_labels)
  a <- melt(a,id.vars=c('row_id','rho','theta','K','R0','status'),
    measure.vars=c(standard,paste0(standard,'_compromise'),
      paste0(standard,'_compromise_3_4'),paste0(standard,'_star')),
    variable.name='variant',value.name='value')
  a[,boundary:=fcase(grepl('_compromise_3_4$',variant),'3/4 compromise',
    grepl('_compromise$',variant),'2/3 compromise',
    grepl('_star$',variant),'y*',default='sqrt(y*/K)')]
  a[,method_key:=sub('_(compromise_3_4|compromise|star)$','',variant)]
  a[,method:=factor(method_key,levels=names(method_labels),labels=method_labels)]
  a[,boundary:=factor(boundary,levels=c(
    'sqrt(y*/K)','2/3 compromise','3/4 compromise','y*'))]
  setorder(a,rho,theta,K,method,boundary,R0)
  a[,segment:=data.table::rleid(is.finite(value)),by=.(rho,theta,K,method,boundary)]
  a[is.finite(value)]
}

draw_page <- function(probability,rho_value,theta_value,boundary_value=NULL) {
  a <- make_long(probability)[rho==rho_value & theta==theta_value]
  if(!is.null(boundary_value)) a <- a[boundary==boundary_value]
  s <- sim[rho==rho_value & theta==theta_value]
  if(probability=='unconditional') {
    s[,`:=`(estimate=P_unconditional,low=uncond_low,high=uncond_high)]
    ylabel <- 'Unconditional persistence probability'
    title <- 'Full stochastic validation of unconditional persistence'
  } else {
    s[,`:=`(estimate=P_conditional,low=cond_low,high=cond_high)]
    ylabel <- 'Persistence probability, conditional on not fizzling'
    title <- 'Stochastic validation conditional on escaping early fizzle'
  }
  p <- ggplot()+
    geom_line(data=a,aes(R0-1,value,colour=method,linetype=boundary,
      group=interaction(method,boundary,segment)),linewidth=.75,na.rm=TRUE)+
    geom_errorbar(data=s,aes(R0-1,ymin=low,ymax=high),width=0,colour='grey25',
      linewidth=.42,na.rm=TRUE)+
    geom_point(data=s,aes(R0-1,estimate),colour='#2C7FB8',size=1.9,na.rm=TRUE)+
    facet_wrap(~K,ncol=2,labeller=labeller(K=function(x) paste0('K = ',x)))+
    scale_x_log10(breaks=c(.05,.1,.2,.5,1,2,5),
      labels=c('0.05','0.10','0.20','0.50','1','2','5'))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),expand=expansion(mult=c(.01,.03)))+
    scale_colour_manual(values=c(
      'Exact Kendall + first-order entry'='black',
      'Exact Kendall + second-order entry'='#D95F02',
      'Laplace + first-order entry'='grey35',
      'Laplace + second-order entry'='#7570B3',
      'Early-establishment reference'='grey60'))+
    scale_linetype_manual(values=c('sqrt(y*/K)'='solid',
      '2/3 compromise'='42','y*'='22',
      '3/4 compromise'='44',
      'Early-establishment reference'='solid'),labels=boundary_math_label)+
    labs(title=title,
      subtitle=sprintf('Exact Gillespie CTMC; rho = %.2f, theta = %g; I0 = 1%s',
        rho_value,theta_value,if(is.null(boundary_value)) '' else
          paste0('; matching height = ',boundary_value)),
      x=expression(R[0]-1~'(log scale)'),y=ylabel,
      colour='Analytical method',linetype='Matching height',
      caption=paste0('Points and 95% Wilson intervals are stochastic estimates. ',
        if(is.null(boundary_value))
          'Curves compare sqrt(y*/K), the 2/3 and 3/4 compromises, and y*.'
        else paste0('All analytical curves use matching height ',boundary_value,'.'),'\n',
        'Each curve stops where its matching equation has no admissible root.'))+
    guides(colour=guide_legend(nrow=2,byrow=TRUE,order=1),
      linetype=guide_legend(nrow=1,order=2))+
    theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),
      strip.background=element_rect(fill='grey88'),legend.position='bottom',
      legend.box='vertical',legend.box.just='left',
      plot.title=element_text(size=15),plot.subtitle=element_text(size=10),
      plot.caption=element_text(size=8,hjust=0))
  if(probability=='unconditional') {
    ref <- unique(a[,.(R0)])[order(R0)]
    ref[,value:=1-1/R0]
    p <- p+geom_line(data=ref,aes(R0-1,value,
      colour='Early-establishment reference',
      linetype='Early-establishment reference'),inherit.aes=FALSE,linewidth=.6)
  }
  p
}

write_scan_pdf <- function(probability,file,boundary_value=NULL) {
  grDevices::cairo_pdf(file,width=10.5,height=8.4,onefile=TRUE)
  for(rho_value in rho_values) for(theta_value in theta_values)
    print(draw_page(probability,rho_value,theta_value,boundary_value))
  dev.off()
}

fig09_dir <- file.path(figure_dir,'fig09')
fig10_dir <- file.path(figure_dir,'fig10')
dir.create(fig09_dir,FALSE,TRUE)
dir.create(fig10_dir,FALSE,TRUE)
boundary_files <- c('sqrt(y*/K)'='sqrt','2/3 compromise'='2_3_compromise',
  '3/4 compromise'='3_4_compromise','y*'='ystar')
for(boundary_value in names(boundary_files)) {
  suffix <- unname(boundary_files[boundary_value])
  write_scan_pdf('unconditional',file.path(fig09_dir,
    paste0('P1_unconditional_',suffix,'.pdf')),boundary_value)
  write_scan_pdf('conditional',file.path(fig10_dir,
    paste0('P1_conditional_',suffix,'.pdf')),boundary_value)
}
cat('stochastic scan PDFs complete\n')
