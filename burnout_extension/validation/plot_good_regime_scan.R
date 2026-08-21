source('validation/R/scan.R')
stopifnot(requireNamespace('data.table',quietly=TRUE),
          requireNamespace('ggplot2',quietly=TRUE))
library(data.table)
library(ggplot2)

simulation_file <- 'validation/data/good_regime_scan_results.csv'
analytic_file <- 'validation/data/good_regime_scan_analytic.csv'
output_file <- 'validation/figures/fig11_P1_good_regime_stochastic_validation.pdf'
K_values <- c(1e6,1e7,1e8,1e9)
delta_values <- exp(seq(log(.03),log(5),length.out=181))
grid <- expand.grid(rho=.01,theta=c(0,.5,1),K=K_values,
  R0=1+delta_values,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
row.names(grid) <- NULL
grid <- grid[order(grid$theta,grid$K,grid$R0),]

analytic <- if(file.exists(analytic_file)) fread(analytic_file) else data.table()
key <- function(z) paste(z$rho,z$theta,z$K,format(z$R0,digits=17),sep='|')
if(nrow(analytic)) analytic <- analytic[key(analytic) %in% key(grid)]
done <- if(nrow(analytic)) unique(key(analytic)) else character()
missing <- which(!key(grid) %in% done)
if(length(missing)) {
  constant_grid <- unique(grid[missing,c('theta','R0')])
  constants <- mapply(function(R0,theta) D_regularized(R0,theta),
    constant_grid$R0,constant_grid$theta,SIMPLIFY=FALSE)
  names(constants) <- paste(constant_grid$theta,
    format(constant_grid$R0,digits=17),sep='|')
  one <- function(g) {
    R0<-g$R0;rho<-g$rho;theta<-g$theta;K<-g$K
    D <- constants[[paste(theta,format(R0,digits=17),sep='|')]]
    ystar <- rho*h_theta(1/R0,theta)
    calc <- function(yline,label) {
      m <- try(matching(R0,rho,theta,K,D=D,yline=yline),silent=TRUE)
      if(inherits(m,'try-error')) return(data.frame())
      q <- function(x) if(is.finite(x)) kendall_quantities(x,R0,rho,theta,K,yline) else NULL
      q1<-q(m$x_first);q2<-q(m$x_second)
      data.frame(boundary=label,
        method=c('Exact Kendall + first-order entry',
                 'Exact Kendall + second-order entry',
                 'Laplace + first-order entry',
                 'Laplace + second-order entry'),
        value=c(q1$P1 %||% NA_real_,q2$P1 %||% NA_real_,
                q1$P1_L %||% NA_real_,q2$P1_L %||% NA_real_))
    }
    rbind(transform(calc(NULL,'Standard: sqrt(y*/K)'),rho=rho,theta=theta,K=K,R0=R0),
          transform(calc(ystar,'Alternative: y*'),rho=rho,theta=theta,K=K,R0=R0))
  }
  blocks <- split(missing,ceiling(seq_along(missing)/40L))
  for(ids in blocks) {
    analytic <- rbind(analytic,rbindlist(lapply(ids,function(i) one(grid[i,]))),fill=TRUE)
    fwrite(analytic,analytic_file)
    cat('analytic',max(ids),'/',nrow(grid),'\n')
  }
}

sim <- fread(simulation_file)
stopifnot(nrow(sim)==168L,all(sim$attempts>=3000L),all(sim$unresolved==0L))
analytic[,method:=factor(method,levels=c(
  'Exact Kendall + first-order entry','Exact Kendall + second-order entry',
  'Laplace + first-order entry','Laplace + second-order entry'))]
analytic[,boundary:=factor(boundary,levels=c('Standard: sqrt(y*/K)','Alternative: y*'))]
setorder(analytic,theta,K,method,boundary,R0)
analytic[,segment:=rleid(is.finite(value)),by=.(theta,K,method,boundary)]

draw_page <- function(probability=c('unconditional','conditional'),theta_value) {
  probability <- match.arg(probability)
  a <- copy(analytic[theta==theta_value])
  s <- copy(sim[theta==theta_value])
  if(probability=='unconditional') {
    a[,value:=value*(1-1/R0)]
    s[,`:=`(estimate=P_unconditional,low=uncond_low,high=uncond_high)]
    title <- 'Unconditional persistence probability'
    ylabel <- 'Unconditional persistence probability'
  } else {
    s[,`:=`(estimate=P_conditional,low=cond_low,high=cond_high)]
    title <- 'Persistence probability conditional on escaping early fizzle'
    ylabel <- 'Persistence probability, conditional on not fizzling'
  }
  p <- ggplot()+
    geom_line(data=a[is.finite(value)],aes(R0-1,value,colour=method,
      linetype=boundary,group=interaction(method,boundary,segment)),linewidth=.78)+
    geom_errorbar(data=s,aes(R0-1,ymin=low,ymax=high),width=0,
      colour='grey25',linewidth=.42)+
    geom_point(data=s,aes(R0-1,estimate),colour='#2C7FB8',size=2)+
    facet_wrap(~K,ncol=2,labeller=labeller(K=function(x)
      paste0('K = ',format(as.numeric(x),scientific=TRUE))))+
    scale_x_log10(breaks=c(.03,.05,.1,.3,.5,1,2,5),
      labels=c('.03','.05','.10','.30','.50','1','2','5'))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),
      expand=expansion(mult=c(.01,.03)))+
    scale_colour_manual(values=c('Exact Kendall + first-order entry'='black',
      'Exact Kendall + second-order entry'='#D95F02',
      'Laplace + first-order entry'='grey45',
      'Laplace + second-order entry'='#7570B3'))+
    scale_linetype_manual(values=c('Standard: sqrt(y*/K)'='solid',
                                   'Alternative: y*'='22'))+
    labs(title=title,subtitle=sprintf('rho = 0.01, theta = %g, I0 = 1',theta_value),
      x=expression(R[0]-1~'(log scale)'),y=ylabel,
      colour='Analytical method',linetype='Matching height',
      caption=paste0('Points and 95% Wilson intervals: adaptive tau-leaping (epsilon = 0.01); 3,000 attempts per point, expanded to 10,000 when needed.\n',
        'Curves use the analytical D_theta expression. Missing curve segments indicate no admissible matching root.'))+
    guides(colour=guide_legend(nrow=2,byrow=TRUE,order=1),
      linetype=guide_legend(nrow=1,order=2))+
    theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),
      strip.background=element_rect(fill='grey88'),legend.position='bottom',
      legend.box='vertical',legend.box.just='left',plot.title=element_text(size=15),
      plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0))
  if(probability=='unconditional') {
    ref <- unique(a[,.(K,R0)]);ref[,value:=1-1/R0]
    p <- p+geom_line(data=ref,aes(R0-1,value,group=K),colour='grey65',linewidth=.55)
  }
  p
}

dir.create(dirname(output_file),FALSE,TRUE)
cairo_pdf(output_file,width=10.5,height=8.4,onefile=TRUE)
for(theta_value in c(0,.5,1)) {
  print(draw_page('unconditional',theta_value))
  print(draw_page('conditional',theta_value))
}
dev.off()
cat('wrote',output_file,'\n')
