source('validation/R/theory.R')
stopifnot(requireNamespace('data.table',quietly=TRUE),
          requireNamespace('ggplot2',quietly=TRUE))
library(data.table)
library(ggplot2)

sim <- fread('validation/data/stochastic_scan_results.csv')
grid <- unique(sim[,.(rho,theta,K,R0,point_id,P_conditional,cond_low,cond_high)])

# Cache D_theta by (R0,theta); no stochastic calculation is performed here.
pars <- unique(grid[,.(R0,theta)])
pars[,`:=`(D=vapply(seq_len(.N),function(i) D_regularized(R0[i],theta[i]),0.),
           C=vapply(seq_len(.N),function(i) C_explicit(R0[i],theta[i]),0.),
           c_L=vapply(seq_len(.N),function(i) laplace_correction(R0[i],theta[i]),0.))]
grid <- merge(grid,pars,by=c('R0','theta'))

pred <- t(vapply(seq_len(nrow(grid)),function(i) {
  z <- bi_next_quantities(grid$R0[i],grid$rho[i],grid$theta[i],grid$K[i],D=grid$D[i])
  c(z$P_conditional,z$log_B)
},numeric(8)))
pred <- as.data.table(pred)
setnames(pred,c('P_leading','P_D_only','P_Laplace_only','P_combined',
                               'logB_leading','logB_D_only','logB_Laplace_only','logB_combined'))
grid <- cbind(grid,pred)
for(nm in c('leading','D_only','Laplace_only','combined'))
  grid[,(paste0('ae_',nm)):=abs(P_conditional-get(paste0('P_',nm)))]
grid[,`:=`(improvement=ae_leading-ae_combined,
           near_threshold=R0<=1.2)]

# Directly test the saddle expansion against the exact, boundary-independent
# scaled Kendall integral with lower endpoint x_f.
kg <- unique(grid[,.(rho,theta,R0,c_L)])
kg[,`:=`(logJ_exact=NA_real_,logJ_leading=NA_real_,logJ_next=NA_real_)]
for(i in seq_len(nrow(kg))) {
  R0 <- kg$R0[i]; rho <- kg$rho[i]; theta <- kg$theta[i]
  xf <- x_final(R0); xs <- 1/R0; hs <- h_theta(xs,theta)
  f <- function(X) vapply(X,function(xx) if(xx>=1) 0 else
    exp(action_diff(xs,xx,R0,theta)/rho)/(rho*h_theta(xx,theta)),0.)
  J <- integrate(f,xf,1,rel.tol=2e-9,abs.tol=0,subdivisions=1200L,
                 stop.on.error=TRUE)$value
  JL <- sqrt(2*pi/(R0*rho*hs))
  kg[i,`:=`(logJ_exact=log(J),logJ_leading=log(JL),
            logJ_next=log(JL)+log1p(rho*c_L))]
}
kg[,`:=`(relerr_J_leading=expm1(logJ_leading-logJ_exact),
         relerr_J_next=expm1(logJ_next-logJ_exact),
         abslog_J_improvement=abs(logJ_leading-logJ_exact)-abs(logJ_next-logJ_exact),
         near_threshold=R0<=1.2)]
fwrite(grid,'validation/data/BI_next_validation.csv')
fwrite(kg,'validation/data/BI_next_Kendall_validation.csv')

theme_bi <- theme_bw(base_size=9)+theme(panel.grid.minor=element_blank(),
  strip.background=element_rect(fill='grey90'),legend.position='bottom',
  plot.title=element_text(size=11,face='bold'),plot.subtitle=element_text(size=8.5))
cols <- c(Leading='#009E73',Combined='#D55E00')

show <- melt(grid,measure.vars=c('P_leading','P_combined'),
  variable.name='method',value.name='prediction')
show[,method:=factor(method,c('P_leading','P_combined'),c('Leading','Combined'))]
p1 <- ggplot(show,aes(P_conditional,prediction,colour=method))+
  geom_abline(slope=1,intercept=0,colour='grey55',linewidth=.45)+
  geom_segment(aes(x=cond_low,xend=cond_high,y=prediction,yend=prediction),
               colour='grey75',linewidth=.25)+
  geom_point(alpha=.67,size=1.15)+facet_grid(theta~rho,labeller=label_both)+
  scale_colour_manual(values=cols)+coord_equal(xlim=c(0,1),ylim=c(0,1))+
  labs(title='A  Conditional persistence: cached stochastic estimate vs BI',
       x='Stochastic conditional persistence',y='Analytical prediction',colour=NULL)+theme_bi

p2 <- ggplot(grid,aes(R0-1,improvement,colour=factor(rho)))+
  geom_hline(yintercept=0,colour='grey45',linewidth=.4)+geom_point(alpha=.7,size=1.25)+
  facet_wrap(~theta,nrow=1,labeller=label_both)+scale_x_log10()+
  scale_colour_brewer(palette='Dark2')+
  labs(title='B  Absolute-error improvement (positive favors combined)',
       x=expression(R[0]-1~'(log scale)'),y='|error leading| - |error combined|',colour=expression(rho))+theme_bi

effects <- melt(grid,measure.vars=c('ae_D_only','ae_Laplace_only','ae_combined'),
  variable.name='method',value.name='ae')
effects[,method:=factor(method,c('ae_D_only','ae_Laplace_only','ae_combined'),
  c('D_theta only','c_L only','Combined'))]
effects[,delta:=ae-ae_leading]
p3 <- ggplot(effects,aes(factor(rho),delta,fill=method))+
  geom_hline(yintercept=0,colour='grey45',linewidth=.4)+
  geom_boxplot(outlier.size=.35,linewidth=.35,position=position_dodge(width=.78))+
  facet_wrap(~theta,nrow=1,labeller=label_both)+
  scale_fill_manual(values=c('#0072B2','#CC79A7','#D55E00'))+
  labs(title='C  Individual correction effects and cancellation',
       subtitle='Negative change in absolute error is beneficial',x=expression(rho),
       y='Corrected absolute error - leading absolute error',fill=NULL)+theme_bi

kj <- melt(kg,measure.vars=c('relerr_J_leading','relerr_J_next'),
  variable.name='method',value.name='relative_error')
kj[,method:=factor(method,c('relerr_J_leading','relerr_J_next'),c('Leading','Next-order'))]
p4 <- ggplot(kj,aes(R0-1,relative_error,colour=method,linetype=factor(rho)))+
  geom_hline(yintercept=0,colour='grey55',linewidth=.35)+geom_line(linewidth=.55)+
  facet_wrap(~theta,nrow=1,labeller=label_both,scales='free_y')+scale_x_log10()+
  scale_colour_manual(values=c(Leading='#009E73','Next-order'='#D55E00'))+
  labs(title='D  Laplace approximation vs exact Kendall integral',
       subtitle='Relative error in the scaled integral; line type identifies rho',
       x=expression(R[0]-1~'(log scale)'),y='(approximation / exact) - 1',colour=NULL,linetype=expression(rho))+theme_bi

dir.create('validation/figures',FALSE,TRUE)
pdf_file <- 'validation/figures/fig14_BI_next_order_validation.pdf'
cairo_pdf(pdf_file,width=12,height=10.5,onefile=TRUE)
grid::grid.newpage()
lay <- grid::grid.layout(2,2,widths=grid::unit(c(1,1),'null'),heights=grid::unit(c(1,1),'null'))
grid::pushViewport(grid::viewport(layout=lay))
print(p1,vp=grid::viewport(layout.pos.row=1,layout.pos.col=1))
print(p2,vp=grid::viewport(layout.pos.row=1,layout.pos.col=2))
print(p3,vp=grid::viewport(layout.pos.row=2,layout.pos.col=1))
print(p4,vp=grid::viewport(layout.pos.row=2,layout.pos.col=2))
dev.off()

summ <- grid[,.(n=.N,MAE_leading=mean(ae_leading),MAE_D=mean(ae_D_only),
  MAE_Laplace=mean(ae_Laplace_only),MAE_combined=mean(ae_combined),
  median_improvement=median(improvement),fraction_improved=mean(improvement>0)),
  by=.(scope=ifelse(near_threshold,'R0<=1.2','R0>1.2'),rho)]
print(summ)
print(grid[,.(n=.N,MAE_leading=mean(ae_leading),MAE_combined=mean(ae_combined),
  fraction_improved=mean(improvement>0)),by=.(scope=ifelse(near_threshold,'R0<=1.2','R0>1.2'))])
print(kg[,.(mean_abs_rel_leading=mean(abs(relerr_J_leading)),
             mean_abs_rel_next=mean(abs(relerr_J_next)),
             fraction_next_better=mean(abslog_J_improvement>0)),by=rho])
cat('next-order BI validation complete:',pdf_file,'\n')
