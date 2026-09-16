library(data.table);library(ggplot2)
sim<-fread('2d_ode_validation/data/psi1_largeK_tau_results.csv')
curves<-fread('2d_ode_validation/data/psi1_largeK_2d_ode_curves.csv')[!is.na(P_kendall)]
dir.create('2d_ode_validation/figures',FALSE,TRUE)
theme_2d<-theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),
  strip.background=element_rect(fill='grey88'),plot.title=element_text(size=15),
  plot.subtitle=element_text(size=10),plot.caption=element_text(size=8,hjust=0),legend.position='bottom')
method_name<-'Numerical finite-horizon Kendall'
k_lab<-function(x)paste0('K = ',format(as.numeric(x),scientific=TRUE))

draw_probability<-function(kind,rho_value,theta_value){
  if(kind=='conditional'){
    measure<-'P_kendall';ttl<-'Conditional persistence probability'
    yl<-'Persistence probability, conditional on establishment'
    s<-copy(sim[rho==rho_value&theta==theta_value]);s[,`:=`(estimate=P_conditional,low=cond_low,high=cond_high)]
  }else{
    measure<-'P_kendall_unconditional';ttl<-'Unconditional persistence probability'
    yl<-'Unconditional persistence probability';s<-copy(sim[rho==rho_value&theta==theta_value])
    s[,`:=`(estimate=P_unconditional,low=uncond_low,high=uncond_high)]
  }
  a<-copy(curves[rho==rho_value&theta==theta_value])
  a[,`:=`(value=get(measure),method=method_name)]
  stub<-CJ(K=sort(unique(s$K)),method=method_name);stub[,`:=`(R0=NA_real_,value=NA_real_)]
  p<-ggplot()+geom_line(data=stub,aes(R0-1,value,colour=method),linewidth=.9,na.rm=TRUE)+
    geom_line(data=a,aes(R0-1,value,colour=method),linewidth=.9,na.rm=TRUE)+
    geom_errorbar(data=s,aes(R0-1,ymin=low,ymax=high),width=0,colour='grey25',linewidth=.42)+
    geom_point(data=s,aes(R0-1,estimate),colour='#2C7FB8',size=2)+
    facet_wrap(~K,ncol=2,labeller=labeller(K=k_lab))+
    scale_x_log10(breaks=c(.05,.1,.2,.5,1,2,5))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),expand=expansion(mult=c(.01,.03)))+
    scale_colour_manual(values=setNames('#009E73',method_name),limits=method_name,drop=FALSE)+
    labs(title=ttl,subtitle=sprintf('Full 2D ODE; psi = 1, rho = %.2f, theta = %g; I0 = 1',rho_value,theta_value),
      x=expression(R[0]-1~'(log scale)'),y=yl,colour=NULL,
      caption='Blue points and grey bars: new adaptive tau-leaping estimates (epsilon = 0.01) and 95% Wilson intervals.\nGreen: numerical finite-horizon Kendall. Curves stop where the finite-horizon interval does not exist.')+theme_2d
  if(kind=='unconditional'){ref<-unique(a[,.(K,R0)]);ref[,value:=1-1/R0]
    p<-p+geom_line(data=ref,aes(R0-1,value,group=K),colour='grey60',linetype='dashed',linewidth=.55)}
  p
}
write_pages<-function(kind,file){cairo_pdf(file,10.5,8.4,onefile=TRUE)
  for(rho in c(.01,.02,.05,.10))for(theta in c(0,.5,1))print(draw_probability(kind,rho,theta));dev.off()}
write_pages('unconditional','2d_ode_validation/figures/fig12_psi1_largeK_2d_ode_unconditional.pdf')
write_pages('conditional','2d_ode_validation/figures/fig13_psi1_largeK_2d_ode_conditional.pdf')
cat('psi=1 large-K Figures 12 and 13 complete\n')
