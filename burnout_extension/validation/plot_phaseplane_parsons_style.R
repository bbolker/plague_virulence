# Same graphical elements as Parsons et al. (2024), Fig. 1D,E.
# Generalized model, own illustrative parameters and own boundary-layer threshold.
source('validation/R/theory.R');source('validation/R/native_simulation.R');library(deSolve)
R0<-1.5;rho<-.01;theta<-.5;K<-1e6;I0<-3;tf<-240
xs<-1/R0;ystar<-rho*h_theta(xs,theta);ybl<-sqrt(ystar/K)
rhs<-function(t,z,p)list(c(rho*h_theta(z[1],theta)-R0*z[1]*exp(z[2]),R0*z[1]-1))
z<-as.data.frame(ode(c(x=1-I0/K,logy=log(I0/K)),seq(0,tf,by=.01),rhs,NULL,
 method='lsoda',rtol=2e-11,atol=c(1e-13,1e-12),maxsteps=1e6));z$y<-exp(z$logy)
n<-nrow(z);peak<-which(z$x[-n]>=xs & z$x[-1]<xs)[1]
entry<-which(seq_len(n-1)>peak & z$y[-n]>=ybl & z$y[-1]<ybl)[1]
stopifnot(is.finite(entry))
w<-(ybl-z$y[entry])/(z$y[entry+1]-z$y[entry]);xin<-z$x[entry]+w*(z$x[entry+1]-z$x[entry])
g<-rho*h_theta(z$x,theta)-R0*z$x*z$y
crit<-sort(unique(c(which(z$x[-n]>=xs & z$x[-1]<xs),which(z$x[-n]<xs & z$x[-1]>=xs),which(g[-n]*g[-1]<0))))
dll<-compile_native_simulation('validation/src/phaseplane_sample_paths.c');dyn.load(dll)
set.seed(20261107);np<-20L;dt<-.1;nt<-as.integer(tf/dt+1)
paths<-.C('phaseplane_paths',as.double(c(R0,rho,theta,K,I0)),np,nt,as.double(dt),sx=double(np*nt),iy=double(np*nt))
X<-matrix(paths$sx,nrow=nt);Y<-matrix(paths$iy,nrow=nt)
out<-'validation/figures/previews';dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(capabilities('aqua'))quartzFonts(STIX=quartzFont(c('STIXGeneral-Regular','STIXGeneral-Bold','STIXGeneral-Italic','STIXGeneral-BoldItalic')))
library(ggplot2)
library(data.table)
panel_levels<-c('(a) Linear scale','(b) Log scale')
transform_panel<-function(d){
 rbindlist(lapply(seq_along(panel_levels),function(k){
  q<-copy(as.data.table(d));q[,panel:=factor(panel_levels[k],levels=panel_levels)]
  q[,yp:=if(k==1L)y else log10(pmax(y,1e-8))];q
 }))
}
ode_plot<-transform_panel(data.frame(x=z$x,y=z$y,group='ODE',kind='ODE trajectory'))
stoch_plot<-rbindlist(lapply(seq_len(np),function(j){
 zero<-which(Y[,j]==0);ids<-seq_len(if(length(zero))zero[1]else nrow(Y))
 data.table(x=X[ids,j],y=Y[ids,j],group=paste0('path',j),kind='Stochastic trajectories')
}))
stoch_plot<-transform_panel(stoch_plot)
xx<-seq(.38,1,length.out=1000)
null_plot<-transform_panel(data.frame(x=xx,y=rho*h_theta(xx,theta)/(R0*xx),group='nullcline',kind='Susceptible nullcline'))
entry_plot<-transform_panel(data.frame(x=xin,y=ybl))
rect_plot<-data.table(panel=factor(panel_levels,levels=panel_levels),xmin=.38,xmax=1.01,
 ymin=c(0,-8),ymax=c(ybl,log10(ybl)),kind='Boundary layer')
limits<-data.table(panel=factor(rep(panel_levels,each=2),levels=panel_levels),x=rep(c(.38,1.01),2),yp=c(0,.085,-8,log10(.5)))
p<-ggplot()+
 geom_rect(data=rect_plot,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=kind),colour=NA)+
 geom_path(data=stoch_plot,aes(x,yp,group=group,colour=kind),linewidth=.3,alpha=.4)+
 geom_line(data=null_plot,aes(x,yp,colour=kind),linewidth=.45)+
 geom_path(data=ode_plot,aes(x,yp,colour=kind),linewidth=.9)+
 geom_point(data=entry_plot,aes(x,yp,shape='Entry'),fill='#FFE34D',colour='black',size=3,stroke=.45)+
 geom_blank(data=limits,aes(x,yp))+
 facet_wrap(~panel,nrow=1,scales='free_y')+
 scale_x_continuous(breaks=seq(.4,1,.1),limits=c(.38,1.01),expand=c(0,0))+
 scale_y_continuous(breaks=function(v)if(min(v)< -1) -8:-1 else seq(0,.08,.02),
  labels=function(v)if(any(v<0))parse(text=paste0('10^',v))else format(v,nsmall=2),expand=c(0,0))+
 scale_colour_manual(values=c('ODE trajectory'='#EE3838','Stochastic trajectories'='grey65','Susceptible nullcline'='#CE54C7'),
  breaks=c('ODE trajectory','Stochastic trajectories','Susceptible nullcline'),name=NULL)+
 scale_fill_manual(values=c('Boundary layer'='#FFF8D8'),labels=expression(paste('Boundary layer: ',y<y[BL])),name=NULL)+
 scale_shape_manual(values=c('Entry'=21),labels=expression(paste('Entry: (',x['in'],', ',y[BL],')')),name=NULL)+
 guides(colour=guide_legend(order=1,override.aes=list(alpha=1)),
        fill=guide_legend(order=2),shape=guide_legend(order=3,override.aes=list(fill='#FFE34D')))+
 labs(x='x',y='y',
  subtitle=bquote(.(intToUtf8(0x211B))[0]==.(R0)~','~rho==.(rho)~','~theta==.(theta)~','~K==10^6~','~I[0]==.(I0)),
  caption=NULL)+
 theme_bw(base_size=10.5)+
 theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(linewidth=.25,colour='grey88'),
  strip.background=element_rect(fill='grey93',colour='grey55'),strip.text=element_text(face='bold'),
  legend.position='bottom',legend.box='horizontal',legend.text=element_text(size=9),
  plot.margin=margin(5.5,7,5.5,5.5),plot.caption=element_text(hjust=0,size=8),
  plot.subtitle=element_text(family=if(capabilities('aqua'))'STIX'else 'serif',size=10),
  panel.spacing=grid::unit(1.1,'lines'))
if(capabilities('aqua')){
 quartz(file=file.path(out,'phaseplane_framework_preview.pdf'),type='pdf',width=9,height=4.9);print(p);dev.off()
 png(file.path(out,'phaseplane_framework_preview.png'),type='quartz',width=1800,height=980,res=200);print(p);dev.off()
}else{
 ggsave(file.path(out,'phaseplane_framework_preview.pdf'),p,width=9,height=4.9,device=cairo_pdf)
 ggsave(file.path(out,'phaseplane_framework_preview.png'),p,width=9,height=4.9,dpi=200)
}

write.csv(data.frame(R0,rho,theta,K,I0,tf,xs,ystar,ybl,xin,seed=20261107),
 'validation/data/phaseplane_parsons_style_parameters.csv',row.names=FALSE)
saveRDS(list(ODE=z,stochastic=list(x=X,y=Y,dt=dt),parameters=c(R0=R0,rho=rho,theta=theta,K=K,I0=I0)),
 'validation/data/phaseplane_parsons_style_paths.rds')
cat('threshold',ybl,'entry',xin,'20 fixed-seed exact sample paths\n')

if ('--paper-export' %in% commandArgs(trailingOnly=TRUE)) {
 stopifnot(file.copy(file.path(out,'phaseplane_framework_preview.pdf'),
  'validation/figures/paper/phaseplane_framework.pdf',overwrite=TRUE))
}
