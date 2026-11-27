# Same graphical elements as Parsons et al. (2024), Fig. 1D,E.
# Generalized model, own illustrative parameters and own boundary-layer threshold.
source('validation/R/theory.R');library(deSolve)
R0<-1.5;rho<-.01;theta<-.5;K<-1e6;I0<-3;tf<-240
xs<-1/R0;ystar<-rho*h_theta(xs,theta);ybl<-sqrt(ystar/K)
rhs<-function(t,z,p)list(c(rho*h_theta(z[1],theta)-R0*z[1]*exp(z[2]),R0*z[1]-1))
z<-as.data.frame(ode(c(x=1-I0/K,logy=log(I0/K)),seq(0,tf,by=.01),rhs,NULL,
 method='lsoda',rtol=2e-11,atol=c(1e-13,1e-12),maxsteps=1e6));z$y<-exp(z$logy)
n<-nrow(z);peak<-which(z$x[-n]>=xs & z$x[-1]<xs)[1]
g<-rho*h_theta(z$x,theta)-R0*z$x*z$y
out<-'validation/figures/paper';dir.create(out,recursive=TRUE,showWarnings=FALSE)
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
xx<-seq(.38,1,length.out=1000)
null_plot<-transform_panel(data.frame(x=xx,y=rho*h_theta(xx,theta)/(R0*xx),group='nullcline',kind='Susceptible nullcline'))
rect_plot<-data.table(panel=factor(panel_levels,levels=panel_levels),xmin=.38,xmax=1.01,
 ymin=c(0,-8),ymax=c(ybl,log10(ybl)),kind='Boundary layer')
limits<-data.table(panel=factor(rep(panel_levels,each=2),levels=panel_levels),x=rep(c(.38,1.01),2),yp=c(0,.085,-8,log10(.5)))
# Mark the first susceptible minimum on the logarithmic panel only.
turn<-which(seq_len(n-1)>peak & g[-n]<0 & g[-1]>=0)[1]
phi<-seq(0,2*pi,length.out=200)
corner_plot<-data.table(panel=factor('(b) Log scale',levels=panel_levels),
 x=z$x[turn]+.035*cos(phi),yp=log10(z$y[turn])+.6*sin(phi),kind='Corner region')
p<-ggplot()+
 geom_rect(data=rect_plot,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=kind),colour=NA)+
 geom_polygon(data=corner_plot,aes(x,yp,fill=kind),colour=NA)+
 geom_line(data=null_plot,aes(x,yp,colour=kind),linewidth=.45)+
 geom_path(data=ode_plot,aes(x,yp,colour=kind),linewidth=.9)+
 geom_blank(data=limits,aes(x,yp))+
 facet_wrap(~panel,nrow=1,scales='free_y')+
 scale_x_continuous(breaks=seq(.4,1,.1),limits=c(.38,1.01),expand=c(0,0))+
 scale_y_continuous(breaks=function(v)if(min(v)< -1) -8:-1 else seq(0,.08,.02),
  labels=function(v)if(any(v<0))parse(text=paste0('10^',v))else format(v,nsmall=2),expand=c(0,0))+
 scale_colour_manual(values=c('ODE trajectory'='#EE3838','Susceptible nullcline'='#CE54C7'),
  breaks=c('ODE trajectory','Susceptible nullcline'),name=NULL)+
 scale_fill_manual(values=c('Boundary layer'='#FFF8D8','Corner region'='#FBE5D6'),breaks=c('Boundary layer','Corner region'),labels=c('Boundary layer: $y < y_{\\mathrm{BL}}$','Corner region'),name=NULL)+
 guides(colour=guide_legend(order=1,override.aes=list(alpha=1)),
        fill=guide_legend(order=2))+
 labs(x='x',y='y',
  subtitle=bquote(.(intToUtf8(0x211B))[0]==.(R0)~','~rho==.(rho)~','~theta==.(theta)~','~K==10^6~','~I[0]==.(I0)),
  caption=NULL)+
 theme_bw(base_size=10.5)+
 theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
  strip.background=element_rect(fill='grey93',colour='grey55'),strip.text=element_text(face='bold'),
  legend.position='bottom',legend.box='horizontal',legend.text=element_text(size=9),
  plot.margin=margin(5.5,7,5.5,5.5),plot.caption=element_text(hjust=0,size=8),
  plot.subtitle=element_text(family=if(capabilities('aqua'))'STIX'else 'serif',size=10),
  panel.spacing=grid::unit(1.1,'lines'))
# Export the maintained two-panel figure with native LaTeX typography.
source('validation/R/paper_graphics.R')
for (i in seq_along(p$scales$scales)) {
  sc <- p$scales$scales[[i]]
  if (is.expression(sc$labels)) p$scales$scales[[i]]$labels <- paper_text(sc$labels)
  if (is.function(sc$labels)) p$scales$scales[[i]]$labels <- local({
    f <- sc$labels; function(x) paper_text(f(x))
  })
}
p <- p + labs(subtitle='$\\mathcal{R}_0=1.5,\\quad \\rho=0.01,\\quad \\theta=0.5,\\quad K=10^6,\\quad I_0=3$') +
 theme(plot.subtitle=element_text(family=''))
tex <- file.path(out,'phaseplane_framework.tex')
paper_tikz_open(tex,9,4.9,pointsize=10)
grid::grid.draw(paper_grob(ggplot2::ggplotGrob(p)))
dev.off()
paper_compile(tex)
