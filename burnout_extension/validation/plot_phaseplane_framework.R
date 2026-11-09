# Framework preview, inspired by Parsons et al. (2024), Fig. 1D,E.
# Uses this manuscript's own threshold sqrt(y_star/K), not y_star.
source('validation/R/theory.R')
library(deSolve)
R0<-1.5;rho<-.01;theta<-.5;K<-1e6
xs<-1/R0;xf<-x_final(R0);ystar<-rho*h_theta(xs,theta);ybl<-sqrt(ystar/K)
y0<-1e-10;x0<-1-R0/(R0-1+rho)*y0
rhs<-function(t,z,p)list(c(rho*h_theta(z[1],theta)-R0*z[1]*exp(z[2]),R0*z[1]-1))
z<-as.data.frame(ode(c(x=x0,logy=log(y0)),seq(0,1200,by=.025),rhs,NULL,
 method='lsoda',rtol=2e-11,atol=c(1e-13,1e-12),maxsteps=1e6))
z$y<-exp(z$logy);n<-nrow(z)
peak<-which(z$x[-n]>=xs & z$x[-1]<xs)[1]
trough<-which(seq_len(n-1)>peak & z$x[-n]<xs & z$x[-1]>=xs)[1]
entry<-which(seq_len(n-1)>peak & z$y[-n]>=ybl & z$y[-1]<ybl)[1]
exit<-which(seq_len(n-1)>trough & z$y[-n]<ybl & z$y[-1]>=ybl)[1]
stopifnot(all(is.finite(c(peak,trough,entry,exit))))
g<-rho*h_theta(z$x,theta)-R0*z$x*z$y
turn<-which(seq_len(n-1)>peak & g[-n]<0 & g[-1]>=0)[1]
w<-(ybl-z$y[entry])/(z$y[entry+1]-z$y[entry]);xin<-z$x[entry]+w*(z$x[entry+1]-z$x[entry])
end<-which(seq_len(n)>exit & z$y>3*ybl)[1]
tr<-z[seq_len(end),]
corner<-z[turn,];ymin<-min(z$y[peak:trough]);xmin<-min(tr$x)
cat(sprintf('R0=%g rho=%g theta=%g K=%g\nx_f=%.8g x_star=%.8g y_star=%.8g y_BL=%.8g x_in=%.8g\ncorner x=%.8g y=%.8g; y_min=%.8g\n',R0,rho,theta,K,xf,xs,ystar,ybl,xin,corner$x,corner$y,ymin))
out<-'validation/figures/previews';dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(data.frame(R0,rho,theta,K,xf,xs,ystar,ybl,xin,xturn=corner$x,yturn=corner$y,ymin),
 'validation/data/phaseplane_framework_parameters.csv',row.names=FALSE)
write.csv(tr,'validation/data/phaseplane_framework_trajectory.csv',row.names=FALSE)
if(capabilities('aqua')) quartzFonts(STIX=quartzFont(c('STIXGeneral-Regular','STIXGeneral-Bold','STIXGeneral-Italic','STIXGeneral-BoldItalic')))
blue<-'#0072B2';orange<-'#D55E00';gold<-'#F0C44B';cyan<-'#30A3A8';gray<-'#777777'
arrow_on<-function(j,logscale=FALSE,colour=blue){
 k<-min(j+30,nrow(tr));if(j<1||k<=j)return()
 arrows(tr$x[j],tr$y[j],tr$x[k],tr$y[k],length=.075,angle=22,col=colour,lwd=1.6)
}
plot_it<-function(){
 par(mfrow=c(1,2),mar=c(4.1,4.6,2.3,1),oma=c(2.7,0,0,0),mgp=c(2.7,.7,0),las=1,
     cex.axis=.85,cex.lab=1, xaxs='i',yaxs='i',family='serif')
 xx<-seq(.36,1,length.out=1800);ff<-F0(xx,R0);nn<-rho*h_theta(xx,theta)/(R0*xx)
 # Full first-wave geometry.
 plot(NA,xlim=c(.36,1.02),ylim=c(0,.085),xlab='Susceptible fraction, x',ylab='Infectious fraction, y')
 rect(.36,0,1.02,ybl,col='#FFF3C7',border=NA)
 abline(v=xs,col=cyan,lty=3,lwd=1.1)
 lines(xx,nn,col='#B571B0',lwd=1.1)
 lines(xx,ifelse(ff>0,ff,NA),col=gray,lty=2,lwd=1.3)
 lines(tr$x,tr$y,col=blue,lwd=1.8)
 points(xin,ybl,pch=21,bg=gold,col='#8C6D1F',cex=1.15)
 points(corner$x,corner$y,pch=21,bg=orange,col=orange,cex=.9)
 arrow_on(which.min(abs(tr$time-(tr$time[peak]-4))))
 arrow_on(which.min(abs(tr$time-(tr$time[peak]+7))))
 mtext('(a) First epidemic wave',side=3,line=.6,adj=0,font=2,cex=.95)
 legend('topright',legend=c('Deterministic trajectory','Leading-order trajectory F(x)',
  'Susceptible nullcline','x = x*','Minimum of x'),
  col=c(blue,gray,'#B571B0',cyan,orange),lty=c(1,2,1,3,NA),
  pch=c(NA,NA,NA,NA,21),pt.bg=c(NA,NA,NA,NA,orange),
  lwd=c(1.8,1.3,1.1,1.1,1),bty='n',cex=.77)
 text(.46,.035,'Outer region',cex=.9,adj=0)
 jo<-peak-1L+which.min(abs(tr$x[peak:turn]-.55))
 arrows(.56,.037,tr$x[jo],tr$y[jo],length=.065,col='grey35',lwd=.8)
 # Dashed ellipse is a visual locator, not a boundary of asymptotic validity.
 phi<-seq(0,2*pi,length.out=200)
 lines(corner$x+.035*cos(phi),.008+.0075*sin(phi),col=orange,lty=3,lwd=.9)
 text(.38,.023,'Corner region',cex=.88,adj=0,col=orange)
 arrows(.43,.020,.43,.0145,length=.065,col=orange,lwd=.9)
 # Horizontal positions marked below the axis, not used as hard region boundaries.
 axis(1,at=c(xf,xs),labels=FALSE,tck=-.018)
 mtext(expression(x[f]),side=1,line=1.35,at=xf,cex=.85,col=gray)
 mtext(expression(x['*']),side=1,line=1.35,at=xs,cex=.85,col=cyan)
 # Same trajectory with log prevalence to resolve the stochastic region.
 yminplot<-10^floor(log10(ymin/2));ymaxplot<-.1
 plot(NA,xlim=c(.36,1.02),ylim=c(yminplot,ymaxplot),log='y',
      xlab='Susceptible fraction, x',ylab='Infectious fraction, y (log scale)',yaxt='n')
 ticks<-10^seq(ceiling(log10(yminplot)),floor(log10(ymaxplot)),by=1)
 axis(2,at=ticks,labels=as.expression(lapply(log10(ticks),function(k)bquote(10^.(k)))),las=1)
 rect(.36,yminplot,1.02,ybl,col='#FFF3C7',border=NA)
 abline(v=xs,col=cyan,lty=3,lwd=1.1)
 abline(h=ybl,col='#B58B20',lty=2,lwd=1.2)
 lines(xx,nn,col='#B571B0',lwd=1.1)
 lines(tr$x,tr$y,col=blue,lwd=1.8)
 points(xin,ybl,pch=21,bg=gold,col='#8C6D1F',cex=1.45)
 points(corner$x,corner$y,pch=21,bg=orange,col=orange,cex=.95)
 points(xs,z$y[trough],pch=21,bg='white',col=blue,cex=.9)
 arrow_on(which.min(abs(tr$time-(tr$time[peak]+13))))
 arrow_on(which.min(abs(tr$time-(tr$time[trough]-20))))
 arrow_on(which.min(abs(tr$time-(tr$time[trough]+35))))
 mtext('(b) Entry into the low-prevalence region',side=3,line=.6,adj=0,font=2,cex=.95)
 text(.99,ybl*1.7,expression(y[BL]==sqrt(y['*']/K)),adj=1,cex=.87,col='#8C6D1F')
 text(.41,ybl/5,'Stochastic boundary layer',adj=0,cex=.88,col='#775A16')
 text(.41,ybl/12,expression(y<y[BL]),adj=0,cex=.88,col='#775A16')
 text(.40,ybl*8,expression(paste('Entry: (',x["in"],', ',y[BL],')')),adj=0,cex=.86)
 arrows(.46,ybl*5,xin,ybl*1.08,length=.065,col='#8C6D1F',lwd=.9)
 lines(corner$x+.035*cos(phi),10^(log10(corner$y)+.6*sin(phi)),col=orange,lty=3,lwd=.9)
 text(.73,.012,'Corner region',adj=0,cex=.86,col=orange)
 arrows(.72,.010,.475,.012,length=.065,col=orange,lwd=.9)
 text(.79,ymin*2.4,'First trough',adj=0,cex=.85,col=blue)
 arrows(.78,ymin*2.0,xs+.005,z$y[trough],length=.065,col=blue,lwd=.8)
 mtext(bquote(.(intToUtf8(0x211B))[0]==.(R0)~','~rho==.(rho)~','~theta==.(theta)~','~K==10^6),side=1,outer=TRUE,line=1.1,cex=.95,family=if(capabilities('aqua'))'STIX' else 'serif')
}
if(capabilities('aqua')){
 quartz(file=file.path(out,'phaseplane_framework_preview.pdf'),type='pdf',width=11.5,height=5.4)
 plot_it();dev.off()
 png(file.path(out,'phaseplane_framework_preview.png'),type='quartz',width=2300,height=1080,res=200)
 plot_it();dev.off()
}else{
 pdf(file.path(out,'phaseplane_framework_preview.pdf'),width=11.5,height=5.4);plot_it();dev.off()
 png(file.path(out,'phaseplane_framework_preview.png'),width=2300,height=1080,res=200);plot_it();dev.off()
}
