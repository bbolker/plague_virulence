# Preview only: the same transition points are shown for all three methods.
library(data.table)
library(ggplot2)
d<-fread('validation/data/direct_entry_comparison_original_grid.csv',integer64='double')
d<-d[ODE_status=='OK' & is.finite(P_BI) & is.finite(P_first_K) &
       is.finite(P_second_K) & P_conditional>.01 & P_conditional<.99]
stopifnot(nrow(d)==42L)
methods<-c('P_BI','P_first_K','P_second_K')
labels<-c('BI','First-order entry + Kendall','Second-order entry + Kendall')
z<-melt(d,id.vars=c('R0','rho','theta','K','P_conditional','cond_low','cond_high'),
        measure.vars=methods,variable.name='method',value.name='prediction')
z[,method:=factor(method,levels=methods,labels=labels)]
z[,theta_label:=factor(theta,levels=seq(0,1,.25))]
mae<-z[,.(MAE=mean(abs(prediction-P_conditional))),by=method]
mae[,text:=sprintf('Mean absolute error = %.3f',MAE)]
pal<-c('#0072B2','#56B4E9','#009E73','#E69F00','#D55E00')
p<-ggplot(z,aes(P_conditional,prediction,colour=theta_label))+
 geom_abline(slope=1,intercept=0,colour='grey45',linetype='dashed',linewidth=.55)+
 geom_segment(aes(x=cond_low,xend=cond_high,y=prediction,yend=prediction),
              linewidth=.45,alpha=.65)+
 geom_point(size=2.3,alpha=.9)+
 geom_text(data=mae,aes(x=.04,y=.96,label=text),inherit.aes=FALSE,
           hjust=0,vjust=1,size=3.2,colour='grey20')+
 facet_wrap(~method,nrow=1)+coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+
 scale_x_continuous(breaks=seq(0,1,.25))+
 scale_y_continuous(breaks=seq(0,1,.25))+
 scale_colour_manual(values=pal,name=expression(theta))+
 labs(x='Simulated conditional persistence probability',
      y='Predicted conditional persistence probability',
      subtitle='Same 42 transition points in every panel; dashed line: prediction = simulation',
      caption='rho = 0.01; K = 10⁶, 10⁷, 10⁸, 10⁹. Simulation probabilities between 0.01 and 0.99. Horizontal bars: 95% Wilson intervals.')+
 theme_bw(base_size=11)+
 theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(colour='grey92',linewidth=.25),
       legend.position='bottom',strip.background=element_rect(fill='grey95'),
       plot.subtitle=element_text(size=10),plot.caption=element_text(hjust=0,size=8.5),
       panel.spacing=grid::unit(1.7,'lines'))
out<-'validation/figures/previews';dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(capabilities('aqua')){
 quartz(file=file.path(out,'entry_transition_comparison.pdf'),type='pdf',width=11,height=4.65)
 print(p);dev.off()
 png(file.path(out,'entry_transition_comparison.png'),width=2200,height=930,res=200,type='quartz')
 print(p);dev.off()
}else{
 ggsave(file.path(out,'entry_transition_comparison.pdf'),p,width=11,height=4.65,device=cairo_pdf)
 ggsave(file.path(out,'entry_transition_comparison.png'),p,width=11,height=4.65,dpi=200)
}
fwrite(d,file.path(out,'entry_transition_comparison_points.csv'))
print(mae)
