library(data.table)
bas<-readRDS('validation/data/broad_structure_base.rds');rhos<-10^seq(-6,3,length.out=301)
bb<-exp(seq(log(1e-10),log(80),length.out=50000));logH<-ifelse(bb<1e-4,log1p(bb/2+bb^2/6),log(expm1(bb))-log(bb))
rows<-list();j<-0
for(z in bas){
 z<-z[y<=log(1e7)];y<-z$y;r<-exp(y);q<-r*(r-1);th<-z$theta[1]
 for(rv in rhos){
 Lr<-z$rC+.5*r/(r-1)-th/2-z$rA/rv;s<--(r-1)*Lr
 ix<-which(s>1.000001);if(length(ix)<3)next
 bc<-exp(approx(logH,log(bb),xout=log(s[ix]),rule=2)$y)
 logb<-z$logC+.5*(log(rv)+log(r-1)-th*y-log(2*pi))-z$A/rv
 kcrit<-log(bc)-logb[ix]
 pieces<-split(seq_along(ix),cumsum(c(TRUE,diff(ix)>1)))
 peaks<-0;troughs<-0
 for(ind in pieces){if(length(ind)<6)next
  v<-diff(kcrit[ind]);sg<-sign(v);ch<-which(diff(sg)!=0)
  # Only count stationary changes that persist over more than one mesh interval.
  for(k in ch){if(k<2||k>length(v)-2)next;if(all(sg[(k-1):k]>0)&&all(sg[(k+1):(k+2)]<0))peaks<-peaks+1;if(all(sg[(k-1):k]<0)&&all(sg[(k+1):(k+2)]>0))troughs<-troughs+1}
 }
 j<-j+1;rows[[j]]<-data.table(rho=rv,theta=th,components=length(pieces),Kcurve_peaks=peaks,Kcurve_troughs=troughs)
 }
}
a<-rbindlist(rows);fwrite(a,'validation/data/critical_K_geometry_exploration.csv');print(a[,.N,by=.(components,Kcurve_peaks,Kcurve_troughs)][order(components,Kcurve_peaks)]);print(a[Kcurve_peaks>1|Kcurve_troughs>1][1:15])
