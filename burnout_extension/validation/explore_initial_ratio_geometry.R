source('validation/R/theory.R');source('validation/R/extrema.R');source('validation/R/shape_jets.R');library(data.table)
rr<-seq(1.1,4,length.out=450);ths<-seq(0,1,.025);out<-list();j<-0
for(th in ths){for(r in rr){
 z<-shape_derivatives4(r,th);a1<-unname(action_barrier_derivatives(r,th)['d1']);q<-r*(r-1);qp<-2*r-1
 c0<-z['C1']+1/(2*(r-1))-th/(2*r);c1<-z['C2']-1/(2*(r-1)^2)+th/(2*r^2);c2<-z['C3']+1/(r-1)^3-th/r^3;c3<-z['C4']-3/(r-1)^4+3*th/r^4
 U<--q*c0;Up<--qp*c0-q*c1;Upp<--2*c0-2*qp*c1-q*c2;Uppp<--6*c1-3*qp*c2-q*c3
 V<-q*a1;Vp<-qp*a1+q*z['A2'];Vpp<-2*a1+2*qp*z['A2']+q*z['A3'];Vppp<-6*z['A2']+3*qp*z['A3']+q*z['A4']
 U1<-q*Up;V1<-q*Vp;U2<-q*(qp*Up+q*Upp);V2<-q*(qp*Vp+q*Vpp);U3<-q*(qp^2+2*q)*Up+3*q^2*qp*Upp+q^3*Uppp;V3<-q*(qp^2+2*q)*Vp+3*q^2*qp*Vpp+q^3*Vppp
 if(U1<0 && V1>0){j<-j+1;out[[j]]<-data.table(R0=r,theta=th,U=as.numeric(U),V=as.numeric(V),U1=as.numeric(U1),V1=as.numeric(V1),W=as.numeric(U1*V2-U2*V1),Wt=as.numeric(U1*V3-U3*V1),M=as.numeric(V*U1+(1-U)*V1),Mt=as.numeric(V*U2+(1-U)*V2))}
 };cat('theta',th,'done\n')}
z<-rbindlist(out);fwrite(z,'validation/data/initial_ratio_geometry.csv');roots<-list();k<-0
for(th in unique(z$theta)){a<-z[theta==th][order(R0)];for(name in c('W','M')){ix<-which(diff(sign(a[[name]]))!=0);for(i in ix){fraction<--a[[name]][i]/(a[[name]][i+1]-a[[name]][i]);k<-k+1;roots[[k]]<-data.table(theta=th,name=name,R0=a$R0[i]+fraction*(a$R0[i+1]-a$R0[i]),slope=a[[paste0(name,'t')]][i]+fraction*(a[[paste0(name,'t')]][i+1]-a[[paste0(name,'t')]][i]))}}}
x<-rbindlist(roots);print(x);fwrite(x,'validation/data/initial_ratio_zero_orientations.csv')
