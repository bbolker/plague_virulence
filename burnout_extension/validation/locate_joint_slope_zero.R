source('validation/R/theory.R');source('validation/R/extrema.R');source('validation/R/shape_jets.R');library(data.table)
F<-function(p){r<-p[1];th<-p[2];cd<-C_log_derivatives(r,th);ad<-action_barrier_derivatives(r,th);q<-r*(r-1);qp<-2*r-1;c0<-cd['d1']+.5/(r-1)-th/(2*r);c1<-cd['d2']-.5/(r-1)^2+th/(2*r^2);c(U1=unname(-q*(qp*c0+q*c1)),V1=unname(q*(qp*ad['d1']+q*ad['d2'])))}
p<-c(1.83,.12)
for(it in 1:15){h<-c(1e-5,1e-5);J<-sapply(1:2,function(j){a<-b<-p;a[j]<-a[j]+h[j];b[j]<-b[j]-h[j];(F(a)-F(b))/(2*h[j])});step<-solve(J,F(p));p<-p-step;if(max(abs(step))<1e-11)break}
z<-data.table(R0=p[1],theta=p[2],S=-log(x_final(p[1])),U1=F(p)[1],V1=F(p)[2]);print(z);fwrite(z,'validation/data/joint_slope_zero_center.csv')
