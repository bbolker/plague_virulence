source('validation/R/theory.R');source('validation/R/extrema.R');source('validation/R/shape_jets.R');library(data.table)
rr<-c(seq(1.1,12,length.out=100),exp(seq(log(12),log(300),length.out=80)));ths<-seq(0,1,.025);out<-list();j<-0
for(th in ths){for(r in rr){z<-shape_derivatives4(r,th);j<-j+1
 cv<-unname(z['C2']+2*z['C1']/r)
 jv<-cv-1/(2*(r-1)^2)+th/(2*r^2)+1/(r*(r-1))-th/r^2
 out[[j]]<-data.table(R0=r,theta=th,C1=unname(z['C1']),C_x_curvature=cv,b_x_curvature=jv)
};cat(th,'done\n')}
x<-rbindlist(out);fwrite(x,'validation/data/inverse_R_logconcavity_checks.csv');print(x[order(-C_x_curvature)][1:10]);cat('C positive',nrow(x[C_x_curvature>=0]),'b positive',nrow(x[b_x_curvature>=0]),'of',nrow(x),'\n')
