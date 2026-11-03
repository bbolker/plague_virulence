source('validation/R/theory.R');source('validation/R/extrema.R');library(data.table)
expanded<-function(r,th){f<-x_final(r);delta<-1-r*f;aa<-(1-f)*f^(th-1);df<-f*(f-1)/delta;ga<-(1-th+th*f)/delta
 v<-function(u){F<-1-u+log(u)/r;(1-u)*u^th/(aa*u^2*F)*((1/r-u)*(log(u)/(r^2*F)-ga)-1/r^2)-df/(u-f)^2}
 ee<-min(1e-4*(1-f),.003*f);ss<-ee*c(1,1.4,2,2.8,4);coef<-lm.fit(cbind(1,ss,ss^2),v(f+ss))$coefficients
 endpoint<-coef[1]*ee+coef[2]*ee^2/2+coef[3]*ee^3/3
 -(1-th)*(1-f)/delta-1/(r*delta)+(1-f)/(2*delta^2)+endpoint+integrate(v,f+ee,1,rel.tol=1e-7,abs.tol=1e-7,subdivisions=1000)$value}
z<-CJ(R0=c(1.2,2,3,5,10),theta=c(0,.25,.5,.75,1));z[,expanded:=mapply(expanded,R0,theta)];z[,automatic:=mapply(function(r,t)unname(C_log_derivatives(r,t)['d1']),R0,theta)];z[,error:=expanded-automatic];print(z);fwrite(z,'validation/data/expanded_C_derivative_checks.csv');stopifnot(max(abs(z$error))<1e-5)
