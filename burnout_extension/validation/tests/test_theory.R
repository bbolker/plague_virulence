source('validation/R/theory.R')
for(th in c(0,.5,1)) for(x in c(.13,.37,.71)) {
  e<-1e-6; got<-action_diff(x-e,x+e,3,th)/(2*e)
  stopifnot(abs(got-action_derivative(x,3,th))<2e-7)
}
expected<-c('2'=0.07968023,'3'=0.14061515,'5'=0.16570580)
for(nm in names(expected)) stopifnot(abs(C_explicit(as.numeric(nm),.5)-expected[nm])<3e-6)
z<-matching(3,.01,.5,1e4,D=-.07962)
stopifnot(abs(action_diff(z$x_second,1/3,3,.5)/.01-z$Lambda_second)<2e-8)
cat('theory unit tests passed\n')
