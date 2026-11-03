source('validation/R/stochastic.R')
check <- function(S,I,established,persist,rho=.2) {
  ans<-classify_first_trough(S,I,2,rho,0,100)
  stopifnot(identical(unname(ans['established']),established),
            identical(unname(ans['persist']),persist))
}
# A genuine recovery is not undone by extinction in a later trough.
check(c(99,49,51,55,49),c(1,20,1,15,15),TRUE,TRUE)
check(c(99,49,51,55,49,49),c(1,20,1,15,15,0),TRUE,TRUE)
# Extinction before recovery, and early fizzle, remain failures.
check(c(99,49,51,51),c(1,20,1,0),TRUE,FALSE)
check(c(99,99),c(1,0),FALSE,FALSE)
# Time/event limits without extinction or completed recovery are unresolved.
check(c(99,49,51),c(1,20,1),TRUE,NA)
check(c(99,98),c(1,2),FALSE,NA)
# A brief second crossing without the required susceptible turn is insufficient.
check(c(99,49,51,49),c(1,20,1,1),TRUE,NA)
# At the exact threshold, use the same crossing convention as one_ctmc().
check(c(99,50,49,50,55,50),c(1,20,20,1,15,15),TRUE,TRUE)
# Exercise the tau wrapper with supplied output; do not simulate.
f<-one_adaptive_tau;bb<-body(f)
ix<-which(vapply(as.list(bb),function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name('z')),logical(1)))
stopifnot(length(ix)==1L)
bb[[ix]]<-quote(z<-trajectory);body(f)<-bb
e<-new.env(parent=globalenv());e$requireNamespace<-function(...)TRUE
e$trajectory<-cbind(S=c(99,49,51,55,49,49),I=c(1,20,1,15,15,0))
environment(f)<-e
stopifnot(isTRUE(f(2,.2,0,100)['persist']))
cat('stochastic classification checks passed (no stochastic simulations run)\n')
