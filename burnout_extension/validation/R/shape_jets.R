# Taylor coefficients through fourth order; no finite differences.
.ord<-4L
jconst<-function(x)c(list(x),rep(list(0),.ord))
jvar<-function(x){a<-jconst(x);a[[2]]<-1;a}
jadd<-function(a,b)Map(`+`,a,b)
jscale<-function(a,s)lapply(a,function(x)x*s)
jsub<-function(a,b)jadd(a,jscale(b,-1))
jmul<-function(a,b)lapply(0:.ord,function(n)Reduce(`+`,lapply(0:n,function(k)a[[k+1]]*b[[n-k+1]])))
jinv<-function(a){b<-jconst(1/a[[1]]);for(n in 1:.ord)b[[n+1]]<--Reduce(`+`,lapply(1:n,function(k)a[[k+1]]*b[[n-k+1]]))/a[[1]];b}
jdiv<-function(a,b)jmul(a,jinv(b))
jexp<-function(a){b<-jconst(exp(a[[1]]));for(n in 1:.ord)b[[n+1]]<-Reduce(`+`,lapply(1:n,function(k)k*a[[k+1]]*b[[n-k+1]]))/n;b}
jexpm1<-function(a){b<-jexp(a);b[[1]]<-expm1(a[[1]]);b}
jlog<-function(a){iv<-jinv(a);b<-jconst(log(a[[1]]));for(n in 1:.ord)b[[n+1]]<-Reduce(`+`,lapply(1:n,function(k)k*a[[k+1]]*iv[[n-k+1]]))/n;b}
jpow<-function(a,p)jexp(jscale(jlog(a),p))
shape_derivatives4<-function(r,theta,nquad=256){
 quad<-statmod::gauss.quad(nquad,kind='legendre');tt<-(quad$nodes+1)/2;ww<-quad$weights/2
 R<-jvar(r);one<-jconst(1);z<-jconst(log(x_final(r)))
 for(it in 1:5){f<-jexp(z);fun<-jadd(jsub(one,f),jdiv(z,R));slope<-jsub(jinv(R),f);z<-jsub(z,jdiv(fun,slope))}
 f<-jexp(z);S<-jscale(z,-1);s<-jscale(S,tt);u<-jexp(jadd(z,s));F<-jadd(jsub(f,u),jdiv(s,R));dec<-jexp(jscale(s,theta-1))
 first<-jdiv(jmul(jmul(jsub(jinv(R),u),jsub(one,u)),dec),jmul(jsub(one,f),F))
 second<-jinv(jscale(jexpm1(jscale(s,-1)),-1));iv<-lapply(jsub(first,second),function(v)sum(ww*v));J<-jmul(S,iv)
 lc<-jadd(jsub(jadd(jlog(jsub(one,f)),jlog(jsub(jinv(R),f))),z),J)
 A2<-jsub(jdiv(jpow(R,theta-2),jsub(R,one)),jdiv(jmul(jpow(f,1-theta),jadd(jconst(1-theta),jscale(f,theta))),jsub(one,jmul(R,f))))
 c(logC=lc[[1]],C1=lc[[2]],C2=2*lc[[3]],C3=6*lc[[4]],C4=24*lc[[5]],A2=A2[[1]],A3=A2[[2]],A4=2*A2[[3]])
}
