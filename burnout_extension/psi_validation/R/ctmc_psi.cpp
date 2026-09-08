#include <Rcpp.h>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
IntegerVector simulate_counts_cpp(double R0,double rho,double theta,double psi,
                                  int K,int attempts,int seed,int max_events=20000000){
  int established=0,persistent=0,unresolved=0;
  #pragma omp parallel for reduction(+:established,persistent,unresolved) schedule(dynamic)
  for(int run=0;run<attempts;++run){
    std::mt19937_64 rng((uint64_t)(uint32_t)seed+0x9e3779b97f4a7c15ULL*(run+1));
    std::uniform_real_distribution<double> uniform(0.0,1.0);
    int S=K-1,I=1;bool peak=false,trough=false,driftpos=false,turned=false,done=false;
    for(int ev=0;ev<max_events;++ev){
      double gold=R0*((double)S/K)*std::pow((double)K/(S+I),psi)-1.0;
      double inf=R0*S*(double)I/K*std::pow((double)K/(S+I),psi);
      double rem=I,rec=S<K?K*rho*(1.0-(double)S/K)*std::pow((double)S/K,theta):0.0;
      double total=inf+rem+rec;if(total<=0)break;
      double u=uniform(rng)*total;
      if(u<inf){--S;++I;}else if(u<inf+rem)--I;else ++S;
      double gnew=I>0?R0*((double)S/K)*std::pow((double)K/(S+I),psi)-1.0:NA_REAL;
      if(!peak&&gold>0&&I>0&&gnew<=0)peak=true;
      if(peak&&!trough&&gold<0&&I>0&&gnew>=0)trough=true;
      if(trough){
        double drift=K*rho*(1.0-(double)S/K)*std::pow((double)S/K,theta)-
          R0*S*(double)I/K*std::pow((double)K/(S+I),psi);
        if(drift>0)driftpos=true;if(driftpos&&drift<0)turned=true;
      }
      if(I==0){if(peak)++established;done=true;break;}
      if(trough&&turned&&gold>0&&gnew<=0){++established;++persistent;done=true;break;}
    }
    if(!done)++unresolved;
  }
  return IntegerVector::create(_["attempts"]=attempts,_["established"]=established,
    _["persistent"]=persistent,_["unresolved"]=unresolved);
}
