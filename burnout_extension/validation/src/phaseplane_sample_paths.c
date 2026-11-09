#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>
/* Exact Gillespie sample paths observed on a regular plotting grid. */
void phaseplane_paths(double *p,int *npaths,int *nobs,double *dt,double *sx,double *iy){
 double R=p[0],rho=p[1],th=p[2],K=p[3],I0=p[4];
 GetRNGstate();
 for(int j=0;j<*npaths;j++){
  double S=K-I0,I=I0,t=0,next=0;int event=0;
  for(int m=0;m<*nobs;m++){
   double target=m*(*dt);
   if(m==0){sx[j*(*nobs)+m]=S/K;iy[j*(*nobs)+m]=I/K;}
   else{
    while(1){
     if(next==0){
      double inf=R*S*I/K,rem=I,rec=S<K?K*rho*(1-S/K)*pow(S/K,th):0;
      double total=inf+rem+rec;
      if(total<=0){next=INFINITY;break;}
      next=t+exp_rand()/total;double u=unif_rand()*total;
      event=u<inf?0:(u<inf+rem?1:2);
     }
     if(next>target)break;
     t=next;next=0;
     if(event==0){S--;I++;}else if(event==1){I--;}else S++;
    }
    sx[j*(*nobs)+m]=S/K;iy[j*(*nobs)+m]=I/K;
   }
   if(m%500==0)R_CheckUserInterrupt();
  }
 }
 PutRNGstate();
}
