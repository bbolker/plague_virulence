
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>
void exact_counts(double *p,int *n,double *out) {
  double R=p[0],rho=p[1],theta=p[2],K=p[3],xs=1/R;
  GetRNGstate();
  for(int k=0;k<*n;k++) {
    double S=K-1,I=1;
    int peak=0,trough=0,positive=0,turned=0,zero=0,done=0;
    for(long ev=0;ev<20000000;ev++) {
      if(ev%100000==0) R_CheckUserInterrupt();
      double old=S,inf=R*S*I/K,rem=I;
      double rec=S<K?K*rho*(1-S/K)*pow(S/K,theta):0;
      double total=inf+rem+rec;
      if(total<=0) break;
      double u=unif_rand()*total;
      if(u<inf){S--;I++;}else if(u<inf+rem){I--;}else{S++;}
      if(S==0)zero=1;
      if(!peak && old/K>xs && S/K<=xs && I>0)peak=1;
      if(peak && !trough && old/K<xs && S/K>=xs && I>0)trough=1;
      if(trough){
        double g=K*rho*(1-S/K)*pow(S/K,theta)-R*S*I/K;
        if(g>0)positive=1;
        if(positive && g<0)turned=1;
      }
      out[5]++;
      if(I==0){done=1;break;}
      if(trough && turned && old/K>xs && S/K<=xs && I>0){out[1]++;done=1;break;}
    }
    out[0]+=peak;out[2]+=!done;out[3]+=zero;out[4]+=trough;
  }
  PutRNGstate();
}