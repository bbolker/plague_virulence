/* Explicit Poisson tau-leap with a small-population event budget.
 * Full SSA is selected when a leap would contain fewer than eight events.
 * At small S or I the expected number of events affecting that population
 * is capped separately, to control hidden zero-hitting/replenishment order.
 * Model rates and first-trough classifier match stochastic.R.
 */
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>

void guarded_counts(double *p,int *n,double *out) {
  double R=p[0],rho=p[1],theta=p[2],K=p[3],eps=p[4];
  double budget=p[5],small=p[6],tf=p[7],xs=1/R;
  GetRNGstate();
  for(int k=0;k<*n;k++) {
    double S=K-1,I=1,t=0;
    int peak=0,trough=0,positive=0,turned=0,zero=0,done=0;
    for(long step=0;step<20000000;step++) {
      if(step%10000==0)R_CheckUserInterrupt();
      if(I==0 || (theta>0 && S==0)){zero|=(S==0);done=1;break;}
      double old=S,inf=R*S*I/K,rem=I;
      double rec=S<K?K*rho*(1-S/K)*pow(S/K,theta):0;
      double total=inf+rem+rec;
      if(total<=0 || t>=tf)break;
      double bs=fmax(eps*S,1),bi=fmax(eps*I,1),tau=tf-t;
      if(rec!=inf)tau=fmin(tau,bs/fabs(rec-inf));
      if(inf+rec>0)tau=fmin(tau,bs*bs/(inf+rec));
      if(inf!=rem)tau=fmin(tau,bi/fabs(inf-rem));
      if(inf+rem>0)tau=fmin(tau,bi*bi/(inf+rem));
      if(theta>0 && S<=small && inf+rec>0)tau=fmin(tau,budget/(inf+rec));
      if(I<=small && inf+rem>0)tau=fmin(tau,budget/(inf+rem));
      if(tau*total<8) {
        double dt=exp_rand()/total;
        if(t+dt>tf)break;
        t+=dt;double u=unif_rand()*total;
        if(u<inf){S--;I++;}else if(u<inf+rem){I--;}else{S++;}
        out[5]++;out[7]++;
      } else {
        int accepted=0;
        for(int retry=0;retry<50;retry++) {
          double a=rpois(inf*tau),b=rpois(rem*tau),c=rpois(rec*tau);
          double nextS=S-a+c,nextI=I+a-b;
          if(nextS>=0 && nextS<=K && nextI>=0) {
            S=nextS;I=nextI;t+=tau;out[5]+=a+b+c;out[6]++;
            accepted=1;break;
          }
          tau/=2;out[8]++;
        }
        if(!accepted)error("Failed to find a nonnegative tau leap");
      }
      if(S==0)zero=1;
      if(!peak && old/K>xs && S/K<=xs && I>0)peak=1;
      if(peak && !trough && old/K<xs && S/K>=xs && I>0)trough=1;
      if(trough) {
        double g=K*rho*(1-S/K)*pow(S/K,theta)-R*S*I/K;
        if(g>0)positive=1;
        if(positive && g<0)turned=1;
      }
      if(I==0 || (theta>0 && S==0)){done=1;break;}
      if(trough && turned && old/K>xs && S/K<=xs && I>0){out[1]++;done=1;break;}
    }
    out[0]+=peak;out[2]+=!done;out[3]+=zero;out[4]+=trough;
  }
  PutRNGstate();
}
