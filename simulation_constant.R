library("burnout")
# Initialize parameters
T <- 100                    # Number of time steps
P <- numeric(T)             # Proportion of patches with infected rats
NI <- numeric(T)              # Average rat population per infected patch
NS <- numeric(T)             # Average rat population per susceptible patch

# Initial conditions
P[1] = 0.01                  
NI[1] = 100
NS[1] = 100

# Model parameters
B0 = 0.8                     # Infection rate between patches
r = 0.5                     # growth rate of rat population
D = 0.5                    # Death rate due to pathogen
epsilon=0.05
R0=10
ki=0.01
ks=0.01
c0=0.2
K=300               # Total number of patches
N = 1000                # Total rat population

# Simulation loop
for (t in 1:(T-1)) {
  # Current state
  p = P[t]
  ni = NI[t]
  ns = NS[t]
  
  
  if(p<1/N){
    p=0
    ni=0
  }else{
  # migration and infection
  B=B0*(1-exp(-ki*ni))*(1-exp(-ks*ns))
  delta_p=B*p*(1-p)
  ni=(p*ni+delta_p*ns)/(p+delta_p)
  p = p + delta_p 

  c=min(c0,1/(1/p+1/(1-p)))
  delta_N <- c*(ns-ni)
  ni <- ni + delta_N/p
  ns <- ns - delta_N/(1-p)
  
  #fizzle and burnout
  P1=P1_prob(R0,epsilon,k=1,N=ni)
  delta_p<-(1-P1)*p
  ns<-((1-p)*ns+delta_p*ni)/(1-p+delta_p)
  p <-p-delta_p
  
  #epidemic
  z=final_size(R0)
  ni <- ni * (1 - z*D)  
  }
  
  # birth
  ns=ns+r*ns*(1-ns/K)
  ni=ni+r*ni*(1-ni/K)
 
   NI[t+1]=ni
   NS[t+1]=ns
   P[t+1]=p

}

plot(1:T,P[1:T],type="l",main="patch state",xlab="time step",ylab="infected proportion")
plot(1:T,NI[1:T],type="l",main="average population size of infected patches",xlab="time step",ylab="NI")
plot(1:T,NS[1:T],type="l",main="average population size of susceptible patches",xlab="time step",ylab="NS")
