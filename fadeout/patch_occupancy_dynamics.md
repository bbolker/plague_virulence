# Patch-Occupancy Dynamics

Let

$$
p(t)
$$

denote the proportion of patches that are infected and persistently occupied at time $t$.

The approximation is based on four assumptions:

1. The initial epidemic and first burnout occur rapidly relative to the later recolonization dynamics.
2. A patch that survives the first burnout persists thereafter.
3. Persistent patches rapidly approach the deterministic endemic equilibrium $I^*$.
4. Extinct patches rapidly recover to $S=K$.

---

## 1. Initial occupancy after the first burnout

Let

$$
P_1 = P_1(R_0,K)
$$

be the probability that an infected patch survives the first burnout.

If all patches are initially infected, then the initial epidemic is approximated as an instantaneous event. Immediately after the first burnout,

$$
p(0)=P_1.
$$

For fixed $R_0$ and $K$, $P_1$ is treated as a constant. It can be predicted using the emulator fitted to the single-patch simulation results.

---

## 2. Infection pressure from occupied patches

Assume that each occupied patch rapidly reaches the endemic equilibrium number of infected hosts,

$$
I^*=I^*(R_0,K).
$$

The fraction of occupied patches is $p(t)$, so the infection pressure on an empty patch is assumed to be proportional to

$$
p(t)I^*.
$$

Let $\alpha$ denote the between-patch transmission parameter. The colonization rate of an empty patch is approximated by

$$
c(t)=\alpha I^*p(t).
$$

This expression assumes that the between-patch transmission term is normalized so that the number of patches does not appear explicitly.

---

## 3. Successful colonization

A colonization event does not necessarily produce a persistent infected patch.

After infection is introduced into an empty patch, the new outbreak survives the first burnout with probability $P_1$. Therefore, the rate of successful colonization is

$$
c(t)P_1.
$$

The fraction of patches available for colonization is

$$
1-p(t).
$$

Hence the occupancy dynamics are

$$
\frac{dp}{dt}=(1-p(t))c(t)P_1.
$$

Substituting

$$
c(t)=\alpha I^*p(t)
$$

gives

$$
\frac{dp}{dt}=\alpha I^*P_1p(t)(1-p(t)).
$$

Define

$$
\lambda=\alpha I^*P_1.
$$

Then

$$
\frac{dp}{dt}=\lambda p(t)(1-p(t)).
$$

This is a logistic equation with carrying capacity equal to 1.

---

## 4. Analytical solution

Starting from

$$
\frac{dp}{dt}=\lambda p(1-p),
$$

separate variables:

$$
\frac{dp}{p(1-p)}=\lambda\,dt.
$$

Using

$$
\frac{1}{p(1-p)}=\frac{1}{p}+\frac{1}{1-p},
$$

we obtain

$$
\int\left(\frac{1}{p}+\frac{1}{1-p}\right)dp
=\int\lambda\,dt.
$$

Therefore,

$$
\log p-\log(1-p)=\lambda t+C,
$$

or

$$
\log\left(\frac{p}{1-p}\right)=\lambda t+C.
$$

Exponentiating,

$$
\frac{p}{1-p}=Ae^{\lambda t},
$$

where $A=e^C$. Solving for $p$,

$$
p(t)=\frac{1}{1+A^{-1}e^{-\lambda t}}.
$$

Using the initial condition

$$
p(0)=P_1,
$$

we find

$$
A^{-1}=\frac{1-P_1}{P_1}.
$$

Thus the predicted occupancy trajectory is

$$
p(t)=\frac{1}{1+\frac{1-P_1}{P_1}\exp\left(-\lambda t\right)}.
$$

Substituting

$$
\lambda=\alpha I^*P_1
$$

gives

$$
p(t)=\frac{1}{1+\frac{1-P_1}{P_1}\exp\left(-\alpha I^*P_1t\right)}.
$$

