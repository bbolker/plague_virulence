# Complete critical-K classification conditional on two model inequalities

This theorem closes the abstract root-count/partition step. Its two global
hypotheses have NOT yet been proved for the full burnout amplitude on the
entire parameter domain, so it is not a completed theorem about the paper's
BI probability.

Let R>1, t=log(1-1/R), B=K b(R), s=-d_t log b, and
H(B)=(exp B-1)/B. Put D(B)=B H'(B)>0. The derivative of the unconditional
probability has the sign of E=H(B)-s. On s<=1 it is positive. On the domain
D_s={R:s>1}, define

Kc(R)=H^{-1}(s(R))/b(R).

There E<0 iff K<Kc, and E>0 iff K>Kc.
Assume the established endpoint signs (positive near R=1 and negative
sufficiently far out), s<1 near R=1, s->infinity, and Kc->infinity as R->infinity.

The two additional hypotheses are:
1. I={t:s>1,s_t<0} is empty or one interval.
2. Q=s_t s_ttt-s_tt^2<=0 on I.

## At most one decreasing segment of Kc

On s_t>=0, d_t log Kc=s+s_t/D>0. On I set h=-s_t>0 and Y(s)=h(t(s)).
The second hypothesis gives Y''=(log h)''/h<=0. Define
F(s)=s D(H^{-1}(s)). Then

sign(d_t Kc)=sign(F(s)-Y(s)).

The previously proved F''(s)>0 implies (Y-F)''<0. Its positive set is thus
an interval. Consequently Kc has at most one decreasing segment; its ends,
if interior, are a single local maximum followed by a single local minimum.

The domain D_s has at most two components. Every component that ends at a
finite R must contain a terminal segment with s_t<0 and s>1. Two such finite
components would give two separated parts of I. Thus there is at most one
finite component followed by the unbounded final component. Boundary
contacts with s=1 do not cause probability stationary points because they
have Kc=0, while K>0.

## Parameter partition

If D_s is connected and Kc has a strict local maximum K_H and subsequent
strict local minimum K_L, then 0<K_L<K_H and:

- K_L<K<K_H: exactly three simple stationary points, ordered maximum,
  minimum, maximum;
- K<K_L or K>K_H: exactly one simple maximum;
- K=K_L or K=K_H: one simple maximum plus one stationary inflection point
  (the probability derivative has a double zero).

If D_s has two components, the bounded first component has a unique maximum
K_H and endpoint values zero. The final component is strictly increasing
from zero to infinity. The same classification holds with K_L=0. In
particular, 0<K<K_H gives three ordinary extrema, and K>K_H gives one.

If Kc is strictly increasing with no degenerate contact, there is one simple
maximum for every K>0. The remaining boundary case is a collapse of the
maximum/minimum pair of Kc: Y-F has maximum zero at one point. Its second
derivative is strictly negative, so Kc' has a zero of exactly order two.
At K=Kc there, E has a zero of exactly order three; the probability has a
strict quartic maximum. This excludes higher-order degeneracy under these
hypotheses. For other K there is one ordinary maximum.

These statements follow by counting horizontal-line intersections with Kc.
The sign rule determines the order and type of each probability stationary
point. A strict concavity contact can only have order two, establishing the
stated degeneracy orders.

Thus the hypotheses imply a complete partition via the surfaces K=K_L and
K=K_H, including degenerate boundary cases. The dependence of these surfaces
on rho need not be monotone: the three certified local folds at fixed
theta=.55,K=exp(15/4) illustrate why the partition cannot be replaced by one
rho threshold. Proving hypotheses 1 and 2 for the exact U,V remains the main
unfinished task.
