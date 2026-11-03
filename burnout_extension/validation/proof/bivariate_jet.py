"""Formal bivariate Taylor arithmetic with outward-rounded Arb coefficients.
Coefficients are ordinary Taylor coefficients. This module has no remainder
semantics: callers must use separate analytic/Cauchy remainder bounds.
"""
from flint import arb

class BiJet:
    order=16
    theta_order=12
    @classmethod
    def keys(cls):return [(n-j,j) for n in range(cls.order+1) for j in range(min(n,cls.theta_order)+1)]
    def __init__(self,data=0):
        if isinstance(data,dict):self.c={ij:arb(v) for ij,v in data.items() if sum(ij)<=self.order and ij[1]<=self.theta_order and v!=0}
        else:self.c={} if data==0 else {(0,0):arb(data)}
    def at(self,i,j):return self.c.get((i,j),arb(0))
    def __add__(self,b):
        if not isinstance(b,BiJet):b=BiJet(b)
        c=self.c.copy()
        for ij,v in b.c.items():c[ij]=c.get(ij,arb(0))+v
        return BiJet(c)
    __radd__=__add__
    def __neg__(self):return BiJet({ij:-v for ij,v in self.c.items()})
    def __sub__(self,b):return self+-b if isinstance(b,BiJet) else self+(-arb(b))
    def __rsub__(self,b):return -self+b
    def __mul__(self,b):
        if not isinstance(b,BiJet):return BiJet({ij:v*b for ij,v in self.c.items()})
        c={}
        for (i,j),u in self.c.items():
            for (k,l),v in b.c.items():
                if i+j+k+l<=self.order and j+l<=self.theta_order:
                    ij=(i+k,j+l);c[ij]=c.get(ij,arb(0))+u*v
        return BiJet(c)
    __rmul__=__mul__
    def inv(self):
        a0=self.at(0,0)
        if a0.contains(0):raise ZeroDivisionError('constant jet contains zero')
        out={(0,0):1/a0};nonconstant=[(ij,v) for ij,v in self.c.items() if ij!=(0,0)]
        for i,j in self.keys()[1:]:
            val=sum(v*out.get((i-k,j-l),arb(0)) for (k,l),v in nonconstant if k<=i and l<=j)
            out[(i,j)]=-val/a0
        return BiJet(out)
    def __truediv__(self,b):return self*b.inv() if isinstance(b,BiJet) else self*(1/arb(b))
    def __rtruediv__(self,b):return self.inv()*b
    def exp(self):
        out={(0,0):self.at(0,0).exp()};parts=[(ij,v) for ij,v in self.c.items() if ij!=(0,0)]
        for i,j in self.keys()[1:]:
            out[(i,j)]=sum((k+l)*v*out.get((i-k,j-l),arb(0)) for (k,l),v in parts if k<=i and l<=j)/(i+j)
        return BiJet(out)
    def log(self):
        euler=BiJet({ij:sum(ij)*v for ij,v in self.c.items() if ij!=(0,0)})
        p=euler*self.inv()
        return BiJet({(0,0):self.at(0,0).log(),**{ij:v/sum(ij) for ij,v in p.c.items() if ij!=(0,0)}})
    def Ds(self):return BiJet({(i-1,j):i*v for (i,j),v in self.c.items() if i>0})
    def truncate(self,n):return {ij:v for ij,v in self.c.items() if sum(ij)<=n}
    def evaluate(self,ds,dt,n):
        value=arb(0)
        for j in range(n,-1,-1):
            part=arb(0)
            for i in range(n-j,-1,-1):part=part*ds+self.at(i,j)
            value=value*dt+part
        return value
