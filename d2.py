import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fast import *
import c2
mp.mp.dps=40
def tau1_star(name,cf,x0,av):
    gf,Gf,_=FUN[name]
    def bvec(x):
        g=gf(x,av); G=Gf(x,av); c=cf(x)
        dc=[mp.diff(lambda z: cf(z)[m],x) for m in R3]
        return [c[1]*dc[m]+sum(G[m][n][l]*c[n]*c[l] for n in R3 for l in R3) for m in R3]
    def dvec(x):
        g=gf(x,av); G=Gf(x,av); c=cf(x); b=bvec(x)
        db=[mp.diff(lambda z: bvec(z)[m],x) for m in R3]
        return [c[1]*db[m]+sum(G[m][n][l]*c[n]*b[l] for n in R3 for l in R3) for m in R3]
    g=gf(x0,av); dot=lambda p,q: sum(g[i,j]*p[i]*q[j] for i in R3 for j in R3)
    c=cf(x0); b=bvec(x0); d=dvec(x0)
    nuf=lambda x: sum(gf(x,av)[i,j]*cf(x)[i]*cf(x)[j] for i in R3 for j in R3)
    nup=lambda x: cf(x)[1]*mp.diff(nuf,x)
    n1=nup(x0); n2=cf(x0)[1]*mp.diff(nup,x0)
    B0=dot(b,b); B1=2*dot(b,d); D0=dot(d,d)
    val=4*(-4*B0**3+4*B0**2*n2-2*B0*B1*n1-B0*n2**2+B1*n1*n2-D0*n1**2)/n1**4
    return val,(n1,n2,B0,B1,D0),dot(c,b),n1/2
print('a w lam | r_null | tau1*^2 (analytic) | tau1^2 numeric at dr=1e-8 | check c.b=nu\'/2')
for (av,wv,lv) in [('0.5','0.1','-0.45'),('0.5','0.1','-0.2'),('0.9','0.1','-0.3'),('1.0','0.09','-0.2'),('1.5','0.1','-0.3'),('0.3','0.2','-0.3'),('0','0.1','-0.3')]:
    av=mp.mpf(av); w=mp.mpf(wv); lam=mp.mpf(lv)
    rl,f=c2.null_radius('BL','BL',av,w,lam)
    cf=cspiral('BL',av,w,lam,ref='BL')
    val,pars,cb,half=tau1_star('BL',cf,rl,av)
    k2,t2,N2=fs('BL',cf,rl+mp.mpf('1e-8'),av)
    print('%s %s %s | %s | %s | %s | c.b=%s nu\'/2=%s'%(mp.nstr(av,3),wv,lv,mp.nstr(rl,8),mp.nstr(val,10),mp.nstr(t2,10),mp.nstr(cb,8),mp.nstr(half,8)),flush=True)
