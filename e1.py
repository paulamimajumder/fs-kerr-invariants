import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from fast import *
mp.mp.dps=40
def Nfun(name,cf,av):
    gf=FUN[name][0]
    return lambda x: sum(gf(x,av)[i,j]*cf(x)[i]*cf(x)[j] for i in R3 for j in R3)
def tau1_star(name,cf,x0,av):
    gf,Gf,_=FUN[name]
    def bvec(x):
        G=Gf(x,av); c=cf(x); dc=[mp.diff(lambda z: cf(z)[m],x) for m in R3]
        return [c[1]*dc[m]+sum(G[m][n][l]*c[n]*c[l] for n in R3 for l in R3) for m in R3]
    def dvec(x):
        G=Gf(x,av); c=cf(x); b=bvec(x); db=[mp.diff(lambda z: bvec(z)[m],x) for m in R3]
        return [c[1]*db[m]+sum(G[m][n][l]*c[n]*b[l] for n in R3 for l in R3) for m in R3]
    g=gf(x0,av); dot=lambda p,q: sum(g[i,j]*p[i]*q[j] for i in R3 for j in R3)
    b=bvec(x0); d=dvec(x0); nuf=Nfun(name,cf,av); nup=lambda x: cf(x)[1]*mp.diff(nuf,x)
    n1=nup(x0); n2=cf(x0)[1]*mp.diff(nup,x0); B0=dot(b,b); B1=2*dot(b,d); D0=dot(d,d)
    return 4*(-4*B0**3+4*B0**2*n2-2*B0*B1*n1-B0*n2**2+B1*n1*n2-D0*n1**2)/n1**4, n1
cases=[('0.5','0.1','-0.45','1.6','+'),('1.5','0.3','-0.45','2.09','-')]
print('--- v4 Table II cases (KS-form spirals): my numerics vs v4')
for (a_,w_,l_,g0,side) in cases:
    av=mp.mpf(a_); w=mp.mpf(w_); lam=mp.mpf(l_); cf=cspiral('KS',av,w,lam)
    f=Nfun('KS',cf,av); rn=mp.findroot(f,mp.mpf(g0))
    sgn=1 if side=='+' else -1
    d=mp.mpf('1e-9')*sgn
    k2,t2,N2=fs('KS',cf,rn+d,av)
    ts,n1=tau1_star('KS',cf,rn,av)
    print('a=%s r_null=%s | kappa*S^1.5=%s (=|nu\'|/2=%s) | tau1 numeric=%s | tau1* analytic=%s'%(a_,mp.nstr(rn,14),mp.nstr(mp.sqrt(k2)*N2**1.5,10),mp.nstr(abs(n1)/2,10),mp.nstr(mp.sqrt(t2),10),mp.nstr(mp.sqrt(ts),10)))
print('v4: a=0.5: r_null=1.60624052119998 kappa S^1.5=0.0235790882193 tau1=1.42715212120464 ; a=1.5: 2.09459810503911, 0.0669359397948, 3.38803372805812')
# ---- v4 exact polynomials, a=0.5
P=[8069871680,-26772037120,174273128,-156465270888,1833728643912,-3433636613679,62941115872140,-244350791458500,214616875380000,267021212950000,-389737982600000,47045881000000]
R=[370177600,49254400,-43184511776,48627189168,-30703052415,4890095876382,-10003404494220,13645866063000,-352080083490000,1019540279300000,-831474044200000,94091762000000]
Pf=lambda x: mp.polyval([mp.mpf(c) for c in P],x); Rf=lambda x: mp.polyval([mp.mpf(c) for c in R],x)
av=mp.mpf('0.5'); cf=cspiral('KS',av,mp.mpf('0.1'),mp.mpf('-0.45'))
print('--- v4 Eq.(22): tau1^2 = P^2/(4 r^4 R^2) vs my FS code (KS, a=0.5, w=0.1, lam=-0.45)')
for x in ['2.5','2.0','1.8','1.7','1.62']:
    x=mp.mpf(x); t2=fs('KS',cf,x,av)[1]
    print('r=%s  v4 formula tau1^2=%s | code tau1^2=%s'%(mp.nstr(x,5),mp.nstr(Pf(x)**2/(4*x**4*Rf(x)**2),12),mp.nstr(t2,12)))
print('R(r_null)=',mp.nstr(Rf(mp.mpf('1.60624052119998')),8),' limit tau1=',mp.nstr(abs(Pf(mp.mpf('1.60624052119998')))/(2*mp.mpf('1.60624052119998')**2*abs(Rf(mp.mpf('1.60624052119998')))),12))
