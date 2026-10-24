import sympy as sp, mpmath as mp, time
mp.mp.dps=30
r,th,a=sp.symbols('r theta a')
Dl=r**2+a**2-2*r; S=r**2+a**2*sp.cos(th)**2; s2=sp.sin(th)**2
gB=sp.zeros(4)
gB[0,0]=1-2*r/S; gB[1,1]=-S/Dl; gB[2,2]=-S; gB[0,3]=gB[3,0]=2*r*a*s2/S
gB[3,3]=-(r**2+a**2+2*r*a**2*s2/S)*s2
e1,e2=sp.Rational(3,10),sp.Rational(1,5)
TR={'BL':(0,0,0,0),'KS':(2*r/Dl,0,a/Dl,0),
    'KT':(2*r/Dl+e1*sp.cos(th)**2,-2*e1*r*sp.cos(th)*sp.sin(th),a/Dl+e2*sp.cos(th)**2,-2*e2*r*sp.cos(th)*sp.sin(th))}
FN={}
for nm,(Tr,Tt,Pr,Pt) in TR.items():
    J=sp.eye(4); J[0,1]=-Tr; J[0,2]=-Tt; J[3,1]=-Pr; J[3,2]=-Pt
    g=J.T*gB*J
    FN[nm]=(sp.lambdify((r,th,a),g,'mpmath'),sp.lambdify((r,th,a),g.diff(r),'mpmath'),sp.lambdify((r,th,a),g.diff(th),'mpmath'),
            sp.lambdify((r,th,a),list(TR[nm]),'mpmath'))
R4=range(4)
def geo(nm,x,t_,av):
    g=FN[nm][0](x,t_,av); gi=g**-1; d1=FN[nm][1](x,t_,av); d2=FN[nm][2](x,t_,av)
    dg=[None,d1,d2]
    def D(k,i,j): return dg[k][i,j] if k in (1,2) else 0
    G=[[[sum(gi[m,s]*(D(n,s,l)+D(l,s,n)-D(s,n,l)) for s in R4)/2 for l in R4] for n in R4] for m in R4]
    return g,G
def master(nm,av,w,lam,mu,form):
    """tangent comps in chart nm. form='KS': (1,r lam w,mu,w) defined in KS chart; form='BL': defined in BL chart"""
    def c(x,t_):
        ur=x*lam*w
        if form=='BL': cb=[1,ur,mu,w]
        else:
            Tr,Tt,Pr,Pt=FN['KS'][3](x,t_,av); cb=[1-Tr*ur,ur,mu,w-Pr*ur]
        Tr,Tt,Pr,Pt=FN[nm][3](x,t_,av)
        return [cb[0]+Tr*ur+Tt*mu, ur, mu, cb[3]+Pr*ur+Pt*mu]
    return c
def kin(nm,cf,x,t_,av):
    g,G=geo(nm,x,t_,av); c=cf(x,t_)
    N2=sum(g[i,j]*c[i]*c[j] for i in R4 for j in R4); u=[q/mp.sqrt(N2) for q in c]
    return g,G,u,N2
def ddir(F,x,t_,u):
    return [u[1]*mp.diff(lambda y,z:F(y,z)[m],(x,t_),(1,0))+u[2]*mp.diff(lambda y,z:F(y,z)[m],(x,t_),(0,1)) for m in R4]
def fs4(nm,cf,x,t_,av):
    def uf(y,z): return kin(nm,cf,y,z,av)[2]
    def accf(y,z):
        g,G,u,N2=kin(nm,cf,y,z,av); du=ddir(uf,y,z,u)
        return [du[m]+sum(G[m][n][l]*u[n]*u[l] for n in R4 for l in R4) for m in R4]
    def e1f(y,z):
        g=FN[nm][0](y,z,av); ac=accf(y,z); k=mp.sqrt(-sum(g[i,j]*ac[i]*ac[j] for i in R4 for j in R4)); return [q/k for q in ac]
    g,G,u,N2=kin(nm,cf,x,t_,av); ac=accf(x,t_); k2=-sum(g[i,j]*ac[i]*ac[j] for i in R4 for j in R4); k=mp.sqrt(k2)
    e=e1f(x,t_); de=ddir(e1f,x,t_,u)
    De=[de[m]+sum(G[m][n][l]*u[n]*e[l] for n in R4 for l in R4) for m in R4]
    V=[De[m]-k*u[m] for m in R4]; t2=-sum(g[i,j]*V[i]*V[j] for i in R4 for j in R4)
    return k2,t2,N2
if __name__=='__main__':
    av=mp.mpf('0.5'); w=mp.mpf('0.1'); lam=mp.mpf('-0.45'); mu=mp.mpf('0.03'); t0=time.time()
    print('--- off-equatorial worldline (theta=1.0 rad, polar drift mu=0.03): chart independence, KS-defined tangent')
    for (x,t_) in [('2.5','1.0'),('1.8','1.0')]:
        for nm in (['BL','KS','KT'] if x=='2.5' else ['KS','KT']):
            k2,t2,N2=fs4(nm,master(nm,av,w,lam,mu,'KS'),mp.mpf(x),mp.mpf(t_),av)
            print('r=%s chart=%s N2=%s kappa2=%s tau1^2=%s (%.0fs)'%(x,nm,mp.nstr(N2,14),mp.nstr(k2,14),mp.nstr(t2,14),time.time()-t0),flush=True)
    # null-limit kappa law off equator (BL-form)
    cB=master('BL',av,w,lam,mu,'BL'); th0=mp.mpf('1.0')
    nuf=lambda y,z: sum(FN['BL'][0](y,z,av)[i,j]*cB(y,z)[i]*cB(y,z)[j] for i in R4 for j in R4)
    xs=None
    xx=mp.mpf(6)
    while xx>1:
        if nuf(xx,th0)>0 and nuf(xx-mp.mpf('0.01'),th0)<=0: xs=mp.findroot(lambda y:nuf(y,th0),(xx-mp.mpf('0.01'),xx),solver='illinois');break
        xx-=mp.mpf('0.01')
    print('--- off-equator BL-form null point at theta=1.0: r_null=',mp.nstr(xs,10))
    c0=cB(xs,th0)
    for d in ['1e-4','1e-6']:
        x=xs+mp.mpf(d); k2,t2,N2=fs4('BL',cB,x,th0,av)
        nup=cB(x,th0)[1]*mp.diff(lambda y:nuf(y,th0),x)+cB(x,th0)[2]*mp.diff(lambda z:nuf(x,z),th0)
        print('  dr=%s kappa*N^3=%s  |nu\'|/2=%s  tau1=%s'%(d,mp.nstr(mp.sqrt(k2)*N2**1.5,8),mp.nstr(abs(nup)/2,8),mp.nstr(mp.sqrt(t2),8)),flush=True)
