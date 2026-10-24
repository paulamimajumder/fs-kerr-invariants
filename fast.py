import sympy as sp, mpmath as mp
mp.mp.dps=40
r,a=sp.symbols('r a')
D=r**2+a**2-2*r   # M=1
gBL=sp.Matrix([[1-2/r,0,2*a/r],[0,-r**2/D,0],[2*a/r,0,-(r**2+a**2+2*a**2/r)]])
def pull(F,G):
    J=sp.eye(3); J[0,1]=-F; J[2,1]=-G
    return J.T*gBL*J
DOF=sp.sqrt(2*r*(r**2+a**2))/D; DOG=a*sp.sqrt(2*r/(r**2+a**2))/D
CH={'BL':(gBL,None,None),'KS':(pull(2*r/D,a/D),2*r/D,a/D),
    'IK':(pull((r**2+a**2)/D,a/D),(r**2+a**2)/D,a/D),
    'DO':(pull(DOF,DOG),DOF,DOG)}
FUN={}
def build(name):
    g=CH[name][0]; gi=g.inv(); dg=g.diff(r)
    G=[[[sum(gi[m,s]*((dg[s,l] if n==1 else 0)+(dg[s,n] if l==1 else 0)-(dg[n,l] if s==1 else 0))/2 for s in range(3)) for l in range(3)] for n in range(3)] for m in range(3)]
    FUN[name]=(sp.lambdify((r,a),g,'mpmath'),sp.lambdify((r,a),G,'mpmath'),sp.lambdify((r,a),gi,'mpmath'))
for nm in CH: build(nm)
R3=range(3)
def kin(name,cf,x,av):
    gf,Gf,_=FUN[name]; g=gf(x,av); G=Gf(x,av); c=cf(x)
    N2=sum(g[i,j]*c[i]*c[j] for i in R3 for j in R3)
    u=[ci/mp.sqrt(N2) for ci in c]
    return g,G,u,N2
def fs(name,cf,x,av):
    """kappa^2, tau1^2, N^2 for equatorial curve with tangent comps cf(r) (functions of r only)"""
    def accf(y):
        g,G,u,N2=kin(name,cf,y,av)
        du=[mp.diff(lambda z: kin(name,cf,z,av)[2][m],y) for m in R3]
        return [u[1]*du[m]+sum(G[m][n][l]*u[n]*u[l] for n in R3 for l in R3) for m in R3],g
    def e1f(y):
        ac,g=accf(y); k=mp.sqrt(-sum(g[i,j]*ac[i]*ac[j] for i in R3 for j in R3))
        return [q/k for q in ac]
    g,G,u,N2=kin(name,cf,x,av)
    ac,_=accf(x); k2=-sum(g[i,j]*ac[i]*ac[j] for i in R3 for j in R3); k=mp.sqrt(k2)
    e=e1f(x); de=[mp.diff(lambda z: e1f(z)[m],x) for m in R3]
    De=[u[1]*de[m]+sum(G[m][n][l]*u[n]*e[l] for n in R3 for l in R3) for m in R3]
    V=[De[m]-k*u[m] for m in R3]
    t2=-sum(g[i,j]*V[i]*V[j] for i in R3 for j in R3)
    return k2,t2,N2
def curve(name,w,lam,chart_ref='KS'):
    """worldline fixed by (1, r lam w, 0, w) in chart_ref, expressed in chart `name`"""
    def c(x,av=None):
        return None
    return c
def cspiral(name,av,w,lam,ref='KS'):
    def c(x):
        Dv=x**2+av**2-2*x
        def FGn(nm):
            if nm=='BL': return mp.mpf(0),mp.mpf(0)
            if nm=='KS': return 2*x/Dv, av/Dv
            if nm=='IK': return (x**2+av**2)/Dv, av/Dv
            if nm=='DO': return mp.sqrt(2*x*(x**2+av**2))/Dv, av*mp.sqrt(2*x/(x**2+av**2))/Dv
        Fr,Gr=FGn(ref); Fx,Gx=FGn(name); ur=x*lam*w
        return [1+(Fx-Fr)*ur, ur, w+(Gx-Gr)*ur]
    return c
if __name__=='__main__':
    import time
    av=mp.mpf('0.5'); w=mp.mpf('0.1'); lam=mp.mpf('-0.45')
    t0=time.time()
    for nm,x in [('KS','2.5'),('BL','2.5'),('IK','1.8'),('KS','1.8')]:
        print(nm,x,[mp.nstr(v,15) for v in fs(nm,cspiral(nm,av,w,lam),mp.mpf(x),av)],round(time.time()-t0,1),'s')
    print('earlier (KS 2.5): N2 0.34996875 k2 0.0692327352084786 t2 0.0136248281803597; (1.8): 0.111371222222222 1.70263947002147 0.0639969290288769')
