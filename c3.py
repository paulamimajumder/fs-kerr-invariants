import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fast import *
mp.mp.dps=40
name='KS'
def eps_cross(g,u,w1):
    gl=lambda v:[sum(g[i,j]*v[j] for j in R3) for i in R3]
    ul=gl(u); wl=gl(w1)
    det=abs(mp.det(g)); s=mp.sqrt(det)
    e=lambda i,j,k: (1 if (i,j,k) in [(0,1,2),(1,2,0),(2,0,1)] else (-1 if (i,j,k) in [(0,2,1),(2,1,0),(1,0,2)] else 0))
    return [sum(e(m,n,l)*ul[n]*wl[l] for n in R3 for l in R3)/s for m in R3]
def frame(cf,x,av):
    g,G,u,N2=kin(name,cf,x,av)
    l=[mp.mpf(1),mp.mpf(-1),mp.mpf(0)]  # ingoing principal null dir in KS (t',r,phi')
    lu=sum(g[i,j]*l[i]*u[j] for i in R3 for j in R3)
    w1=[(l[i]-lu*u[i])/abs(lu) for i in R3]
    return g,G,u,N2,w1
def beta(cf,x,av):
    """gyro precession rate relative to the regular frame w1 = projection of ingoing principal null dir; returns beta, w2 norm"""
    g,G,u,N2,w1=frame(cf,x,av)
    w2=eps_cross(g,u,w1)
    dot=lambda p,q: sum(g[i,j]*p[i]*q[j] for i in R3 for j in R3)
    dw1=[mp.diff(lambda z: frame(cf,z,av)[4][m],x) for m in R3]
    Dw1=[u[1]*dw1[m]+sum(G[m][n][l]*u[n]*w1[l] for n in R3 for l in R3) for m in R3]
    return dot(w2,Dw1),dot(w2,w2),dot(w1,w1),dot(u,w1),w1,w2,u,g,u[1]
def theta_e1(cf,x,av):
    """angle of FS e1 relative to w1 in the (w1,w2) plane"""
    g,G,u,N2,w1=frame(cf,x,av); w2=eps_cross(g,u,w1)
    dot=lambda p,q: sum(g[i,j]*p[i]*q[j] for i in R3 for j in R3)
    def e1f(y):
        gg,GG,uu,_=kin(name,cf,y,av)
        du=[mp.diff(lambda z: kin(name,cf,z,av)[2][m],y) for m in R3]
        ac=[uu[1]*du[m]+sum(GG[m][n][l]*uu[n]*uu[l] for n in R3 for l in R3) for m in R3]
        k=mp.sqrt(-sum(gg[i,j]*ac[i]*ac[j] for i in R3 for j in R3)); return [q/k for q in ac]
    def th(y):
        gg,GG,uu,_,ww1=frame(cf,y,av); ww2=eps_cross(gg,uu,ww1); e=e1f(y)
        d=lambda p,q: sum(gg[i,j]*p[i]*q[j] for i in R3 for j in R3)
        return mp.atan2(-d(e,ww2),-d(e,ww1))
    return th(x), mp.diff(th,x), u[1]
av=mp.mpf('0.5'); w=mp.mpf('0.1'); lam=mp.mpf('-0.45')
cf=cspiral('KS',av,w,lam)
print('=== consistency: beta (gyro rate rel. w1) vs dtheta/dtau +- tau1 ; spiral a=.5 w=.1 lam=-.45 (KS chart)')
for x in ['3.0','2.0','1.8']:
    x=mp.mpf(x); b,n2,n1,uw,*_=beta(cf,x,av); th,dth,ur=theta_e1(cf,x,av)
    k2,t2,N2=fs('KS',cf,x,av); t1=mp.sqrt(t2)
    print('r=',x,' beta=',mp.nstr(b,8),' dtheta/dr*u^r=',mp.nstr(dth*ur,8),' tau1=',mp.nstr(t1,8),' |beta-dth|=',mp.nstr(abs(b-dth*ur),8),' |beta+dth|?',mp.nstr(abs(b+dth*ur),6),'  norms w2,w1,u.w1:',mp.nstr(n2,4),mp.nstr(n1,4),mp.nstr(uw,3))
