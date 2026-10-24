import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fast import *
import c3, c4
from c3 import eps_cross
mp.mp.dps=35
av=mp.mpf('0.5')
# ---------- (A) cross-check of manuscript appendix A.1 (BL) kappa against FS code
def appendixA(x,av,w,lam,m=1):
    D=lambda y: y**2-2*m*y+av**2
    FK=lambda y: 1-2*m/y+4*av*m*w/y+(-av**2-2*av**2*m/y-y**2)*w**2-lam**2*y**4*w**2/D(y)
    Fp=mp.diff(FK,x); F=FK(x); d=D(x)
    PK=-lam*x*w**2*Fp/(2*F**2)+2*av*lam*m*w/(x*d*F)+2*lam*(-av**2*m*x**2+x**4*(-2*m+x))*w**2/(x**3*d*F)
    QK=-lam*x*w*Fp/(2*F**2)+2*lam*m*(av**2+x**2)*w/(x*d*F)-2*av*lam*m*(av**2+3*x**2)*w**2/(x*d*F)
    BK=lambda y: lam*y*w/mp.sqrt(FK(y)); Bp=mp.diff(BK,x)
    Xi=m*d/(x**4*F)-2*av*m*d*w/(x**4*F)+lam**2*x*(av**2-m*x)*w**2/(d*F)-(-av**2*m*x**2+x**5)*d*w**2/(x**6*F)+Bp**2
    k2=-(av**2+2*av**2*m/x+x**2)*PK**2-4*av*m/x*PK*QK-(1-2*m/x)*QK**2+x**2/d*Xi
    return k2
print('--- (A) Appendix A.1 kappa^2 (manuscript) vs FS code (BL, Eq.13-form curve)')
for (a_,w_,l_,x_) in [('0.5','0.1','-0.45','3.0'),('1.0','0.09','-0.2','2.5'),('0.5','0.1','-0.2','2.2')]:
    a1=mp.mpf(a_); w1=mp.mpf(w_); l1=mp.mpf(l_); x1=mp.mpf(x_)
    k2c=fs('BL',cspiral('BL',a1,w1,l1,ref='BL'),x1,a1)[0]
    print(a_,w_,l_,x_,'| appendix k2=',mp.nstr(appendixA(x1,a1,w1,l1),12),'| code k2=',mp.nstr(k2c,12))
# ---------- (E) Marck frame comparison, equatorial geodesic plunge (a=0.5)
rp=1+mp.sqrt(1-av**2)
ri,E,L=c4.ri,c4.E,c4.L
cg=c4.plunge
Kc=(L-av*E)**2
def marck(x):
    P=E*(x**2+av**2)-av*L
    return mp.sqrt(Kc)/x**2*(P/(Kc+x**2)+av*(L-av*E)/Kc)
def frameM(x):
    g,G,u,N2=kin('KS',cg,x,av)
    lin=[mp.mpf(1),mp.mpf(-1),mp.mpf(0)]
    Dv=x**2+av**2-2*x
    lout=[x**2+av**2+2*x,Dv,2*av]
    dot=lambda p,q: sum(g[i,j]*p[i]*q[j] for i in R3 for j in R3)
    A=dot(lout,u); B=dot(lin,u)
    w=[A*lin[i]-B*lout[i] for i in R3]
    nrm=mp.sqrt(abs(dot(w,w))); w=[q/nrm for q in w]
    return g,G,u,w,dot(lout,lout)
def rate(fr,x):
    g,G,u,w,_=fr(x)
    w2=eps_cross(g,u,w); dot=lambda p,q: sum(g[i,j]*p[i]*q[j] for i in R3 for j in R3)
    dw=[mp.diff(lambda z: fr(z)[3][m],x) for m in R3]
    Dw=[u[1]*dw[m]+sum(G[m][n][l]*u[n]*w[l] for n in R3 for l in R3) for m in R3]
    return dot(w2,Dw)
print('--- (E) regular Marck-type axis (unit vector in Carter t-r plane orthogonal to u) vs Marck formula; plunge a=0.5')
for x in ['3.0','2.2','1.8660254','1.5','1.0']:
    x=mp.mpf(x); g,G,u,w,lo=frameM(x)
    bm=rate(frameM,x); b1=c3.beta(cg,x,av)[0]
    print('r=%s | rate rel. Carter-plane axis=%s | Marck Psi_dot=%s | beta (principal-null frame)=%s | l_out null check=%s'%(mp.nstr(x,6),mp.nstr(bm,10),mp.nstr(marck(x),10),mp.nstr(b1,10),mp.nstr(lo,3)),flush=True)
# ---------- (B) extremal Kerr, Killing observers
print('--- (B) extremal Kerr a=1: A^2 at omega=1/2 and light surfaces for omega<1/2')
for x in ['1.001','1.1','2.0','5.0']:
    x=mp.mpf(x); A2=1-(x**2+1)*mp.mpf('0.25')-2*(1-mp.mpf('0.5'))**2/x
    print('r=%s A^2=%s  -(r-1)^2(r+2)/(4r)=%s'%(x,mp.nstr(A2,8),mp.nstr(-(x-1)**2*(x+2)/(4*x),8)))
import sympy as sp
rr,ww=sp.symbols('r w',positive=True)
aa=1
A2s=1-(rr**2+aa**2)*ww**2-2*(1-aa*ww)**2/rr
Nk=sp.sqrt(rr**2+aa**2-2*rr)*(rr**3*ww**2-(1-aa*ww)**2); Nt=aa*(1+ww**2*(aa**2+3*rr**2))-2*aa**2*ww+ww*rr**2*(rr-3)
for wv in [0.45,0.49,0.499,0.4999]:
    f=sp.lambdify(rr,A2s.subs(ww,wv),'mpmath')
    r0=None
    for g0 in (1.0005,1.01,1.1,1.5,3.0,6.0):
        try:
            z=mp.findroot(f,g0)
            if z>1 and abs(f(z))<1e-20: r0=z;break
        except Exception: pass
    if r0 is None: print('w=%s no light surface found'%wv); continue
    print('w=%s: outermost-found light surface r_l=%s ; kappa*A^2=%s tau1*A^2=%s (closed forms, /r^3)'%(wv,mp.nstr(r0,10),mp.nstr(abs(Nk.subs({ww:wv,rr:r0}).evalf(20))/r0**3,8),mp.nstr(abs(Nt.subs({ww:wv,rr:r0}).evalf(20))/r0**3,8)))
