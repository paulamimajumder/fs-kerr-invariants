import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from fast import *
import sympy as sp
import fast
mp.mp.dps=40
# --- regularity of the Doran-type chart at the horizon (equatorial block)
rr,aa=fast.r,fast.a
Dl=rr**2+aa**2-2*rr
F=sp.sqrt(2*rr*(rr**2+aa**2))/Dl; G=aa*sp.sqrt(2*rr/(rr**2+aa**2))/Dl
J=sp.eye(3); J[0,1]=-F; J[2,1]=-G
gD=sp.simplify(J.T*fast.gBL*J)
fD=sp.lambdify((rr,aa),gD,'mpmath')
fK=fast.FUN['KS'][0]; fI=fast.FUN['IK'][0]
av=mp.mpf('0.5'); rp=1+mp.sqrt(1-av**2)
print('--- metric components of the three charts near r+ (a=0.5): max |g_ij| at r+ +- 1e-9')
for nm,f in [('KS',fK),('IK',fI),('DO',fD)]:
    vals=[]
    for s in (1,-1):
        g=f(rp+s*mp.mpf('1e-9'),av); vals.append(max(abs(g[i,j]) for i in R3 for j in R3))
    print(nm,[mp.nstr(v,6) for v in vals],' det at r+ + 1e-9 =',mp.nstr(mp.det(f(rp+mp.mpf('1e-9'),av)),8))
# --- regenerate v4-style figures with my code
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def nrm(cf,av): 
    gf=FUN['KS'][0]; return lambda x: sum(gf(x,av)[i,j]*cf(x)[i]*cf(x)[j] for i in R3 for j in R3)
cases=[('0.5','0.1','-0.45',1,'a/M=0.5','C0'),('1.5','0.3','-0.45',-1,'a/M=1.5','C1')]
fig1,ax1=plt.subplots(figsize=(5.5,4)); fig2,ax2=plt.subplots(figsize=(5.5,4)); fig3,ax3=plt.subplots(figsize=(5.5,4))
for (a_,w_,l_,side,lab,col) in cases:
    av=mp.mpf(a_); cf=cspiral('KS',av,mp.mpf(w_),mp.mpf(l_)); f=nrm(cf,av)
    rn=mp.findroot(f,mp.mpf('1.6') if side==1 else mp.mpf('2.09'))
    # Q0 and tau1*: reuse analytic via numerical limit
    ds=np.logspace(-9,-1.7,22); S=[];K=[];T=[]
    for d in ds:
        k2,t2,N2=fs('KS',cf,rn+side*mp.mpf(float(d)),av); S.append(float(N2));K.append(float(mp.sqrt(k2)));T.append(float(mp.sqrt(t2)))
    S=np.array(S);K=np.array(K);T=np.array(T)
    Q0=K[0]*S[0]**1.5
    ax1.loglog(S,K,'o',ms=3,color=col,label=lab); ax1.loglog(S,Q0*S**-1.5,'--',color=col)
    ax2.semilogx(ds,T,'-',color=col,label=lab); ax2.axhline(T[0],ls='--',color=col)
    rs=rn+side*np.linspace(0,0.2,60); Sv=[float(f(mp.mpf(float(x)))) for x in rs]
    ax3.plot(rs,Sv,color=col,label=lab); ax3.axvline(float(rn),ls=':',color=col)
    print(lab,'r_null=',mp.nstr(rn,12),' Q0~',Q0,' tau1(smallest delta)=',T[0])
ax1.set_xlabel(r'$S=g(v,v)$'); ax1.set_ylabel(r'$\kappa M$'); ax1.legend(title='dashed: $|Q_0|S^{-3/2}$'); fig1.tight_layout(); fig1.savefig('fig_curvature_growth.png',dpi=150)
ax2.set_xlabel(r'$|r-r_{\rm null}|/M$'); ax2.set_ylabel(r'$\tau_1 M$'); ax2.invert_xaxis(); ax2.legend(title='dashed: limit'); fig2.tight_layout(); fig2.savefig('fig_torsion_limits.png',dpi=150)
ax3.set_xlabel('r/M'); ax3.set_ylabel(r'$S=g(v,v)$'); ax3.legend(title='dotted: null radius'); fig3.tight_layout(); fig3.savefig('fig_timelike_domain.png',dpi=150)
