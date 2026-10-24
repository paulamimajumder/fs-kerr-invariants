import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from fast import *
    import c3, c4, c2, d2
    from c3 import beta
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
mp.mp.dps=40
# ---------- Figure B: null-limit scaling
av=mp.mpf('0.5'); w=mp.mpf('0.1'); lam=mp.mpf('-0.45')
rl,f=c2.null_radius('BL','BL',av,w,lam); cf=cspiral('BL',av,w,lam,ref='BL')
val,pars,cb,half=d2.tau1_star('BL',cf,rl,av); t1s=float(mp.sqrt(val)); nup=abs(float(pars[0]))
Ns=[];ks=[];ts=[]
for e in np.logspace(-9,-1,17):
    k2,t2,N2=fs('BL',cf,rl+mp.mpf(float(e)),av); Ns.append(float(mp.sqrt(N2))); ks.append(float(mp.sqrt(k2))); ts.append(float(mp.sqrt(t2)))
Ns=np.array(Ns);ks=np.array(ks);ts=np.array(ts)
fig,ax=plt.subplots(1,2,figsize=(10,4))
ax[0].loglog(Ns,ks,'ko',label=r'$\kappa$ (numerical)'); ax[0].loglog(Ns,nup/(2*Ns**3),'r--',label=r'$|S^\prime|/(2N^3)$')
ax[0].set_xlabel(r'$N=\sqrt{S}$'); ax[0].set_ylabel(r'$\kappa$'); ax[0].legend()
ax[1].semilogx(Ns,ts,'ko',label=r'$\tau_1$ (numerical)'); ax[1].axhline(t1s,color='r',ls='--',label=r'$\tau_1^*$ (Proposition 1)')
ax[1].set_xlabel(r'$N$'); ax[1].set_ylabel(r'$\tau_1$'); ax[1].legend()
fig.suptitle(r'BL-form spiral, $a=0.5,\ \omega=0.1,\ \lambda=-0.45$, approaching its null point $r_\ell=%.4f$'%float(rl)); plt.tight_layout()
plt.savefig('fig_null_limit.png',dpi=150)
print('tau1* =',t1s,' |nu\'|=',nup)
# ---------- Figure C: regular-frame rates
av=mp.mpf('0.5'); rp=float(1+mp.sqrt(1-av**2)); cg=c4.plunge
Z1=None
rs=np.linspace(4.2,0.8,22)
bp=[];mk=[]
for x in rs:
    x=mp.mpf(float(x)); bp.append(float(beta(cg,x,av)[0]))
    P=c4.E*(x**2+av**2)-av*c4.L; K=(c4.L-av*c4.E)**2
    mk.append(float(mp.sqrt(K)/x**2*(P/(K+x**2)+av*(c4.L-av*c4.E)/K)))
cs=cspiral('KS',av,mp.mpf('0.1'),mp.mpf('-0.45'))
r2=np.linspace(3.0,1.7,14); t1=[];b2=[]
for x in r2:
    x=mp.mpf(float(x)); t1.append(float(mp.sqrt(fs('KS',cs,x,av)[1]))); b2.append(float(beta(cs,x,av)[0]))
fig,ax=plt.subplots(1,2,figsize=(10,4))
ax[0].plot(rs,bp,'b-o',ms=3,label=r'$\beta$ (principal-null frame)'); ax[0].plot(rs,mk,'g-s',ms=3,label=r'Marck $\dot\Psi$ (Carter-plane axis)')
ax[0].axvline(rp,color='k',ls='--',label=r'$r_+$'); ax[0].set_xlabel('r/M'); ax[0].set_ylabel('rate per proper time'); ax[0].set_title('geodesic plunge from the ISCO, $a=0.5$'); ax[0].legend(fontsize=8)
ax[1].plot(r2,t1,'r-o',ms=3,label=r'$\tau_1$ (FS)'); ax[1].plot(r2,b2,'b-o',ms=3,label=r'$\beta$')
ax[1].axvline(rp,color='k',ls='--'); ax[1].set_xlabel('r/M'); ax[1].set_title(r'KS-defined spiral, $a=0.5,\omega=0.1,\lambda=-0.45$'); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('fig_frame_rates.png',dpi=150)
# Marck accumulated angle
def marck(x):
    P=c4.E*(x**2+av**2)-av*c4.L; K=(c4.L-av*c4.E)**2
    return mp.sqrt(K)/x**2*(P/(K+x**2)+av*(c4.L-av*c4.E)/K)
mp.mp.dps=20
Mk=mp.quad(lambda x: marck(x)/abs(cg(x)[1]),mp.linspace(mp.mpf(1),mp.mpf(3),4))
print('Marck accumulated angle plunge r:3->1 =',mp.nstr(Mk,5))
