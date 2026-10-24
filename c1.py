import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fast import *
import numpy as np
M=1
# ---- (5a) lambda->0 vs Eq.(12) kappa^2 and Killing continuity
def eq12(x,av,w):
    Dv=x**2+av**2-2*x
    f=1-(x**2+av**2)*w**2-2*(av*w-1)**2/x
    br=((av*w-1)**2-x**3*w**2)
    return Dv*br**2/(x**6*f**2)
print('--- lambda->0 vs Eq.(12)')
for (av,w,x) in [('0.5','0.15','3.0'),('1.0','-0.1','3.5'),('1.5','0.25','2.0')]:
    av=mp.mpf(av); w=mp.mpf(w); x=mp.mpf(x)
    k0,t0,N0=fs('BL',cspiral('BL',av,w,mp.mpf(0)),x,av)
    k1,t1,N1=fs('BL',cspiral('BL',av,w,mp.mpf('1e-8')),x,av)
    print('a=%s w=%s r=%s | Eq12 k2=%s | code lam=0 k2=%s tau1^2=%s | lam=1e-8 k2=%s tau1^2=%s'%(av,w,x,mp.nstr(eq12(x,av,w),12),mp.nstr(k0,12),mp.nstr(t0,10),mp.nstr(k1,12),mp.nstr(t1,10)))
# ---- (5b) Schwarzschild static
print('--- a=0 static: kappa^2 vs 1/(r^4(1-2/r)),  tau1')
for x in ['3','5']:
    x=mp.mpf(x); k,t,N=fs('BL',cspiral('BL',mp.mpf(0),mp.mpf(0),mp.mpf(0)),x,mp.mpf(0))
    print(x,mp.nstr(k,12),mp.nstr(1/(x**4*(1-2/x)),12),'tau1^2=',mp.nstr(t,5))
# ---- (1) figure: same parameters, BL-form vs KS-form
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
av=mp.mpf('0.5'); w=mp.mpf('0.1'); lam=mp.mpf('-0.45'); rp=1+mp.sqrt(1-av**2)
N2BL=lambda x: kin('BL',cspiral('BL',av,w,lam,ref='BL'),x,av)[3]
rl=mp.findroot(N2BL,2.0)
def sweep(name,ref,rs):
    out=[]
    for x in rs:
        try:
            k2,t2,N2=fs(name,cspiral(name,av,w,lam,ref=ref),mp.mpf(x),av)
            if N2>0 and k2>0 and t2>0: out.append((x,float(mp.sqrt(k2)),float(mp.sqrt(t2)),float(N2)))
        except Exception as e: pass
    return np.array(out)
rsBL=np.concatenate([float(rl)+np.logspace(-4,np.log10(2.0),60)])
rsKS=np.linspace(3.2,1.62,120)
B=sweep('BL','BL',rsBL); K=sweep('KS','KS',rsKS)
print('r_l(BL-form)=',mp.nstr(rl,8),' r+=',mp.nstr(rp,8),' KS-form null near',K[-1,0] if len(K) else None)
fig,ax=plt.subplots(1,2,figsize=(10,4))
for i,(idx,lab) in enumerate([(1,r'$\kappa$'),(2,r'$\tau_1$')]):
    ax[i].semilogy(B[:,0],B[:,idx],'r-',label='BL-form tangent, Eq. (13)')
    ax[i].semilogy(K[:,0],K[:,idx],'b-',label='KS-form tangent, Eq. (24)')
    ax[i].axvline(float(rp),color='k',ls='--',label=r'$r_+$'); ax[i].axvline(float(rl),color='r',ls=':',label=r'BL curve null ($r_\ell$)')
    ax[i].set_xlabel('r/M'); ax[i].set_ylabel(lab); ax[i].set_xlim(1.55,4.0)
ax[0].legend(fontsize=8); fig.suptitle(r'$a=0.5,\ \omega=0.1,\ \lambda=-0.45$: two different worldlines')
plt.tight_layout(); plt.savefig('fig_two_worldlines.png',dpi=150)
np.save('B.npy',B); np.save('K.npy',K)
