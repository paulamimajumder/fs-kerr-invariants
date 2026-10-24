import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fast import *
import c3
from c3 import beta
mp.mp.dps=40
av=mp.mpf('0.5'); rp=1+mp.sqrt(1-av**2)
# ISCO (prograde), Bardeen
Z1=1+(1-av**2)**(mp.mpf(1)/3)*((1+av)**(mp.mpf(1)/3)+(1-av)**(mp.mpf(1)/3)); Z2=mp.sqrt(3*av**2+Z1**2)
ri=3+Z2-mp.sqrt((3-Z1)*(3+Z1+2*Z2))
den=mp.sqrt(1-3/ri+2*av/ri**mp.mpf(1.5))
E=(1-2/ri+av/ri**mp.mpf(1.5))/den; L=mp.sqrt(ri)*(1-2*av/ri**mp.mpf(1.5)+av**2/ri**2)/den
print('r_isco=',mp.nstr(ri,10),'E=',mp.nstr(E,10),'L=',mp.nstr(L,10))
def plunge(x):
    Dv=x**2+av**2-2*x; P=E*(x**2+av**2)-av*L
    Rr=P**2-Dv*(x**2+(L-av*E)**2)
    ur=-mp.sqrt(Rr)/x**2
    ut=((x**2+av**2)*P/Dv-av*(av*E-L))/x**2
    up=(av*P/Dv-(av*E-L))/x**2
    return [ut+2*x/Dv*ur, ur, up+av/Dv*ur]   # KS components
cg=plunge
g=FUN['KS'][0]
for x in ['4.0','2.5','1.9']:
    x=mp.mpf(x); c=cg(x); gm=g(x,av); print('norm check r=%s: u.u=%s'%(x,mp.nstr(sum(gm[i,j]*c[i]*c[j] for i in R3 for j in R3),12)))
print('--- geodesic plunge from ISCO (E,L of ISCO): kappa^2 (should be ~0) and regular-frame precession rate beta')
for x in ['4.5','3.0','2.2','1.95','1.8660254','1.8','1.5','1.0']:
    x=mp.mpf(x)
    if x>=ri: continue
    k2,t2,N2=fs('KS',cg,x,av)
    b=beta(cg,x,av)[0]
    print('r=%s kappa^2=%s  beta=%s  |dr/dtau|=%s'%(mp.nstr(x,8),mp.nstr(k2,3),mp.nstr(b,10),mp.nstr(abs(cg(x)[1]),6)),flush=True)
