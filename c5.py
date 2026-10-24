import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fast import *
import c3, c4
from c3 import beta
mp.mp.dps=30
# ---- (1) the bump in Fig.6 (a=2.5, w=.3, lam=-.15), KS chart
av=mp.mpf('2.5'); w=mp.mpf('0.3'); lam=mp.mpf('-0.15'); cf=cspiral('KS',av,w,lam)
best=None
rows=[]
for x in [mp.mpf(v)/1000 for v in range(600,1001,20)]:
    k2,t2,N2=fs('KS',cf,x,av); rows.append((x,mp.sqrt(k2),mp.sqrt(t2)))
kmin=min(rows,key=lambda t:t[1]); tmax=max(rows,key=lambda t:t[2])
print('--- Fig.6 a=2.5 bump: kappa_min=%s at r=%s ; tau1_max=%s at r=%s ; kappa*tau1 at tau1_max = %s'%(mp.nstr(kmin[1],6),mp.nstr(kmin[0],4),mp.nstr(tmax[2],6),mp.nstr(tmax[0],4),mp.nstr(tmax[1]*tmax[2],6)))
# refine: minimise kappa
xm=mp.findroot(lambda y: mp.diff(lambda z: fs('KS',cf,z,av)[0],y),kmin[0])
k2,t2,N2=fs('KS',cf,xm,av); print('refined: r_min=%s kappa^2_min=%s tau1=%s'%(mp.nstr(xm,8),mp.nstr(k2,6),mp.nstr(mp.sqrt(t2),6)))
# ---- (2) accumulated angles
print('--- accumulated angles (proper-time integrals)')
av=mp.mpf('0.5'); w=mp.mpf('0.1'); lam=mp.mpf('-0.45'); cf=cspiral('KS',av,w,lam)
mp.mp.dps=20
def tau_int(fun,r0,r1,n=14):
    # integral over r from r0 down to r1 of fun(r)*dtau/dr ; dtau/dr = N/|c^r|
    f=lambda x: fun(x)*mp.sqrt(kin('KS',cf,x,av)[3])/abs(cf(x)[1])
    return mp.quad(f,mp.linspace(r1,r0,4))
T1=tau_int(lambda x: mp.sqrt(fs('KS',cf,x,av)[1]),mp.mpf(3),mp.mpf('1.7'))
B1=tau_int(lambda x: beta(cf,x,av)[0],mp.mpf(3),mp.mpf('1.7'))
dtau=tau_int(lambda x: 1,mp.mpf(3),mp.mpf('1.7'))
print('KS spiral a=.5 w=.1 lam=-.45, r: 3.0 -> 1.7 (crosses r+=1.866): proper time=%s  int tau1 dtau=%s  int beta dtau=%s'%(mp.nstr(dtau,6),mp.nstr(T1,6),mp.nstr(B1,6)))
# plunge
av=mp.mpf('0.5')
cg=c4.plunge
def tau_int_g(fun,r0,r1):
    f=lambda x: fun(x)/abs(cg(x)[1])
    return mp.quad(f,mp.linspace(r1,r0,4))
Bg=tau_int_g(lambda x: beta(cg,x,av)[0],mp.mpf(3),mp.mpf(1)); dtg=tau_int_g(lambda x:1,mp.mpf(3),mp.mpf(1))
print('geodesic plunge r: 3.0 -> 1.0: proper time=%s  int beta dtau=%s (rad)'%(mp.nstr(dtg,6),mp.nstr(Bg,6)))
