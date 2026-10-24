import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import fast
from fast import *
from c3 import beta
import sympy as sp
mp.mp.dps=30
# ---- beta around the a=2.5 bump
av=mp.mpf('2.5'); w=mp.mpf('0.3'); lam=mp.mpf('-0.15'); cf=cspiral('KS',av,w,lam)
print('--- a=2.5 spiral: tau1 spike vs beta (regular frame) around r=0.78')
for x in ['0.6','0.7','0.78','0.86','1.0']:
    x=mp.mpf(x); k2,t2,N2=fs('KS',cf,x,av); b=beta(cf,x,av)[0]
    print('r=%s kappa=%s tau1=%s beta=%s'%(x,mp.nstr(mp.sqrt(k2),5),mp.nstr(mp.sqrt(t2),5),mp.nstr(b,6)))
# ---- Killing near-horizon universality: extremal Kerr and Kerr-Newman (BL, equatorial, M=1)
def add_chart(name,g):
    gi=g.inv(); dg=g.diff(r)
    G=[[[sum(gi[m,s]*((dg[s,l] if n==1 else 0)+(dg[s,n] if l==1 else 0)-(dg[n,l] if s==1 else 0))/2 for s in range(3)) for l in range(3)] for n in range(3)] for m in range(3)]
    fast.FUN[name]=(sp.lambdify((r,a),g,'mpmath'),sp.lambdify((r,a),G,'mpmath'),sp.lambdify((r,a),gi,'mpmath'))
Q2=sp.Rational(1,4)   # Q^2=0.25
DQ=r**2+a**2+Q2-2*r
gKN=sp.Matrix([[1-(2*r-Q2)/r**2,0,a*(2*r-Q2)/r**2],[0,-r**2/DQ,0],[a*(2*r-Q2)/r**2,0,-(r**2+a**2+a**2*(2*r-Q2)/r**2)]])
add_chart('KN',gKN)
print('--- Kerr-Newman (Q^2=0.25, a=0.5), Killing observer omega=omega_H: check kappa*A -> kappa_s, tau1 -> 0')
av=mp.mpf('0.5'); Q2v=mp.mpf('0.25'); rpl=1+mp.sqrt(1-av**2-Q2v); rmi=1-mp.sqrt(1-av**2-Q2v)
wH=av/(rpl**2+av**2); ks=(rpl-rmi)/(2*(rpl**2+av**2)); print('kappa_s=',mp.nstr(ks,8),' omega_H=',mp.nstr(wH,8))
ck=lambda x:[mp.mpf(1),mp.mpf(0),wH]
for d in ['1e-3','1e-5','1e-7']:
    x=rpl+mp.mpf(d); k2,t2,N2=fs('KN',ck,x,av); print('  dr=%s kappa*A/kappa_s=%s tau1=%s'%(d,mp.nstr(mp.sqrt(k2)*mp.sqrt(N2)/ks,8),mp.nstr(mp.sqrt(t2),5)))
print('--- extremal Kerr a=1: Killing observer omega=omega_H=1/2 (kappa_s=0)')
av=mp.mpf(1); ck=lambda x:[mp.mpf(1),mp.mpf(0),mp.mpf('0.5')]
for d in ['1e-2','1e-3','1e-4']:
    x=1+mp.mpf(d); k2,t2,N2=fs('BL',ck,x,av); print('  dr=%s A=%s kappa=%s tau1=%s kappa*A=%s'%(d,mp.nstr(mp.sqrt(N2),5),mp.nstr(mp.sqrt(k2),6),mp.nstr(mp.sqrt(t2),6),mp.nstr(mp.sqrt(k2*N2),6)))
