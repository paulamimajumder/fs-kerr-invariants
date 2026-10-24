import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from fast import *
    import e1
from e1 import tau1_star, Nfun
mp.mp.dps=50
cases=[('BL','0.5','0.1','-0.45','2.0869','+'),('BL','0.5','0.1','-0.2','1.9551','+'),('BL','0.9','0.1','-0.3','1.7645','+'),
       ('BL','1.0','0.09','-0.2','1.7208','+'),('BL','1.5','0.1','-0.3','1.5190','+'),('BL','0.3','0.2','-0.3','1.7825','+'),
       ('BL','0','0.1','-0.3','2.2101','+'),('KS','0.5','0.1','-0.45','1.6062','+'),('KS','1.5','0.3','-0.45','2.0946','-')]
print('chart a w lam | kappa*S^1.5/|Q0| at the Table III distance | tau1*^2 (s=t)  | tau1*^2 (s=r)  | tau1*^2 (arbitrary reparam h=1+0.3r+0.1r^2)')
for (ch,a_,w_,l_,g0,side) in cases:
    av=mp.mpf(a_); w=mp.mpf(w_); lam=mp.mpf(l_)
    cf=cspiral(ch,av,w,lam,ref=ch)
    f=Nfun(ch,cf,av); rn=mp.findroot(f,mp.mpf(g0))
    d=(mp.mpf('1e-8') if ch=='BL' else mp.mpf('1e-9'))*(1 if side=='+' else -1)
    k2,t2,N2=fs(ch,cf,rn+d,av)
    ts,n1=tau1_star(ch,cf,rn,av)
    ratio=mp.sqrt(k2)*N2**mp.mpf('1.5')/(abs(n1)/2)
    # reparametrisations: v -> h(r) v  (changes s, S, Q, q, p)
    cr=lambda x: [q/cf(x)[1] for q in cf(x)]                    # parameter s = r
    ch2=lambda x: [q*(1+mp.mpf('0.3')*x+mp.mpf('0.1')*x**2) for q in cf(x)]
    tr,_=tau1_star(ch,cr,rn,av); th,_=tau1_star(ch,ch2,rn,av)
    print('%s %s %s %s | %s | %s | %s | %s'%(ch,a_,w_,l_,mp.nstr(ratio,10),mp.nstr(ts,10),mp.nstr(tr,10),mp.nstr(th,10)),flush=True)
