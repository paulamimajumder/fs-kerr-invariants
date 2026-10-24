import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fast import *
mp.mp.dps=50
def null_radius(name,ref,av,w,lam,rmax=8.0,step=0.01):
    cf0=cspiral(name,av,w,lam,ref=ref)
    def f(x):
        g=FUN[name][0](x,av); c=cf0(x)
        return sum(g[i,j]*c[i]*c[j] for i in R3 for j in R3)
    x=mp.mpf(rmax); prev=f(x)
    rr=[]
    while x>0.3:
        xn=x-step; cur=f(xn)
        if prev>0 and cur<=0:
            return mp.findroot(f,(xn,x),solver='illinois'),f
        prev=cur; x=xn
    return None,f
cases=[('0.5','0.1','-0.45'),('0.5','0.1','-0.2'),('0.9','0.1','-0.3'),('1.0','0.09','-0.2'),('1.0','0.3','-0.45'),
       ('1.5','0.01','-0.2'),('1.5','0.1','-0.3'),('2.5','0.05','-0.15'),('0.3','0.2','-0.3'),('0','0.1','-0.3')]
print('BL-form curves approaching their null radius from outside')
print('a w lam | r_null | k*N^3 vs |c^r dN2/dr|/2 (delta=1e-8) | tau1 at N2~1e-4,1e-6,1e-8 | exponent p (tau1~N^p)')
for (av,wv,lv) in cases:
    av=mp.mpf(av); w=mp.mpf(wv); lam=mp.mpf(lv)
    rl,f=null_radius('BL','BL',av,w,lam)
    if rl is None: print(av,wv,lv,'no null radius found'); continue
    cf=cspiral('BL',av,w,lam,ref='BL')
    vals=[]
    for d in ('1e-4','1e-6','1e-8'):
        x=rl+mp.mpf(d); k2,t2,N2=fs('BL',cf,x,av); vals.append((N2,mp.sqrt(k2),mp.sqrt(t2) if t2>0 else mp.mpf('nan')))
    N2,k,t=vals[-1]
    dN2=mp.diff(f,rl+mp.mpf('1e-8')); pred=abs(cf(rl)[1]*dN2)/2
    p=mp.log(vals[2][2]/vals[1][2])/mp.log(mp.sqrt(vals[2][0])/mp.sqrt(vals[1][0]))
    print('%s %s %s | %s | %s vs %s | %s, %s, %s | p=%s'%(mp.nstr(av,3),wv,lv,mp.nstr(rl,8),mp.nstr(k*N2**1.5,6),mp.nstr(pred,6),mp.nstr(vals[0][2],6),mp.nstr(vals[1][2],6),mp.nstr(vals[2][2],6),mp.nstr(p,4)),flush=True)
