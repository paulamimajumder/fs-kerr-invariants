import sys; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fast import *
mp.mp.dps=30
def parts(x,av,w,lam,m=1):
    D=lambda y: y**2-2*m*y+av**2
    FK=lambda y: 1-2*m/y+4*av*m*w/y+(-av**2-2*av**2*m/y-y**2)*w**2-lam**2*y**4*w**2/D(y)
    Fp=mp.diff(FK,x); F=FK(x); d=D(x)
    PK=-lam*x*w**2*Fp/(2*F**2)+2*av*lam*m*w/(x*d*F)+2*lam*(-av**2*m*x**2+x**4*(-2*m+x))*w**2/(x**3*d*F)
    QK=-lam*x*w*Fp/(2*F**2)+2*lam*m*(av**2+x**2)*w/(x*d*F)-2*av*lam*m*(av**2+3*x**2)*w**2/(x*d*F)
    BK=lambda y: lam*y*w/mp.sqrt(FK(y)); Bp=mp.diff(BK,x)
    Xi=m*d/(x**4*F)-2*av*m*d*w/(x**4*F)+lam**2*x*(av**2-m*x)*w**2/(d*F)-(-av**2*m*x**2+x**5)*d*w**2/(x**6*F)+Bp**2
    return PK,QK,Xi,d
def true_acc(x,av,w,lam):
    gf,Gf,_=FUN['BL']; cf=cspiral('BL',av,w,lam,ref='BL')
    g=gf(x,av); G=Gf(x,av)
    def u(y):
        gg=gf(y,av); c=cf(y); N=mp.sqrt(sum(gg[i,j]*c[i]*c[j] for i in R3 for j in R3)); return [q/N for q in c]
    uu=u(x); du=[mp.diff(lambda z:u(z)[m],x) for m in R3]
    return [uu[1]*du[m]+sum(G[m][n][l]*uu[n]*uu[l] for n in R3 for l in R3) for m in R3]
for (a_,w_,l_,x_) in [('0.5','0.1','-0.45','3.0'),('1.0','0.09','-0.2','2.5')]:
    av=mp.mpf(a_); w=mp.mpf(w_); lam=mp.mpf(l_); x=mp.mpf(x_)
    PK,QK,Xi,d=parts(x,av,w,lam); ac=true_acc(x,av,w,lam)
    print('a=%s: true a^t=%s a^r=%s a^phi=%s'%(a_,mp.nstr(ac[0],8),mp.nstr(ac[1],8),mp.nstr(ac[2],8)))
    print('      appendix P_K=%s Q_K=%s  (r^2/Delta)Xi=%s vs (a^r)^2 r^2/Delta=%s'%(mp.nstr(PK,8),mp.nstr(QK,8),mp.nstr(x**2/d*Xi,8),mp.nstr(x**2/d*ac[1]**2,8)))
