import sympy as sp, mpmath as mp
mp.mp.dps=30
r,a,w=sp.symbols('r a omega',positive=True)
D=r**2+a**2-2*r
Nk=sp.sqrt(D)*(r**3*w**2-(1-a*w)**2)       # kappa = |Nk|/(r^3 A^2)
Nt=a*(1+w**2*(a**2+3*r**2))-2*a**2*w+w*r**2*(r-3)   # tau1 = |Nt|/(r^3 A^2)
A2=1-(r**2+a**2)*w**2-2*(1-a*w)**2/r
# (i) horizon, omega=omega_H
for av in [0.3,0.5,0.9]:
    rp=1+mp.sqrt(1-av**2); wH=av/(rp**2+av**2)
    print('a=%s: Nt(r+,wH)=%s'%(av,mp.nstr(Nt.subs({a:av,r:rp,w:wH}).evalf(25),6)))
# (ii) light surface: ratio |Nk|/|Nt| at A^2=0
for (av,wv) in [(0.5,0.15),(1.0,0.15),(1.5,0.25),(0.5,-0.15),(0.9,0.2)]:
    f=sp.lambdify(r,A2.subs({a:av,w:wv}),'mpmath')
    for guess in (1.5,2.0,3.0,4.0,6.0):
        try:
            rl=mp.findroot(f,guess)
            if rl>0.5 and abs(f(rl))<1e-20:
                ratio=abs(Nk.subs({a:av,w:wv,r:rl}).evalf(20)/Nt.subs({a:av,w:wv,r:rl}).evalf(20))
                print('a=%s w=%s light surface r=%s |Nk|/|Nt|=%s'%(av,wv,mp.nstr(rl,8),mp.nstr(ratio,10)));break
        except Exception: pass
# symbolic: is Nk^2 - Nt^2 proportional to A^2 ?
q=sp.simplify(sp.expand(Nk**2-Nt**2))
print('Nk^2-Nt^2 =',sp.factor(q))
print('divisible by (r*A^2)? remainder:',sp.simplify(sp.rem(sp.numer(sp.together(q)), sp.numer(sp.together(r*A2)),r)))
