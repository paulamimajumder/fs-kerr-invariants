import sympy as sp
r,a,w=sp.symbols('r a omega',positive=True)
D=r**2+a**2-2*r
g=sp.Matrix([[1-2/r,0,2*a/r],[0,-r**2/D,0],[2*a/r,0,-(r**2+a**2+2*a**2/r)]])
gi=g.inv(); dg=g.diff(r)
G=[[[sum(gi[m,s]*((dg[s,l] if n==1 else 0)+(dg[s,n] if l==1 else 0)-(dg[n,l] if s==1 else 0))/2 for s in range(3)) for l in range(3)] for n in range(3)] for m in range(3)]
f=sp.Symbol('f')            # f = A^2 = norm of (1,0,omega)
c=[1,0,w]
fexpr=sum(g[i,j]*c[i]*c[j] for i in range(3) for j in range(3))
dot=lambda x,y: sum(g[i,j]*x[i]*y[j] for i in range(3) for j in range(3))
# a^mu = Gamma u u with u=c/A  => a = (Gamma c c)/f
acc=[sp.together(sum(G[m][n][l]*c[n]*c[l] for n in range(3) for l in range(3))/fexpr) for m in range(3)]
k2=sp.factor(sp.together(-dot(acc,acc)))
print('kappa^2 =',k2)
# De1 = Gamma u e1 ; e1 = acc/kappa ; u=c/A ; so De1 = Gamma(c,acc)/(A kappa);  V=De1-kappa*u
# tau1^2 = -(V.V) = -[ (Gamma(c,acc))^2/(f k2) - 2 (Gamma(c,acc).c)/(f k2^(1/2)) *k2^(1/2)... ] compute symbolically with kappa as sqrt
kap=sp.sqrt(k2)
Gca=[sum(G[m][n][l]*c[n]*acc[l] for n in range(3) for l in range(3)) for m in range(3)]
V2=sp.together(-dot([Gca[m]/(sp.sqrt(fexpr)*kap)-kap*c[m]/sp.sqrt(fexpr) for m in range(3)],[Gca[m]/(sp.sqrt(fexpr)*kap)-kap*c[m]/sp.sqrt(fexpr) for m in range(3)]))
t2=sp.factor(sp.simplify(V2))
print('tau1^2 =',t2)
import pickle; pickle.dump((str(k2),str(t2),str(fexpr)),open('killing.pkl','wb'))
# numeric check vs earlier code: a=.5, w=.15, r=3 -> tau1^2 = 0.01496894829
print('check:',sp.N(t2.subs({a:sp.Rational(1,2),w:sp.Rational(15,100),r:3}),12),' expected 0.01496894829')
