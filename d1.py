import sympy as sp
s=sp.symbols('sigma',positive=True)
n1,n2,n3,n4,n5=sp.symbols('n1 n2 n3 n4 n5',real=True)
B0,B1,B2,B3=sp.symbols('B0 B1 B2 B3',real=True)
D0,D1=sp.symbols('D0 D1',real=True)
nu=n1*s+n2*s**2/2+n3*s**3/6+n4*s**4/24+n5*s**5/120
Bf=B0+B1*s+B2*s**2/2+B3*s**3/6        # b.b along the curve
dl=D0+D1*s                              # d.d
p=sp.diff(nu,s); pp=sp.diff(nu,s,2)
k2=p**2/(4*nu**3)-Bf/nu**2
C=pp/(2*nu)-p**2/nu**2
Bp=sp.diff(Bf,s)
W2=dl+9*p**2*Bf/(4*nu**2)+C**2*nu-3*p*Bp/(2*nu)-2*C*(pp/2-Bf)+3*C*p**2/(2*nu)
adot2=W2/nu**3
kd2=(sp.diff(k2,s))**2/(4*k2*nu)          # (dkappa/dtau)^2
t2=(k2**2-kd2-adot2)/k2
t2=sp.together(t2)
lim=sp.limit(t2,s,0,'+')
print('tau1*^2 =',sp.simplify(lim))
