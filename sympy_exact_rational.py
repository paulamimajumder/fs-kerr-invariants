"""
Exact-rational Frenet-Serret computation for the two prescribed equatorial ingoing Kerr-Schild spirals
(M = 1, signature (+,-,-,-), chart (t', r, phi'), tangent v = (1, lambda*omega*r, omega)):

    case 1:  a = 1/2,  omega = 1/10, lambda = -9/20
    case 2:  a = 3/2,  omega = 3/10, lambda = -9/20

Route A : tau_1^2 from the exact identities in terms of S = v.v, q = D_v v, p = D_v q  (no square roots):
              kappa^2 = Q^2/S^3 - q.q/S^2,   tau_1^2 = (kappa^4 - kappa_dot^2 - W.W/S^3)/kappa^2
Route B : tau_1^2 from Eq. (12)  (kappa^4 - kappa_dot^2 - B.B)/kappa^2  with A = Du/dtau, B = DA/dtau built from the
          normalised u = v/sqrt(S) (square roots kept); checked at exact rational radii to 40 digits.
Output : tau_1^2(r) as an exact rational function, the integer polynomials P, R with tau_1 = c*P/(r^2 R),
         the null radius, Q_0, lim kappa*S^(3/2), lim tau_1.
Requires: sympy, mpmath.        Run: python3 sympy_exact_rational.py
"""
import sympy as sp

r = sp.symbols('r', positive=True)

def run(a, w, lam, label, r_guess, radii):
    print('=' * 78); print(label, ' a =', a, ' omega =', w, ' lambda =', lam)
    D = r**2 + a**2 - 2*r
    # equatorial BL block (t, r, phi), Sigma = r^2
    gBL = sp.Matrix([[1 - 2/r, 0, 2*a/r],
                     [0, -r**2/D, 0],
                     [2*a/r, 0, -(r**2 + a**2 + 2*a**2/r)]])
    # ingoing Kerr-Schild chart: dt' = dt + (2r/D) dr, dphi' = dphi + (a/D) dr
    J = sp.eye(3); J[0, 1] = -2*r/D; J[2, 1] = -a/D
    g = (J.T*gBL*J).applyfunc(sp.cancel)
    gi = g.inv().applyfunc(sp.cancel)
    dg = g.diff(r)                                   # only r-derivatives survive on the equator
    Gam = [[[sp.cancel(sum(gi[m, s]*((dg[s, l] if n == 1 else 0) + (dg[s, n] if l == 1 else 0)
                                     - (dg[n, l] if s == 1 else 0))/2 for s in range(3)))
             for l in range(3)] for n in range(3)] for m in range(3)]
    R3 = range(3)
    dot = lambda x, y: sp.cancel(sum(g[i, j]*x[i]*y[j] for i in R3 for j in R3))

    v = [sp.Integer(1), lam*w*r, w]
    dds = lambda f: sp.cancel(lam*w*r*sp.diff(f, r))                       # d/ds along the curve (s = t')
    Dv = lambda X: [sp.cancel(dds(X[m]) + sum(Gam[m][n][l]*v[n]*X[l] for n in R3 for l in R3)) for m in R3]

    S = dot(v, v)
    q = Dv(v); p = Dv(q)
    Q = dot(q, v); qq = dot(q, q); pp = dot(p, p)
    assert sp.cancel(Q - dds(S)/2) == 0                                     # Q = S'/2
    kap2 = sp.cancel(Q**2/S**3 - qq/S**2)                                   # exact identity (Eq. kappa^2)
    C = sp.cancel(dds(Q)/S - 4*Q**2/S**2)
    WW = sp.cancel(pp + 9*Q**2*qq/S**2 + C**2*S - 3*Q*dds(qq)/S - 2*C*(dds(Q) - qq) + 6*C*Q**2/S)
    kd2 = sp.cancel(dds(kap2)**2/(4*kap2*S))                                # kappa_dot^2
    tauA = sp.cancel((kap2**2 - kd2 - WW/S**3)/kap2)                        # Route A: tau_1^2(r)

    num, den = sp.fraction(tauA)
    print('S(r)      =', sp.factor(S))
    print('tau_1^2 numerator degree', sp.degree(num, r), ', denominator degree', sp.degree(den, r))

    # ---- tau_1 = c * P / (r^2 R) : take exact square roots of numerator and denominator
    def sqrt_poly(expr):
        c, fl = sp.factor_list(expr)
        out = sp.Integer(1); sq = sp.sqrt(c)
        for f, e in fl:
            if e % 2: return None
            out *= f**(e//2)
        return sq, sp.expand(out)
    sn, sd = sqrt_poly(num), sqrt_poly(den)
    if sn and sd:
        cn, Pn = sn; cd, Pd = sd
        coef = sp.nsimplify(cn/cd)
        # remove powers of r from the denominator polynomial: den = (r^2 * R)^2 * const
        rp2 = sp.Poly(Pd, r)
        k = min(m[0] for m in rp2.monoms())
        R = sp.expand(sp.cancel(Pd/r**k)); P = Pn
        # normalise to primitive integer polynomials
        cP, P = sp.Poly(P, r).clear_denoms()[0], sp.Poly(P, r).clear_denoms()[1].as_expr()
        gP = sp.gcd_list(sp.Poly(P, r).all_coeffs()); P = sp.expand(P/gP)
        cR, Rp = sp.Poly(R, r).clear_denoms(); R = Rp.as_expr(); gR = sp.gcd_list(sp.Poly(R, r).all_coeffs()); R = sp.expand(R/gR)
        scale = sp.nsimplify(sp.simplify(sp.sqrt(sp.cancel(tauA*(r**2*R)**2/P**2))))
        print('tau_1 = %s * P(r) / (r^%d * R(r)),  with the primitive integer polynomials below' % (scale, k))
        print('P coefficients (descending powers of r):', sp.Poly(P, r).all_coeffs())
        print('R coefficients (descending powers of r):', sp.Poly(R, r).all_coeffs())
        # exact identity check
        assert sp.cancel(tauA - (scale*P/(r**k*R))**2) == 0
        print('exact identity tau_1^2 = (scale*P/(r^k R))^2 verified symbolically')
    else:
        print('numerator/denominator not perfect squares; using the rational function directly')

    # ---- Route B: Eq. (12) with u = v/sqrt(S), square roots kept; evaluated at exact rational radii (no symbolic simplification)
    sq = sp.sqrt(S); u = [vi/sq for vi in v]
    dotn = lambda x, y: sum(g[i, j]*x[i]*y[j] for i in R3 for j in R3)
    Dtau = lambda X: [u[1]*sp.diff(X[m], r) + sum(Gam[m][n][l]*u[n]*X[l] for n in R3 for l in R3) for m in R3]
    A = Dtau(u); Bv = Dtau(A)
    kap2B = -dotn(A, A); AB = dotn(A, Bv); BB = dotn(Bv, Bv)
    tauB = (kap2B**2 - AB**2/kap2B - BB)/kap2B
    import mpmath as mp; mp.mp.dps = 40
    fB = sp.lambdify(r, tauB, 'mpmath'); fA = sp.lambdify(r, tauA, 'mpmath')
    for rv in radii:   # radii must lie on the timelike side (S > 0) of the null point
        x = mp.mpf(rv.p)/rv.q
        a_, b_ = fA(x), fB(x)
        print('  r=%s: Route A = %s | Route B = %s | rel.diff = %s' % (rv, mp.nstr(a_, 22), mp.nstr(b_, 22), mp.nstr(abs(a_-b_)/abs(a_), 3)))

    # ---- null endpoint, Q0, limits
    Snum = sp.fraction(sp.cancel(S))[0]
    roots = [rt for rt in sp.Poly(Snum, r).nroots(n=30) if abs(sp.im(rt)) < 1e-20]
    rn = min(roots, key=lambda z: abs(sp.re(z) - r_guess)); rn = sp.re(rn)
    Q0 = Q.subs(r, rn).evalf(25)
    print('r_null =', sp.N(rn, 25))
    print('Q0 =', Q0, '   |Q0| = lim kappa*S^(3/2) =', abs(Q0))
    print('lim tau_1 =', sp.N(sp.sqrt(tauA.subs(r, rn)), 25))
    return tauA

import sys
if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'both'
    if which in ('1', 'both'):
        run(sp.Rational(1, 2), sp.Rational(1, 10), sp.Rational(-9, 20), 'case 1 (a/M = 1/2)', 1.6, [sp.Rational(5, 2), sp.Rational(9, 4), sp.Rational(17, 8)])
    if which in ('2', 'both'):
        run(sp.Rational(3, 2), sp.Rational(3, 10), sp.Rational(-9, 20), 'case 2 (a/M = 3/2)', 2.09, [sp.Rational(2), sp.Rational(15, 8), sp.Rational(7, 4)])
