"""
Figure: kappa and tau_1 for two DIFFERENT worldlines with the same (a, omega, lambda):
  red : tangent (1, lambda*omega*r, 0, omega) imposed in Boyer-Lindquist (BL) components   [BL-form]
  blue: the same component formula imposed in the ingoing Kerr-Schild (KS) chart          [KS-form]
Equatorial Kerr, M = 1, signature (+,-,-,-).  Needs: sympy, mpmath, numpy, matplotlib.
Run:  python3 plot_two_worldlines.py        (writes fig_two_worldlines.png)
"""
import sympy as sp, mpmath as mp, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

mp.mp.dps = 30                      # high precision: tau_1 needs 2nd derivatives of the tangent

# ---------------- equatorial (t, r, phi) metrics -----------------
r, a = sp.symbols('r a')
Dl = r**2 + a**2 - 2*r              # Delta, with M = 1
gBL = sp.Matrix([[1 - 2/r, 0, 2*a/r],
                 [0, -r**2/Dl, 0],
                 [2*a/r, 0, -(r**2 + a**2 + 2*a**2/r)]])
# ingoing Kerr-Schild chart:  dt' = dt + (2r/Delta) dr,  dphi' = dphi + (a/Delta) dr
J = sp.eye(3); J[0, 1] = -2*r/Dl; J[2, 1] = -a/Dl      # maps new differentials to BL ones
gKS = sp.simplify(J.T*gBL*J)

def build(g):
    """lambdified metric and its analytic Christoffel symbols (only r-derivatives survive on the equator)"""
    gi = g.inv(); dg = g.diff(r)
    G = [[[sum(gi[m, s]*((dg[s, l] if n == 1 else 0) + (dg[s, n] if l == 1 else 0)
                         - (dg[n, l] if s == 1 else 0))/2 for s in range(3))
           for l in range(3)] for n in range(3)] for m in range(3)]
    return sp.lambdify((r, a), g, 'mpmath'), sp.lambdify((r, a), G, 'mpmath')

FUN = {'BL': build(gBL), 'KS': build(gKS)}
R3 = range(3)

# ---------------- curve and Frenet-Serret scalars -----------------
def tangent(av, w, lam):
    """unnormalised tangent (1, lam*w*r, w) in (t, r, phi); same formula, used in either chart"""
    return lambda x: [mp.mpf(1), x*lam*w, w]

def kin(chart, cf, x, av):
    gf, Gf = FUN[chart]
    g, G, c = gf(x, av), Gf(x, av), cf(x)
    S = sum(g[i, j]*c[i]*c[j] for i in R3 for j in R3)           # norm S = g(v,v)
    u = [ci/mp.sqrt(S) for ci in c]
    return g, G, u, S

def fs(chart, cf, x, av):
    """returns kappa^2, tau_1^2, S for the curve with tangent cf(r) in the given chart"""
    def acc(y):                                                    # A^mu = u^r d_r u^mu + Gamma u u
        g, G, u, S = kin(chart, cf, y, av)
        du = [mp.diff(lambda z: kin(chart, cf, z, av)[2][m], y) for m in R3]
        return [u[1]*du[m] + sum(G[m][n][l]*u[n]*u[l] for n in R3 for l in R3) for m in R3], g
    def e1(y):
        A, g = acc(y)
        k = mp.sqrt(-sum(g[i, j]*A[i]*A[j] for i in R3 for j in R3))
        return [q/k for q in A]
    g, G, u, S = kin(chart, cf, x, av)
    A, _ = acc(x)
    k2 = -sum(g[i, j]*A[i]*A[j] for i in R3 for j in R3); k = mp.sqrt(k2)
    e = e1(x)
    de = [mp.diff(lambda z: e1(z)[m], x) for m in R3]
    De = [u[1]*de[m] + sum(G[m][n][l]*u[n]*e[l] for n in R3 for l in R3) for m in R3]
    V = [De[m] - k*u[m] for m in R3]                               # = tau_1 e_2
    t2 = -sum(g[i, j]*V[i]*V[j] for i in R3 for j in R3)
    return k2, t2, S

# ---------------- parameters ----------------
av, w, lam = mp.mpf('0.5'), mp.mpf('0.1'), mp.mpf('-0.45')
rp = 1 + mp.sqrt(1 - av**2)                                        # outer horizon
cf = tangent(av, w, lam)

# radius where the BL-form curve becomes null (S_BL = 0), outside the horizon
def S_BL(x):
    g = FUN['BL'][0](x, av); c = cf(x)
    return sum(g[i, j]*c[i]*c[j] for i in R3 for j in R3)
rl = mp.findroot(S_BL, 2.0)
print('BL-form null radius r_l =', mp.nstr(rl, 8), ' ; r_+ =', mp.nstr(rp, 8))

def sweep(chart, radii):
    out = []
    for x in radii:
        try:
            k2, t2, S = fs(chart, cf, mp.mpf(float(x)), av)
            if S > 0 and k2 > 0 and t2 > 0:
                out.append((x, float(mp.sqrt(k2)), float(mp.sqrt(t2))))
        except Exception:
            pass
    return np.array(out)

# dense grids (log-spaced close to the null radius of the BL-form curve)
r_BL = float(rl) + np.logspace(-4, np.log10(2.0), 200)
r_KS = np.linspace(3.2, 1.62, 200)
B = sweep('BL', r_BL)
K = sweep('KS', r_KS)

# ---------------- plot ----------------
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for i, (idx, lab) in enumerate([(1, r'$\kappa$'), (2, r'$\tau_1$')]):
    ax[i].semilogy(B[:, 0], B[:, idx], 'r-', label='BL-form tangent')
    ax[i].semilogy(K[:, 0], K[:, idx], 'b-', label='KS-form tangent')
    ax[i].axvline(float(rp), color='k', ls='--', label=r'$r_+$')
    ax[i].axvline(float(rl), color='r', ls=':', label=r'BL curve null ($r_\ell$)')
    ax[i].set_xlabel('r/M'); ax[i].set_ylabel(lab); ax[i].set_xlim(1.55, 4.0)
ax[0].legend(fontsize=8)
fig.suptitle(r'$a=0.5,\ \omega=0.1,\ \lambda=-0.45$: two different worldlines')
plt.tight_layout()
plt.savefig('fig_two_worldlines.png', dpi=150)
print('saved fig_two_worldlines.png')
