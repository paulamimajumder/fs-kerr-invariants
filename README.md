# fs-kerr-invariants

Code accompanying the paper *Frenet–Serret invariants of timelike worldlines in Kerr spacetime: chart independence,
horizon crossing and null-endpoint limits* (P. Majumder). It reproduces the tables, figures and exact results quoted in the paper.

Conventions: signature (+,−,−,−), G = c = M = 1; equatorial plane unless stated. Python ≥ 3.9.

```
pip install -r requirements.txt
python3 chart_independence_table.py     # run any script from this folder
```

Scripts use high-precision arithmetic (mpmath, 30–50 digits) because the torsion needs second derivatives of the tangent;
most run in seconds to a few minutes.

| Paper element | Script |
|---|---|
| Table II (equatorial rows: invariants of one worldline in BL, KS, IK, Doran-type charts) | `chart_independence_table.py` |
| Table II (off-equatorial rows, θ-dependent chart) and off-equatorial null-limit law | `d5.py` |
| Killing closed forms, Eqs. (21)–(24); light-surface identity | `c7.py`, `c8.py` |
| Table I (light surfaces); Kerr–Newman and extremal checks, Eqs. (26)–(28) | `c8.py`, `c6.py`, `d3.py` |
| Fig. 3 (BL-form vs KS-form worldlines) | `plot_two_worldlines.py` |
| Proposition 1 (derivation of the torsion limit, Appendix A) | `d1.py` (SymPy series expansion) |
| Table III (null-endpoint limits) and parametrisation-invariance check | `c2.py`, `d2.py`, `e1.py`, `e2.py` |
| Exact rational τ₁²(r), polynomials of Appendix B, Eq. (36) | `sympy_exact_rational.py` |
| Figs. 4, 8 (null-limit scaling; frame rates) | `d6.py` |
| Figs. 5–7 (curvature growth, torsion limits, timelike domain) | `f1.py` |
| Gyroscope rates β, Marck's rate, Table IV, accumulated angles | `c3.py`, `c4.py`, `c5.py`, `d3.py` |
| Spike near a near-geodesic point (a = 2.5 spiral) | `c5.py`, `c6.py` |
| Check of the Appendix A.1 expressions of arXiv:2605.29759v1 | `d4.py` |

`fast.py` is the common library: equatorial (t, r, φ) metrics in BL, Kerr–Schild, ingoing Kerr and Doran-type charts,
analytic Christoffel symbols, and κ², τ₁², S for any curve whose tangent depends on r only.

## Citation
See `CITATION.cff` (and the DOI of the archived release, once available).

## License
MIT (see `LICENSE`).
