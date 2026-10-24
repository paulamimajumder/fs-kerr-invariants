"""Table II (equatorial rows): FS invariants of one worldline evaluated in several charts.
The worldline is fixed by the tangent (1, lambda*omega*r, omega) in the Kerr-Schild chart and transformed to the
other charts with the Jacobian (cspiral(..., ref='KS')).  The last row uses the BL-form tangent (a different curve).
The off-equatorial rows are produced by d5.py."""
from fast import *
mp.mp.dps = 30
rows = [  # a, omega, lambda, r, charts
    ('0.5','0.1','-0.45','2.5',['BL','KS','IK','DO']),
    ('0.5','0.1','-0.45','1.8',['KS','IK','DO']),          # inside the horizon (r+ = 1.866)
    ('1.0','0.3','-0.45','1.2',['BL','KS','IK','DO']),
    ('1.5','0.3','-0.15','1.2',['BL','KS','IK','DO'])]      # a > M: no horizon
for a_, w_, l_, x_, charts in rows:
    av, w, lam, x = map(mp.mpf, (a_, w_, l_, x_))
    vals = {}
    for ch in charts:
        k2, t2, S = fs(ch, cspiral(ch, av, w, lam, ref='KS'), x, av)
        vals[ch] = (S, k2, t2)
    ref = vals[charts[0]]
    spread = max(abs(vals[c][i]-ref[i])/abs(ref[i]) for c in charts for i in range(3))
    print('a=%s w=%s lam=%s r=%s | S=%s kappa^2=%s tau1^2=%s | charts %s | max rel. spread %s' %
          (a_, w_, l_, x_, mp.nstr(ref[0], 15), mp.nstr(ref[1], 15), mp.nstr(ref[2], 15), charts, mp.nstr(spread, 3)))
av, w, lam, x = mp.mpf('0.5'), mp.mpf('0.1'), mp.mpf('-0.45'), mp.mpf('2.5')
k2, t2, S = fs('BL', cspiral('BL', av, w, lam, ref='BL'), x, av)
print('BL-form tangent (different curve): S=%s kappa^2=%s tau1^2=%s' % (mp.nstr(S, 10), mp.nstr(k2, 10), mp.nstr(t2, 10)))
