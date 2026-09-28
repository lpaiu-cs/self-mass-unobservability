import sys, time
import numpy as np
src = open('phase279-nongray.py', encoding='utf-8').read().split('# self-checks')[0]
exec(src)
Teff, lg = 15898.27, 5.73962
G = atmosphere(Teff, lg, gray=True)
N = atmosphere(Teff, lg, gray=False, T_init=G['T'], iters=600)
m, T, rho = N['m'], N['T'], N['rho']; H0 = SIG*Teff**4/(4*np.pi)
ab, sc = opacity(rho, T); J, Hh, B = solve_J(T, ab, sc, m)
Hn = np.r_[J[0:1]*0.5, 0.5*(Hh[1:] + Hh[:-1]), Hh[-1:]]; Ht = np.sum(Hn*w_nu, 1)
Hhalf = np.sum(Hh*w_nu, 1)
div = np.diff(Hhalf)/np.diff(0.5*(m[1:] + m[:-1]))  # flux divergence at interior nodes
kB = np.sum(ab*B*w_nu, 1); kJ = np.sum(ab*J*w_nu, 1)
for i in [0, 1, 5, 20, 40, 60, 80, 100, 120, 140, 160, 180, 199]:
    print('i %3d m %.2e tauR %.2e T %8.1f  H/H0-1 node %+.2e  %s' % (i, m[i], N['tauR'][i], T[i], Ht[i]/H0 - 1, ('half %+.2e  RE residual %+.2e' % (Hhalf[min(i, len(Hhalf)-1)]/H0 - 1, (kJ[i] - kB[i])/kB[i])) ))
