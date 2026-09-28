"""Solver check: with a frequency-independent absorption opacity the non-gray Unsold-Lucy solver must reproduce the Eddington gray
relation T^4 = (3/4) Teff^4 (tau + 2/3) (same top boundary H = J/2)."""
import numpy as np
src = open('phase279-nongray.py', encoding='utf-8').read().split('# self-checks')[0]
exec(src)
def opacity(rho, T, kap=20.):  # override: constant absorption, no scattering
    return np.full((len(T), len(nu)), kap), np.zeros((len(T), len(nu)))
Teff = 15898.27
G = atmosphere(Teff, 5.73962, gray=True)
N = atmosphere(Teff, 5.73962, gray=False, T_init=G['T']*1.03)
Tedd = (0.75*Teff**4*(N['tauR'] + 2/3))**0.25; sel = N['tauR'] < 100
print('iterations', N['iterations'], 'final', N['final'], ' max |T/T_Edd - 1| (tau<100) = %.2e' % np.max(np.abs(N['T'][sel]/Tedd[sel] - 1)))
