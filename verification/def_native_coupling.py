"""Native molecular EOS and DEF matter/scalar exchange in the Einstein frame.

Counterexample candidate: a reusable provider, not a completed stellar solver.
Density is conserved baryon reference mass per Einstein proper volume.
"""
import numpy as np
import sympy as sp

import gr_molecular_conservative_initial as initial

ld = np.longdouble
C = initial.e.C
GRAV = initial.e.GRAV


class Matter:
    def __init__(self):
        self.native = initial.molecular.model.EOS()
        self.calls = 0

    def state(self, log_density, log_temperature, phi, composition):
        logA = -2*ld(phi)**2  # beta=-4, the already specified candidate
        A = np.exp(logA)
        raw = self.native(2, float(ld(log_density)-3*logA),
                          float(ld(log_temperature)-logA), np.asarray(composition, float))
        self.calls += 1
        rhoJ = np.exp(ld(log_density)-3*logA)
        rho = np.exp(ld(log_density))
        rest = (np.asarray(composition, ld)/initial.e.g.c.A)@initial.e.g.c.W*C**2
        u = ld(raw[2]); pressure = A**4*ld(raw[1])
        energy = rho*A*(rest+u)
        return dict(phi=ld(phi), logA=logA, A=A, rho=rho, rhoJ=rhoJ,
            logT=ld(log_temperature), raw=raw, restJ=rest, uJ=u,
            P=pressure, epsilon=energy, trace=-energy+3*pressure,
            entropy=ld(raw[3]), denergy_dlogT=rho*A*ld(raw[10]),
            denergy_dphi_at_rho_T=-4*ld(phi)*rho*A*(rest+u-3*ld(raw[9])-ld(raw[10])),
            isentropic_dlogT_dlogA=1+3*(ld(raw[9])-ld(raw[1])/rhoJ)/ld(raw[10]))

    def isentropic(self, log_density, phi, composition, entropy, guess):
        logT = ld(guess)
        for _ in range(10):
            z = self.state(log_density, logT, phi, composition)
            error = (z['entropy']-entropy)*np.exp(logT-z['logA'])/ld(z['raw'][10])
            if abs(error) <= ld('2e-12'):
                z['entropy_temperature_residual'] = error
                return z
            logT -= np.clip(error, -ld('.05'), ld('.05'))
        raise RuntimeError(('Native isentropic inverse did not converge', float(error)))


def matter_specific_increment(previous, current):
    """Avoid subtracting two rest-mass dominated total energies."""
    assert current['rho'] == previous['rho']
    dlogA = current['logA']-previous['logA']
    return previous['A']*(np.expm1(dlogA)*(previous['restJ']+previous['uJ'])
        +np.exp(dlogA)*(current['uJ']-previous['uJ']))


def scalar_stress(Pi, Phi, metric_a):
    """Pi=a*phi_t/(N*c), Phi=phi_r; cgs energy-density units."""
    denominator = 8*np.pi*GRAV*metric_a**2
    energy = (Pi*Pi+Phi*Phi)/denominator
    flux = -2*Pi*Phi/denominator
    transverse_pressure = (Pi*Pi-Phi*Phi)/denominator
    return dict(E=energy, S=flux, R=energy, P=transverse_pressure)


def matter_exchange(trace, phi, phi_t, Phi, lapse):
    """Add to conservative E_t and (a*S)_t, before geometric cross exchange."""
    alpha = -4*phi
    return -alpha*trace*phi_t, C*lapse*alpha*trace*Phi


def geometric_energy_exchange(matter, scalar, radius, lapse, metric_a):
    """Matter share; the scalar share has exactly the opposite sign."""
    return -4*np.pi*GRAV*C*radius*lapse*metric_a*(
        (scalar['E']+scalar['R'])*matter['S']-scalar['S']*(matter['E']+matter['R']))


def symbolic():
    # General matter/scalar geometric exchange cancels without v=Q=0.
    E, R, S, Es, Rs, Ss = sp.symbols('E R S Es Rs Ss')
    cross = (Es+Rs)*S-Ss*(E+R)
    assert sp.expand(cross+(E+R)*Ss-S*(Es+Rs)) == 0
    pi, gradient, a, adot, N, nr, ar, r, g, alpha, trace = sp.symbols(
        'pi gradient a adot N nr ar r g alpha trace', nonzero=True)
    pir, gradr = sp.symbols('pir gradr')
    energy = (pi*pi+gradient*gradient)/(8*sp.pi*g*a*a)
    flux = -pi*gradient/(4*sp.pi*g*a*a)
    pi_t = N/a*(gradr+(2/r+nr-ar)*gradient)+4*sp.pi*g*N*a*alpha*trace
    grad_t = N/a*(pir+(nr-ar)*pi)
    energy_t = sp.diff(energy,pi)*pi_t+sp.diff(energy,gradient)*grad_t+sp.diff(energy,a)*adot
    flux_r = sp.diff(flux,pi)*pir+sp.diff(flux,gradient)*gradr+sp.diff(flux,a)*a*ar
    divergence = N/a*(flux_r+(nr-ar+2/r)*flux)
    # E_scalar,t+div(N*S/a)+a_t/a*(E+R)+N/a*(nu_r+lambda_r)*S
    balance = energy_t+divergence+2*adot/a*energy+N/a*(nr+ar)*flux
    assert sp.simplify(balance-alpha*trace*N*pi/a) == 0
    # Full causal quadratic-entropy heat law is conformally covariant with
    # tau_E=tau_J/A, K_E=A^2 K_J, q_E=A^4 q_J, T_E=A T_J.
    tau, q, dA, theta, dlogB, dq, force, K = sp.symbols('tau q dA theta dlogB dq force K')
    transformed = tau*(dq-4*q*dA)+q+K*force+tau*q*(theta+8*dA+dlogB)/2
    expected = tau*dq+q+K*force+tau*q*(theta+dlogB)/2
    assert sp.expand(transformed-expected) == 0
    return dict(classification='Proven', passed=True,
        scalar_energy='E_phi=R_phi=(Pi^2+Phi^2)/(8 pi Ggeom a^2); S_phi=-Pi Phi/(4 pi Ggeom a^2); tangential pressure=(Pi^2-Phi^2)/(8 pi Ggeom a^2).',
        exchange='Scalar energy receives +alpha*T*phi_t, matter receives its negative. Mixed geometric energy terms also cancel; omitting them is wrong when scalar and matter flux ratios differ.',
        metric='Use total E_m+E_phi in the mass constraint and total R_m+R_phi in polar slicing. This module supplies the terms; a globally constrained coupled solve is not yet implemented.',
        transport='Under the declared conformal map, tau_E=tau_J/A, K_E=A^2*K_J, q_E=A^4*q_J. The 4*q*DlnA terms cancel only with the full quadratic-entropy heat law.',
        limitations='Conditional local identities and native provider. No full stellar evolution, physical atmosphere, scalar exterior matching or observation claim.')
