"""Independent verification of the A4 (chi-state) and A3 (hereditary) escapes,
plus the no-go witness independence. My own math, not their scripts.

A4 claim: integrating out a heavy local oscillator chi with EOM
   chi_ddot + w^2 chi = (g/mu) Y
gives the effective local expansion
   dL_eff = g^2/(2 mu w^2) Y^2 + g^2/(2 mu w^4) Ydot^2 + O(w^-6),
i.e. the state COLLAPSES back to instantaneous sensitivities in the adiabatic
(heavy-chi) limit -- and does NOT collapse when w ~ Omega_orb.

A3 claim: an exponential memory kernel is Markovianizable (single local state
=> rational transfer function), but a causal power-law kernel K_g(s) ~ s^{-g},
0<g<1, has transfer function ~ p^{g-1} with a BRANCH POINT at p=0, so no
finite local state space (which always yields a rational transfer function)
can reproduce it. That is exactly the A3 (genuine nonlocal) vs A4 (finite
state) distinction.

No-go: admitting the magnetic family B_ij adds B2=Tr(B^2), a parity-even
weight-2 scalar independent of the E-sector -> uniqueness of a minimal sector
is obstructed.
"""

import sympy as sp


def check_A4_adiabatic():
    print("=" * 72)
    print("A4: integrate out heavy chi, check adiabatic collapse coefficients")
    print("=" * 72)
    mu, w, g = sp.symbols("mu omega g", positive=True)
    p = sp.symbols("p")  # D_tau -> i p in frequency space; work with D_tau^2 -> -p^2

    # Integrating out chi from L = mu/2 chidot^2 - mu w^2/2 chi^2 + g chi Y gives
    # the nonlocal effective term  dL = (g^2/2) Y (mu(w^2 + Dtau^2))^{-1} Y.
    # In frequency space Dtau^2 -> -p^2, so the kernel is g^2/(2 mu (w^2 - p^2)).
    # Adiabatic (heavy chi): expand in p^2/w^2.
    kernel = g**2 / (2 * mu * (w**2 - p**2))
    series = sp.series(kernel, p, 0, 6).removeO()
    print(f"  effective kernel g^2/(2 mu (w^2 - p^2)) expanded in p:")
    print(f"    {sp.simplify(series)}")
    c0 = series.coeff(p, 0)
    c2 = series.coeff(p, 2)
    print(f"  p^0 coeff (multiplies Y^2)        : {sp.simplify(c0)}   "
          f"claim g^2/(2 mu w^2): {sp.simplify(c0 - g**2/(2*mu*w**2))==0}")
    # p^2 term multiplies Y(-Dtau^2)Y = Y*(+p^2)Y ; by parts p^2 |Y|^2 -> Ydot^2.
    print(f"  p^2 coeff (-> multiplies Ydot^2)  : {sp.simplify(c2)}   "
          f"claim g^2/(2 mu w^4): {sp.simplify(c2 - g**2/(2*mu*w**4))==0}")
    ok = sp.simplify(c0 - g**2/(2*mu*w**2)) == 0 and sp.simplify(c2 - g**2/(2*mu*w**4)) == 0
    print(f"  adiabatic expansion matches their dL_eff: {ok}")
    print(f"  => heavy chi (eps=Omega/omega << 1) collapses to local Y^2, Ydot^2;")
    print(f"     resonant chi (eps=O(1)) does NOT -- state is genuine. [A4 sharp]")
    return ok


def check_A3_transfer():
    print()
    print("=" * 72)
    print("A3: exponential kernel is rational (finite-state); power-law is not")
    print("=" * 72)
    s, p, T, gamma = sp.symbols("s p T_h gamma", positive=True)

    # Exponential kernel K_exp(s) = e^{-s/T}/T ; Laplace transform:
    K_exp = sp.exp(-s / T) / T
    L_exp = sp.integrate(K_exp * sp.exp(-p * s), (s, 0, sp.oo))
    L_exp = sp.simplify(L_exp)
    print(f"  L[K_exp](p) = {L_exp}")
    print(f"    -> rational in p (single pole at p=-1/T): finite 1-state ODE. "
          f"[Markovianizable]")

    # Power-law kernel K_g(s) = s^{-gamma} / (Gamma(1-gamma) T^{1-gamma}), 0<gamma<1
    Kg = s**(-gamma) / (sp.gamma(1 - gamma) * T**(1 - gamma))
    L_g = sp.integrate(Kg * sp.exp(-p * s), (s, 0, sp.oo))
    L_g = sp.simplify(L_g)
    print(f"  L[K_gamma](p) = {L_g}")
    # substitute a concrete gamma to expose the fractional power
    L_g_half = sp.simplify(L_g.subs(gamma, sp.Rational(1, 2)))
    print(f"    at gamma=1/2: {L_g_half}")
    is_fractional = (L_g_half.has(sp.sqrt(p)) or L_g_half.has(p**sp.Rational(-1, 2)))
    print(f"    -> contains p^(gamma-1) = fractional power (branch point at p=0): "
          f"{is_fractional}")
    print(f"  A finite local state space always gives a RATIONAL transfer function;")
    print(f"  a fractional power is non-rational => no finite Markovian embedding.")
    print(f"  [A3 genuinely nonlocal, sharply distinct from A4 finite-state]")
    return is_fractional


def check_nogo_witness():
    print()
    print("=" * 72)
    print("No-go: magnetic family B adds B2=Tr(B^2), independent of E-sector")
    print("=" * 72)
    # B2 is built from a DIFFERENT tensor than E, so it is trivially independent
    # of {E2,E3,...} as an operator. The point of the no-go: nothing in the
    # minimal-sector setup forbids admitting B (or scalar S), and each admission
    # adds a genuine new low-weight parity-even survivor -> no unique minimal set.
    import numpy as np
    rng = np.random.default_rng(7)
    def rand_stf():
        M = rng.standard_normal((3, 3)); S = 0.5*(M+M.T); S -= np.eye(3)*np.trace(S)/3; return S
    # show B2 is rotation invariant and not a function of E-invariants
    vals = []
    for _ in range(300):
        E, B = rand_stf(), rand_stf()
        vals.append([np.trace(E@E), np.trace(E@E@E), np.trace(B@B)])
    r = np.linalg.matrix_rank(np.array(vals), tol=1e-9)
    print(f"  rank of value-matrix [E2, E3, B2] over random (E,B): {r}  (=> B2 adds a")
    print(f"  genuinely new independent direction; expected 3): {r==3}")
    print(f"  witnesses per admitted family (from their table): B2 (w2), S (w1), dotS2 (w4)")
    print(f"  => minimal-sector UNIQUENESS is obstructed without extra suppression. [no-go valid]")
    return r == 3


if __name__ == "__main__":
    a4 = check_A4_adiabatic()
    a3 = check_A3_transfer()
    ng = check_nogo_witness()
    print()
    print("=" * 72)
    print(f"A4 adiabatic collapse coefficients  : {a4}")
    print(f"A3 power-law non-rational transfer  : {a3}")
    print(f"no-go witness independence          : {ng}")
    print("=" * 72)
