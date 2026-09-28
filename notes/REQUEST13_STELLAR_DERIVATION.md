# Specified-EOS scalar matching

## Equations and conventions

Status: Counterexample candidate. The calculation fixes massless DEF gravity with A(phi)=exp(beta phi^2/2), beta=-4 and zero background scalar, a nonrotating star, and the LAL SLy pressure/energy table. The target ADM/gravitational mass is 1.4378144085 solar masses; it is never passed as baryon mass to a fit. This choice specifies a conditional candidate; it is not an estimate of the actual EOS or gravity parameters of J0337.

Status: Proven. Write ds^2=-exp(2 nu)dt^2+dr^2/f+r^2 dOmega^2, f=1-2m/r, G=c=1. The zero-scalar branch obeys m'=4 pi r^2 epsilon, nu'=(m+4 pi r^3 p)/(r^2 f), p'=-(epsilon+p)nu'. At phi=0 the linear scalar equation decouples from fluid and metric perturbations because A'(0)=0. With time dependence exp(-i omega t), the regular radial field satisfies

    psi'' + (2/r + nu' - lambda') psi'
      + [omega^2 exp(-2 nu) + 4 pi beta(-epsilon+3p)] psi/f = 0,
    lambda' = (4 pi r epsilon - m/r^2)/f.

Status: Proven. At zero frequency the exterior solution is psi=phi_infinity-q log(f)/(2M). Thus q=-R^2 f(R) psi'(R) and phi_infinity=psi(R)+q log(f(R))/(2M). Their ratio is the static charge susceptibility in metres. The beta=0 solution is constant and has zero charge. The center expansion used to start integration is psi=1-B_c r^2/6+O(r^4), psi'=-B_c r/3+O(r^3), where B_c=omega^2 exp(-2nu_c)+4 pi beta(-epsilon_c+3p_c).

Status: Proven. Outside the star u=r psi obeys d^2u/dx^2+[omega^2-f 2M/r^3]u=0 with dx/dr=1/f. For real omega the incoming solution is the complex conjugate of the outgoing solution. The regular interior determines their coefficient ratio S=-A_out/A_in. The reported scattering convention removes the beta=0 propagation phase: H=(S_beta/S_0-1)/(2i omega). Merely subtracting S_0 retains a propagation phase; the first calculation with that convention is preserved and is not used as a local damping match. Elastic flux conservation gives |S_beta/S_0|=1 and hence Im(1/H)=-omega where H is nonzero. This identity concerns this explicitly defined one-channel scattering amplitude; it is not a proof that every effective internal coordinate has a frequency-independent damping coefficient equal to one.

Status: Counterexample candidate. Complex poles are found by matching the regular stellar logarithmic derivative to an outgoing exterior solution continued along r=R+exp(i theta)s. The logarithmic derivative avoids exponential amplitude overflow. Theta=1.8 and 2.0 radians and exterior lengths 2e6 and 4e6 metres are compared. A matched pole is not a count of all modes or a proof that no other unstable mode exists.

## Independent numerical checks

Status: Imported from prior work. The LAL table reader explicitly specifies pressure and energy density in m^-2. The initial monotone PCHIP interpolation and the native LAL pressure-energy interpolant give nearby results. However, the default native LAL TOV call at the same central pressure differs in mass by about 0.878%. This is not explained by the small difference between those two pressure-energy interpolants. Resampling the same native continuum at 4096 and 16384 pressures before the enthalpy-based LAL integration reduces the mass discrepancy to 4.36e-6 and 2.73e-7, respectively, and radius discrepancy to 1.34e-6 and 1.72e-7. This localizes the discrepancy to the sparse-table enthalpy path, rather than changing the EOS to obtain agreement. It is an empirical convergence check, not a rigorously bounded global error.

Status: Imported from prior work. The native-EOS direct integration gives R=11699.2283 m and susceptibility 41566.0079 m. The matched scalar pole is approximately omega=(3423.267-5013.329 i) per second, corresponding to a decay time about 0.1995 ms. All three tested exterior contour/radius choices find the same pole to about 2e-8 relative. The low-frequency radiative scale susceptibility/c is about 0.1387 ms. These results replace the previous absence of actual EOS and mode calculations for this specified candidate; they do not establish an orbital-timescale state.

Status: Conjectural. Extension to nonzero cosmological scalar, a scalarized equilibrium branch, white-dwarf structure, coupled-body radiation, and the full timing force/readout map remains necessary before this candidate yields a physical J0337 inference. A published prototype and this directly solved SLy candidate are distinct, and their coefficients must not be mixed.

Primary source for the prototype and matching definitions: [Khalil et al., PRD 106, 104016](https://arxiv.org/html/2206.13233v2), equations 28--31. Table and native numerical sources: https://git.ligo.org/lscsoft/lalsuite/-/tree/master/lalsimulation/lib (downloaded files and hashes retained). Calculation: verification/stellar_matching.py; raw results: outputs/research-remediation/stellar-matching.json. Classification: loophole progress for a specified physical candidate, with open empirical promotion gates.
