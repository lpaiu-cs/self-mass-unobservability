# Lemma 64: Local Finite-State Construction

## Statement

- Status: Proven. Keep the primitive-family envelope, the fixed-order local operator catalog, and the reduced scalar normal-form basis `Y^I` already established in the positive theorem domain.
- Status: Proven. Modify only the monopole bookkeeping layer by adjoining one explicit local state variable `\chi(\tau)`.
- Status: Proven. Then there are two canonical one-state control models for the `A4` boundary:
  an adiabatic or slaved control,
  and a genuinely dynamical one-state local model.

## Explicit Models

- Status: Proven. The adiabatic or slaved control is

```math
m_A(Y,\chi)=m_0+\alpha Y+\lambda \chi,
\qquad
\chi = \beta Y,
```

or equivalently the limit of the dynamical model in which `\varepsilon_\chi := \Omega_{\rm orb} T_h \ll 1` and the homogeneous mode is absent.

- Status: Proven. The genuinely dynamical one-state local model is

```math
m_A(Y,\chi)=m_0+\alpha Y+\lambda \chi,
\qquad
\dot{\chi} = -\frac{1}{T_h}\chi + \frac{\beta}{T_h}Y.
```

## Why These Models Stay Inside The Intended Boundary Test

1. Status: Proven. Both models are local on the worldline and analytic in the variables that are kept explicitly.
2. Status: Proven. Both models use only the already admitted scalar normal-form coordinate `Y` and one explicit local state variable `\chi`; they introduce no new primitive family class.
3. Status: Proven. Both models leave the fixed-order local scalar operator catalog untouched, because the new ingredient is a state variable in the monopole sector rather than a new local derivative operator family.
4. Status: Proven. The slaved control is the right near-miss benchmark because it shows how a nominal state variable can still collapse back into Y-only bookkeeping.
5. Status: Proven. The dynamical one-state model is the smallest sharp `A4` test because it adds only one local state and one local first-order evolution law.

## Boundary

- Status: Note. This lemma constructs the local state models only.
- Status: Proven. The next step is to distinguish adiabatic eliminability from genuine orbital-timescale state dynamics.
