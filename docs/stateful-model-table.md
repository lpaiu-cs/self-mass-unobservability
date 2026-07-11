# Stateful Model Table

| Model class | Locality | Analyticity | No-state theorem valid? | Finite state-augmented collapse valid? | Bookkeeping replacement data | Theorem layer broken |
| --- | --- | --- | --- | --- | --- | --- |
| Local analytic no-state control `m_A(Y)=m_0+\alpha Y+\beta Y^2` | Yes | Yes | Yes | Yes | Ordinary finite Taylor sensitivities in `Y^I` | None |
| Adiabatic/slaved local-state control `m_A(Y,\chi)=m_0+\alpha Y+\lambda \chi`, `\chi=\beta Y` or `\varepsilon_\chi \ll 1` with no homogeneous mode | Yes | Yes | Yes after elimination | Yes | Y-only sensitivities after eliminating `\chi`, or equivalently finite augmented coordinates with a slaved state | None; not a sharp `A4` escape |
| Genuinely dynamical one-state model `m_A(Y,\chi)=m_0+\alpha Y+\lambda \chi`, `\dot{\chi}=-(1/T_h)\chi+(\beta/T_h)Y` with independent orbital-timescale state data | Yes | Yes | No | Yes | Finite augmented state-space data `(Y^I,\chi^a)` plus local analytic evolution-law data, kept separate from Wilson coefficients | Original no-state reduction and Y-only Lemma 55 fail |

## Readout

- Status: Proven. The adiabatic/slaved control does not count as the sharp `A4` boundary because it collapses back into the Y-only theorem branch.
- Status: Proven. The smallest sharp `A4` counterexample is the one-state genuinely dynamical local model.
- Status: Proven. The same model also supports a finite state-augmented positive branch once the local state variable is kept explicitly.
