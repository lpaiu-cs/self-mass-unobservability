# Request 10.7a: External Nutimo Handoff Packet

## Purpose

Status: Counterexample candidate. This packet turns the Request 10.7 runtime
pilot into an external collaborator contract. It states exactly what a
pre-provisioned `Nutimo` environment must return for the dynamic-chi branch to
move beyond symbolic runtime-worthiness.

Status: Proven. This packet is not a local runtime attempt. It does not build
`Nutimo`, run `Tempo2`, repair dependencies, or produce a `J0337` constraint.

## Target

Status: Imported from prior work. The dynamic-chi target is the shared finite
relaxation transfer

```text
G(z) = c_Y + beta/(1 + tau_chi z).
```

Status: Imported from prior work. The target carrier set is

```text
Omega_in, Omega_out, |Omega_in - Omega_out|.
```

Status: Imported from prior work. Request 10.7 says the external pilot must
run only these gates:

```text
configuration closure
-> named Jacobian rank gate
-> minimal synthetic shared-tau_chi injection
```

## Required Return Artifacts

Status: Counterexample candidate. A useful external pilot should return:

| artifact | purpose |
| --- | --- |
| `configuration_manifest.json` | active parameter blocks, delay flags, specialcase flag, fitted parameters, and target carrier extraction settings |
| `baseline_fit_summary.json` | code/data release identifiers, baseline parameters, residual summary |
| `finite_jacobian.npy` or `finite_jacobian.tsv` | residual derivative matrix for the standard finite fitted parameter set |
| `carrier_projection_rank.json` | rank and singular values for the six-real-coordinate three-carrier projection |
| `dynamic_chi_test_column.tsv` | synthetic shared-`tau_chi` residual column and carrier amplitudes |
| `synthetic_injection_recovery.json` | post-refit survival or collapse of the shared transfer relation |
| `decision_summary.md` | one-page verdict using the Request 10.7 labels |

## Hard Stop Rules

Status: Proven. Stop if the work becomes dependency repair rather than a
projection-class test.

Status: Proven. Stop if `RN_PL`-like harmonic amplitude/phase soak-up is
enabled on the target carriers.

Status: Proven. Stop if the target carrier inventory collapses below the
required three-sample set.

Status: Proven. Stop if the named carrier projection rank is `6/6`.

Status: Proven. Stop if the dynamic-chi test column lies in the finite
standard nuisance span.

## Promotion Rule

Status: Counterexample candidate. Promote to a real external runtime
experiment only if:

- configuration closure passes,
- carrier projection rank is `<=5/6`,
- the dynamic-chi test column is outside the standard finite nuisance span,
- and minimal synthetic shared-`tau_chi` injection survives standard finite
  nuisance refit.

Status: Proven. A pass is not a detection claim. It means the named
implementation is a scientifically meaningful runtime target for the
dynamic-chi positive branch.

Status: Proven. A fail is not an implementation failure. It is a classification
result: the current dynamic-chi route collapses for this named implementation
class.

## Machine Outputs

Status: Proven. The external handoff packet is written to

```text
outputs/tsv/dynamic_chi_external_nutimo_handoff_packet.tsv
outputs/json/dynamic_chi_external_nutimo_handoff_packet.json
```
