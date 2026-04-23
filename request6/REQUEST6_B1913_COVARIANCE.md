# Request 6: B1913+16 Covariance-Aware Clock-Only Extension

## Goal

This stage implements the next step identified by the source scout:
add `PSR B1913+16` to the decoupled clock-only program with a
**covariance-aware non-gamma mass-inference path**.

The target here is limited and deliberate:

- strengthen the decoupled clock-only branch,
- do **not** claim a final tied-model verdict,
- carry the main nuisance that actually matters for this system:
  the Galactic correction to `\dot P_b`.

## Model

No published full PK covariance matrix was used. Instead, this stage uses a
hierarchical nuisance-aware likelihood over the latent set

- `P_b`
- `e`
- `\dot\omega`
- `\dot P_b^{obs}`
- `\Delta \dot P_{b,gal}`

with Gaussian priors taken from
[Weisberg, Nice, Taylor 2010](https://arxiv.org/abs/1011.0718).

The non-gamma mass-inference path is:

1. infer total mass from `\dot\omega`,
2. infer the unordered mass roots from
   `\dot P_b^{intr} = \dot P_b^{obs} - \Delta \dot P_{b,gal}`,
3. marginalize over the two branch assignments:
   - heavy star = timed pulsar,
   - light star = timed pulsar,
4. only then apply the corrected clock-sector observable

```math
\gamma_{\rm SMU} = \gamma_{\rm GR}\left[1+\frac{\eta_p}{1+X_c}\right].
```

This keeps `gamma` out of the mass reconstruction itself while still allowing
the `gamma` likelihood to select the physically relevant branch.

## Main Numerical Result

Running

```bash
python3 request6_b1913_covariance.py
```

currently gives:

```text
GR branch selection: heavy-pulsar = 1.000000

non-gamma predictive:
heavy-branch gamma_GR = 4.0863 ± 0.0989 ms
light-branch gamma_GR = 4.7334 ± 0.1042 ms

effective B1913 row after gamma branch selection:
sbar = 0.163236
sigma_delta = 2.435e-04

lever-arm audit update:
current |kappa_*|_95 = 8.416e+00
with B1913 = 4.877e-02
improvement = 172.55x
```

The key points are:

### 1. `B1913` sharply resolves the branch ambiguity

Once `gamma` is applied, the heavy-pulsar branch is effectively selected with
unit probability in this staging analysis. The resulting branch-selected masses
return to the familiar published values:

- `m_p = 1.439879 ± 0.000215 M_sun`
- `m_c = 1.388504 ± 0.000214 M_sun`

### 2. But `B1913` is not a clean “plug-in precision anchor”

The non-gamma predictive check is the important surprise.

Before `gamma` is applied, the `\dot\omega + \dot P_b` route predicts

- heavy-branch `gamma_GR = 4.0863 ± 0.0989 ms`

instead of the observed

- `gamma_obs = 4.2992 ± 0.0008 ms`.

So in this implementation, `gamma` is not merely making a tiny correction around
an already-sharp non-gamma prediction. It is reweighting the latent nuisance
space, especially the `\Delta \dot P_{b,gal}` tail, toward the region that makes
the branch-selected masses consistent with the observed Einstein delay.

This is exactly why `B1913` helps the clock-only branch but still should not be
marketed as a clean tied-model discriminator.

### 3. The source-scout conclusion survives the more careful fit

The scout predicted that `B1913` should behave like a very strong positive-side
precision source, with an effective `sigma_delta` at roughly the `2e-4` level.

The nuisance-aware implementation finds:

- effective `sigma_delta(B1913) = 2.435e-04`

which is close to the scout estimate.

Likewise, the scout expected `B1913` to push the slope audit to the
few-`10^-2` level rather than to the true `1e-2` target. The present fit gives

- updated `|kappa_*|_95 = 4.877e-02`

which is consistent with that expectation.

## Interpretation

This stage strengthens the project in a useful but limited way.

- It confirms that `B1913` is worth adding.
- It confirms that `B1913` alone is still **not** enough to make Request 6
  decisive.
- It makes the current bottleneck even clearer: the missing ingredient is still
  a low-side source with explicit nuisance control.

So the project logic remains:

1. `B1913` strengthens the decoupled clock-only branch.
2. The next essential upgrade is a low-side source:
   - `J1141-6545` if nuisance cleanliness is the priority,
   - `J1906+0746` if raw timing precision is the priority.
3. Only after that should the `kappa_*` audit be rerun and the role of Request 6
   in the overall novelty case be reconsidered.

## Why This Is Still Not The Tied Verdict

Even after this upgrade, the tied-model warning remains unchanged.

`B1913` uses relativistic orbital dynamics plus a Galactic `\dot P_b` correction
to infer masses. That is a valid and useful clock-side strengthening step, but
it is not the same thing as an internally closed joint free-fall-plus-clock
test. The final tied-vs-decoupled judgement still belongs to the combined
analysis with the free-fall sector.

## Files

- [request6_b1913_covariance.py](request6_b1913_covariance.py#L1)
- [request6_b1913_covariance_summary.json](request6_b1913_covariance_summary.json#L1)
- [request6_b1913_covariance_summary.svg](request6_b1913_covariance_summary.svg)
- [request6_b1913_covariance_posterior_b1913_only.tsv](request6_b1913_covariance_posterior_b1913_only.tsv#L1)
- [request6_b1913_covariance_posterior_combined.tsv](request6_b1913_covariance_posterior_combined.tsv#L1)
