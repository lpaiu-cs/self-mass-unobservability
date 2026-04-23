# Request 2: Internal-Structure Consistency Check

## Setup

To test whether self-mass unobservability can be applied literally to stellar structure, replace the self-gravity strength by

```math
\alpha := 1-\lambda,
\qquad
0 \le \alpha \le 1,
```

inside the hydrostatic equations. At the Newtonian level this is equivalent to

```math
\frac{dP}{dr} = - \alpha G \frac{m(r)\rho(r)}{r^2},
\qquad
\frac{dm}{dr} = 4\pi r^2 \rho(r).
```

So Request 2 asks whether finite, self-gravitating equilibria survive as `alpha -> 0`, i.e. `lambda -> 1`.

## 1. Uniform-Sphere Toy Model

Take a constant-density sphere with total mass `M` and imposed radius `R`. Then

```math
\rho = \frac{3M}{4\pi R^3},
\qquad
m(r)=\frac{4\pi}{3}\rho r^3.
```

Hydrostatic balance gives

```math
P(r)=\frac{2\pi}{3}\alpha G \rho^2 (R^2-r^2),
\qquad
P_c = \frac{3\alpha G M^2}{8\pi R^4}.
```

The gravitational binding energy and pressure integral are

```math
E_g = -\frac{3\alpha G M^2}{5R},
\qquad
\Pi := \int P\,dV = \frac{\alpha G M^2}{5R}.
```

Hence the virial relation is

```math
3\Pi + E_g = 0.
```

As `alpha -> 0`, one finds

```math
P_c \to 0,
\qquad
E_g \to 0,
\qquad
\Pi \to 0.
```

This toy does **not** predict `R(lambda)` because `R` is imposed by hand. What it does show is that once self-gravity is removed, the object ceases to be self-bound. Any finite-radius configuration would need external confinement or a non-gravitational surface stress.

## 2. Newtonian `n=1` Polytrope

Now use a real barotropic model,

```math
P=K\rho^2.
```

For `n=1`, the Lane-Emden solution remains

```math
\rho(r)=\rho_c \frac{\sin(\pi r/R)}{\pi r/R},
```

but the scale length changes because `G -> alpha G`. Keeping the EOS constant `K` and the total mass `M` fixed gives

```math
R(\alpha) = \sqrt{\frac{\pi K}{2\alpha G}} = R_0 \alpha^{-1/2},
```

```math
\rho_c(\alpha)
=
\frac{M\,\alpha^{3/2}G^{3/2}}{\sqrt{2\pi}\,K^{3/2}}
=
\rho_{c,0}\alpha^{3/2},
```

```math
P_c(\alpha)
=
K\rho_c^2
=
\frac{M^2\alpha^3 G^3}{2\pi K^2}
=
P_{c,0}\alpha^3,
```

```math
E_g(\alpha)
=
-\frac{3}{4}\frac{\alpha G M^2}{R}
=
E_{g,0}\alpha^{3/2},
\qquad
\Pi(\alpha) = -\frac13 E_g(\alpha).
```

The script also checks symbolically that the density profile satisfies

```math
\frac{dP}{dr} + \alpha G \frac{m(r)\rho(r)}{r^2}=0
```

once the radius relation above is imposed.

## 3. The `lambda -> 1` Limit

For the `n=1` family at fixed `M` and `K`,

```math
\alpha \to 0^+
\quad\Longrightarrow\quad
R \to \infty,
\quad
\rho_c \to 0,
\quad
P_c \to 0,
\quad
E_g \to 0.
```

So the star does not approach a new compact equilibrium. It decompactifies into a zero-density, zero-binding cloud. This is the precise obstruction: if self-gravity is removed literally from the internal structure equations, the model no longer supports a finite self-gravitating star.

The pathology is not specific to the fixed-`M` path. If one keeps `K` and `rho_c` fixed instead, the `n=1` radius still diverges while the mass blows up. Either way the `alpha -> 0` limit fails to define a sensible compact object.

## 4. Interpretation

Request 2 therefore confirms the conceptual point from the research note:

- Literal self-unobservability inside the star is inconsistent with self-gravitating stellar structure.
- The only viable interpretation is an EFT in which the anomaly lives in the passive-coupling sector, the clock sector, or both.
- In other words, `lambda`-type self-gravity removal is not a structural modification of the star; it can only be a modification of how an already-formed body couples to external gravity or external redshift.

## Outputs

- Analysis code: [request2_internal_structure.py](request2_internal_structure.py)
- This memo: [REQUEST2_INTERNAL_STRUCTURE.md](REQUEST2_INTERNAL_STRUCTURE.md)
- Scaling table: [request2_internal_structure_scalings.tsv](request2_internal_structure_scalings.tsv)
- Scaling figure: [request2_internal_structure_scalings.svg](request2_internal_structure_scalings.svg)
