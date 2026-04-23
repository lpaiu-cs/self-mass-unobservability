# Request 1: COM Decoupling Through Quadrupole Order

## Setup

Start from

```math
L=\int d^3\xi\,\rho(\boldsymbol\xi)\left[
\frac12(\dot{\mathbf X}+\dot{\boldsymbol\xi})^2
-\Phi_{\rm ext}(\mathbf X+\boldsymbol\xi)
-(1-\lambda)\Phi_{\rm self}(\boldsymbol\xi)
\right].
```

Use the COM gauge

```math
\int d^3\xi\,\rho\,\xi^i = 0,
\qquad
\int d^3\xi\,\rho\,\dot\xi^i = 0,
```

and define the internal second moment

```math
I^{ij}:=\int d^3\xi\,\rho\,\xi^i\xi^j.
```

## Quadrupole Expansion

Taylor-expand the external potential about the COM:

```math
\Phi_{\rm ext}(\mathbf X+\boldsymbol\xi)
=
\Phi_{\rm ext}(\mathbf X)
+\xi^i\partial_i\Phi_{\rm ext}(\mathbf X)
+\frac12\xi^i\xi^j\partial_i\partial_j\Phi_{\rm ext}(\mathbf X)
+O(\ell^3\nabla^3\Phi_{\rm ext}).
```

The kinetic term separates as

```math
\int d^3\xi\,\rho\,\frac12(\dot{\mathbf X}+\dot{\boldsymbol\xi})^2
=
\frac12 M\dot{\mathbf X}^2
+\dot X_i\int d^3\xi\,\rho\,\dot\xi^i
+T_{\rm int}
=
\frac12 M\dot{\mathbf X}^2+T_{\rm int}.
```

The potential term becomes

```math
\int d^3\xi\,\rho\,\Phi_{\rm ext}(\mathbf X+\boldsymbol\xi)
=
M\Phi_{\rm ext}(\mathbf X)
+\frac12 I^{ij}\partial_i\partial_j\Phi_{\rm ext}(\mathbf X)
+O(\ell^3\nabla^3\Phi_{\rm ext}),
```

because the dipole term vanishes in the COM gauge. Therefore

```math
L
=
\frac12 M\dot{\mathbf X}^2
-M\Phi_{\rm ext}(\mathbf X)
-\frac12 I^{ij}\partial_i\partial_j\Phi_{\rm ext}(\mathbf X)
+L_{\rm int}(\lambda)
+O(\ell^3\nabla^3\Phi_{\rm ext}),
```

with

```math
L_{\rm int}(\lambda)=T_{\rm int}-(1-\lambda)U_{\rm self}.
```

The crucial point is that `lambda` appears only inside `L_int`, not in the COM-dependent terms.

## COM Equation of Motion

Euler-Lagrange gives

```math
\frac{d}{dt}\frac{\partial L}{\partial \dot X^i}
-\frac{\partial L}{\partial X^i}=0
\quad\Longrightarrow\quad
M\ddot X^i
=
-M\partial_i\Phi_{\rm ext}
-\frac12 I^{jk}\partial_i\partial_j\partial_k\Phi_{\rm ext}
+O(\ell^3\nabla^4\Phi_{\rm ext}).
```

No `lambda` survives. Self-subtraction changes the internal energy bookkeeping, but does not generate a new COM force.

## Relation to the STF Quadrupole

Define the STF quadrupole moment

```math
Q^{ij}:=I^{ij}-\frac13\delta^{ij}I^k{}_k.
```

If the body sits in a vacuum region of the external field, then

```math
\nabla^2\Phi_{\rm ext}=0,
```

so the trace part of `I^{ij}` drops out and

```math
I^{ij}\partial_i\partial_j\Phi_{\rm ext}
=
Q^{ij}\partial_i\partial_j\Phi_{\rm ext}.
```

That is why the quadrupole-order COM Lagrangian may be written as

```math
L
=
\frac12 M\dot{\mathbf X}^2
-M\Phi_{\rm ext}
-\frac12 Q^{ij}\partial_i\partial_j\Phi_{\rm ext}
+L_{\rm int}(\lambda)
```

in the usual external-vacuum setting.

## Spherical Checks

For any spherically symmetric density profile `rho(r)`, one has

```math
I^{ij}=\frac{\delta^{ij}}{3}\,4\pi\int_0^R \rho(r)\,r^4\,dr,
\qquad
Q^{ij}=0.
```

### Uniform sphere

For `rho(r)=rho_0`,

```math
M=\frac{4\pi}{3}\rho_0 R^3,
\qquad
I^{ij}=\frac{MR^2}{5}\delta^{ij},
\qquad
Q^{ij}=0.
```

### Newtonian `n=1` polytrope

Using the Lane-Emden profile

```math
\rho(r)=\rho_c\frac{\sin(\pi r/R)}{\pi r/R},
```

one finds

```math
M=\frac{4\rho_c R^3}{\pi},
\qquad
I^{ij}=\frac{\pi^2-6}{3\pi^2}MR^2\,\delta^{ij},
\qquad
Q^{ij}=0.
```

So both benchmark spherical bodies have vanishing STF quadrupole and therefore follow the same test-body COM law through quadrupole order.

## Files

- Symbolic implementation: [request1_com_decoupling.py](request1_com_decoupling.py)
- This memo: [REQUEST1_COM_DECOUPLING.md](REQUEST1_COM_DECOUPLING.md)
