from __future__ import annotations

import sympy as sp


DIM = 3


def symmetric_symbol_matrix(prefix: str, dim: int = DIM) -> sp.Matrix:
    """Return a symmetric matrix with shared symbols above/below the diagonal."""
    entries: dict[tuple[int, int], sp.Symbol] = {}
    for i in range(dim):
        for j in range(i, dim):
            entries[(i, j)] = sp.symbols(f"{prefix}{i + 1}{j + 1}", real=True)
    return sp.Matrix(
        dim,
        dim,
        lambda i, j: entries[(i, j)] if i <= j else entries[(j, i)],
    )


def matrix_dot(left: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    return sum(left[i] * right[i] for i in range(left.rows))


def gradient(expr: sp.Expr, coords: sp.Matrix) -> sp.Matrix:
    return sp.Matrix([sp.diff(expr, coord) for coord in coords])


def hessian(expr: sp.Expr, coords: sp.Matrix) -> sp.Matrix:
    return sp.Matrix(
        coords.rows,
        coords.rows,
        lambda i, j: sp.diff(expr, coords[i], coords[j]),
    )


def build_symbolic_results() -> dict[str, sp.Expr | sp.Matrix]:
    X = sp.Matrix(sp.symbols("X1:4", real=True))
    Xdot = sp.Matrix(sp.symbols("V1:4", real=True))
    Xddot = sp.Matrix(sp.symbols("A1:4", real=True))

    M, lam = sp.symbols("M lambda", positive=True, real=True)
    T_int, U_self = sp.symbols("T_int U_self", real=True)

    dipole = sp.Matrix(sp.symbols("d1:4", real=True))
    internal_momentum = sp.Matrix(sp.symbols("p1:4", real=True))
    I = symmetric_symbol_matrix("I")

    Phi = sp.Function("Phi")(*X)
    grad_phi = gradient(Phi, X)
    hess_phi = hessian(Phi, X)

    kinetic_general = sp.Rational(1, 2) * M * matrix_dot(Xdot, Xdot)
    kinetic_general += matrix_dot(Xdot, internal_momentum)
    kinetic_general += T_int

    potential_general = M * Phi + matrix_dot(dipole, grad_phi)
    potential_general += sp.Rational(1, 2) * sum(
        I[i, j] * hess_phi[i, j] for i in range(DIM) for j in range(DIM)
    )
    potential_general += (1 - lam) * U_self

    L_general = sp.expand(kinetic_general - potential_general)

    com_subs = {
        dipole[i]: sp.Integer(0) for i in range(DIM)
    } | {
        internal_momentum[i]: sp.Integer(0) for i in range(DIM)
    }
    L_com = sp.expand(L_general.subs(com_subs))
    L_int = T_int - (1 - lam) * U_self

    expected_L_com = (
        sp.Rational(1, 2) * M * matrix_dot(Xdot, Xdot)
        - M * Phi
        - sp.Rational(1, 2)
        * sum(I[i, j] * hess_phi[i, j] for i in range(DIM) for j in range(DIM))
        + L_int
    )

    euler_lagrange = sp.Matrix(
        [sp.simplify(M * Xddot[i] - sp.diff(L_com, X[i])) for i in range(DIM)]
    )
    expected_el = sp.Matrix(
        [
            sp.simplify(
                M * Xddot[i]
                + M * sp.diff(Phi, X[i])
                + sp.Rational(1, 2)
                * sum(
                    I[j, k] * sp.diff(Phi, X[i], X[j], X[k])
                    for j in range(DIM)
                    for k in range(DIM)
                )
            )
            for i in range(DIM)
        ]
    )
    acceleration_rhs = sp.Matrix(
        [
            sp.simplify(
                -sp.diff(Phi, X[i])
                - sp.Rational(1, 2)
                / M
                * sum(
                    I[j, k] * sp.diff(Phi, X[i], X[j], X[k])
                    for j in range(DIM)
                    for k in range(DIM)
                )
            )
            for i in range(DIM)
        ]
    )

    assert sp.simplify(L_com - expected_L_com) == 0
    for i in range(DIM):
        assert sp.simplify(euler_lagrange[i] - expected_el[i]) == 0
        assert not acceleration_rhs[i].has(lam)

    return {
        "X": X,
        "Xdot": Xdot,
        "Xddot": Xddot,
        "M": M,
        "lambda": lam,
        "T_int": T_int,
        "U_self": U_self,
        "dipole": dipole,
        "internal_momentum": internal_momentum,
        "I": I,
        "Phi": Phi,
        "L_general": L_general,
        "L_com": sp.expand(L_com),
        "L_int": L_int,
        "euler_lagrange": euler_lagrange,
        "acceleration_rhs": acceleration_rhs,
    }


def stf_quadrupole(second_moment: sp.Matrix) -> sp.Matrix:
    return sp.simplify(second_moment - sp.eye(second_moment.rows) * sp.trace(second_moment) / 3)


def uniform_sphere_results() -> dict[str, sp.Expr | sp.Matrix]:
    r, R, rho0 = sp.symbols("r R rho_0", positive=True, real=True)
    M_sphere = sp.symbols("M_sphere", positive=True, real=True)

    mass = sp.simplify(4 * sp.pi * sp.integrate(rho0 * r**2, (r, 0, R)))
    i_xx = sp.simplify(
        4 * sp.pi / 3 * sp.integrate(rho0 * r**4, (r, 0, R))
    )
    rho_sub = {rho0: sp.solve(sp.Eq(M_sphere, mass), rho0)[0]}
    second_moment = sp.diag(
        *([sp.simplify(i_xx.subs(rho_sub))] * DIM)
    )

    return {
        "mass": mass,
        "second_moment": second_moment,
        "quadrupole": stf_quadrupole(second_moment),
    }


def n1_polytrope_results() -> dict[str, sp.Expr | sp.Matrix]:
    r, R, rho_c = sp.symbols("r R rho_c", positive=True, real=True)
    M_poly = sp.symbols("M_poly", positive=True, real=True)

    rho = rho_c * sp.sin(sp.pi * r / R) / (sp.pi * r / R)
    mass = sp.simplify(4 * sp.pi * sp.integrate(rho * r**2, (r, 0, R)))
    i_xx = sp.simplify(4 * sp.pi / 3 * sp.integrate(rho * r**4, (r, 0, R)))
    rho_sub = {rho_c: sp.solve(sp.Eq(M_poly, mass), rho_c)[0]}
    second_moment = sp.diag(
        *([sp.simplify(i_xx.subs(rho_sub))] * DIM)
    )

    return {
        "density_profile": sp.simplify(rho),
        "mass": mass,
        "second_moment": second_moment,
        "quadrupole": stf_quadrupole(second_moment),
    }


def format_tensor(label: str, tensor: sp.Matrix) -> str:
    return f"{label} = {sp.sstr(sp.simplify(tensor))}"


def print_summary() -> None:
    symbolic = build_symbolic_results()
    sphere = uniform_sphere_results()
    poly = n1_polytrope_results()

    print("=== General COM-reduced Lagrangian ===")
    print(sp.sstr(symbolic["L_com"]))
    print()

    print("=== Euler-Lagrange equations for the COM coordinates ===")
    for i, eq in enumerate(symbolic["euler_lagrange"], start=1):
        print(f"EL_{i} = {sp.sstr(eq)}")
    print()

    print("=== COM acceleration law (lambda-free) ===")
    for i, acc in enumerate(symbolic["acceleration_rhs"], start=1):
        print(f"a_{i} = {sp.sstr(acc)}")
    print()

    print("=== Uniform sphere check ===")
    print(f"M_sphere = {sp.sstr(sphere['mass'])}")
    print(format_tensor("I_uniform", sphere["second_moment"]))
    print(format_tensor("Q_uniform", sphere["quadrupole"]))
    print()

    print("=== n=1 polytrope check ===")
    print(f"rho_n1(r) = {sp.sstr(poly['density_profile'])}")
    print(f"M_poly = {sp.sstr(poly['mass'])}")
    print(format_tensor("I_n1", poly["second_moment"]))
    print(format_tensor("Q_n1", poly["quadrupole"]))


if __name__ == "__main__":
    print_summary()
