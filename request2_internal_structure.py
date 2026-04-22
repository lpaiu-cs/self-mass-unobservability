from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import sympy as sp


SVG_HEADER = """<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>
  .title {{ font: 700 24px sans-serif; fill: #111827; }}
  .subtitle {{ font: 400 14px sans-serif; fill: #4b5563; }}
  .label {{ font: 600 14px sans-serif; fill: #1f2937; }}
  .tick {{ font: 12px monospace; fill: #374151; }}
  .legend {{ font: 13px sans-serif; fill: #111827; }}
  .note {{ font: 12px sans-serif; fill: #4b5563; }}
  .axis {{ stroke: #111827; stroke-width: 1.4; fill: none; }}
  .grid {{ stroke: #d1d5db; stroke-width: 1; fill: none; }}
  .uniform {{ stroke: #d97706; stroke-width: 2.8; fill: none; stroke-dasharray: 8 5; }}
  .poly {{ stroke: #0f766e; stroke-width: 3.2; fill: none; }}
</style>
"""


@dataclass(frozen=True)
class UniformSphereToy:
    lam: np.ndarray
    alpha: np.ndarray
    radius_ratio: np.ndarray
    density_ratio: np.ndarray
    central_pressure_ratio: np.ndarray
    binding_energy_ratio: np.ndarray
    pressure_integral_ratio: np.ndarray


@dataclass(frozen=True)
class N1PolytropeFamily:
    lam: np.ndarray
    alpha: np.ndarray
    radius_ratio: np.ndarray
    central_density_ratio: np.ndarray
    central_pressure_ratio: np.ndarray
    binding_energy_ratio: np.ndarray
    pressure_integral_ratio: np.ndarray


def alpha_from_lambda(lam: np.ndarray | float) -> np.ndarray | float:
    return 1.0 - lam


def uniform_sphere_symbolic() -> dict[str, sp.Expr]:
    lam, G, M, R = sp.symbols("lambda G M R", positive=True, real=True)
    alpha = 1 - lam
    rho = 3 * M / (4 * sp.pi * R**3)
    pc = sp.simplify(2 * sp.pi * alpha * G * rho**2 * R**2 / 3)
    eg = sp.simplify(-3 * alpha * G * M**2 / (5 * R))
    pressure_integral = sp.simplify(alpha * G * M**2 / (5 * R))
    virial_residual = sp.simplify(3 * pressure_integral + eg)

    return {
        "alpha": alpha,
        "rho": rho,
        "Pc": pc,
        "Eg": eg,
        "Pi": pressure_integral,
        "virial_residual": virial_residual,
    }


def n1_polytrope_symbolic() -> dict[str, sp.Expr]:
    lam, G, M, K = sp.symbols("lambda G M K", positive=True, real=True)
    alpha = 1 - lam
    radius = sp.sqrt(sp.pi * K / (2 * alpha * G))
    rho_c = sp.simplify(M * alpha ** sp.Rational(3, 2) * G ** sp.Rational(3, 2) / (sp.sqrt(2 * sp.pi) * K ** sp.Rational(3, 2)))
    pc = sp.simplify(K * rho_c**2)
    eg = sp.simplify(-3 * alpha * G * M**2 / (4 * radius))
    pressure_integral = sp.simplify(-eg / 3)
    virial_residual = sp.simplify(3 * pressure_integral + eg)

    return {
        "alpha": alpha,
        "R": radius,
        "rho_c": rho_c,
        "Pc": pc,
        "Eg": eg,
        "Pi": pressure_integral,
        "virial_residual": virial_residual,
    }


def n1_hydrostatic_residual_symbolic() -> sp.Expr:
    r, R, K, rho_c, g_eff = sp.symbols("r R K rho_c g_eff", positive=True, real=True)
    x = sp.pi * r / R
    rho = sp.simplify(rho_c * sp.sin(x) / x)
    pressure = sp.simplify(K * rho**2)
    enclosed_mass = sp.simplify(4 * sp.pi * rho_c * R**3 * (sp.sin(x) - x * sp.cos(x)) / sp.pi**3)
    residual = sp.simplify(sp.diff(pressure, r) + g_eff * enclosed_mass * rho / r**2)
    residual = sp.simplify(residual.subs(g_eff, sp.pi * K / (2 * R**2)))
    return sp.trigsimp(sp.factor(residual))


def build_uniform_sphere_family(lam: np.ndarray) -> UniformSphereToy:
    alpha = alpha_from_lambda(lam)
    return UniformSphereToy(
        lam=lam,
        alpha=alpha,
        radius_ratio=np.ones_like(alpha),
        density_ratio=np.ones_like(alpha),
        central_pressure_ratio=alpha,
        binding_energy_ratio=alpha,
        pressure_integral_ratio=alpha,
    )


def build_n1_polytrope_family(lam: np.ndarray) -> N1PolytropeFamily:
    alpha = alpha_from_lambda(lam)
    return N1PolytropeFamily(
        lam=lam,
        alpha=alpha,
        radius_ratio=alpha ** (-0.5),
        central_density_ratio=alpha ** 1.5,
        central_pressure_ratio=alpha**3,
        binding_energy_ratio=alpha ** 1.5,
        pressure_integral_ratio=alpha ** 1.5,
    )


def write_scaling_table(path: str | Path, uniform: UniformSphereToy, poly: N1PolytropeFamily) -> Path:
    output_path = Path(path)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "lambda",
                "alpha",
                "uniform_radius_ratio",
                "uniform_density_ratio",
                "uniform_Pc_ratio",
                "uniform_absEg_ratio",
                "uniform_Pi_ratio",
                "n1_radius_ratio",
                "n1_rho_c_ratio",
                "n1_Pc_ratio",
                "n1_absEg_ratio",
                "n1_Pi_ratio",
            ]
        )
        for idx in range(len(uniform.lam)):
            writer.writerow(
                [
                    f"{uniform.lam[idx]:.6f}",
                    f"{uniform.alpha[idx]:.6e}",
                    f"{uniform.radius_ratio[idx]:.6e}",
                    f"{uniform.density_ratio[idx]:.6e}",
                    f"{uniform.central_pressure_ratio[idx]:.6e}",
                    f"{uniform.binding_energy_ratio[idx]:.6e}",
                    f"{uniform.pressure_integral_ratio[idx]:.6e}",
                    f"{poly.radius_ratio[idx]:.6e}",
                    f"{poly.central_density_ratio[idx]:.6e}",
                    f"{poly.central_pressure_ratio[idx]:.6e}",
                    f"{poly.binding_energy_ratio[idx]:.6e}",
                    f"{poly.pressure_integral_ratio[idx]:.6e}",
                ]
            )
    return output_path


def _polyline_path(xs: np.ndarray, ys: np.ndarray, x_map, y_map) -> str:
    points = [f"{x_map(float(x)):.2f},{y_map(float(y)):.2f}" for x, y in zip(xs, ys)]
    return "M " + " L ".join(points)


def _map_log_y(x0: float, x1: float, y0: float, y1: float, y_min: float, y_max: float):
    log_min = math.log10(y_min)
    log_max = math.log10(y_max)
    span_x = x1 - x0
    span_y = y1 - y0

    def xy(x: float, y: float) -> tuple[float, float]:
        y = max(y, y_min)
        log_y = math.log10(y)
        frac_y = (log_y - log_min) / (log_max - log_min)
        return x0 + span_x * x, y1 - span_y * frac_y

    def y_only(y: float) -> float:
        return xy(0.0, y)[1]

    return xy, y_only


def _panel_frame(parts: list[str], x0: float, y0: float, x1: float, y1: float, title: str) -> None:
    parts.append(f'<rect class="axis" x="{x0:.2f}" y="{y0:.2f}" width="{x1 - x0:.2f}" height="{y1 - y0:.2f}" />')
    parts.append(f'<text class="label" x="{x0:.2f}" y="{y0 - 10:.2f}">{title}</text>')


def _draw_x_ticks(parts: list[str], x0: float, y0: float, x1: float, y1: float) -> None:
    tick_values = [0.0, 0.25, 0.5, 0.75, 0.9, 0.99]
    for value in tick_values:
        x = x0 + (x1 - x0) * value
        parts.append(f'<line class="grid" x1="{x:.2f}" y1="{y0:.2f}" x2="{x:.2f}" y2="{y1:.2f}" />')
        parts.append(f'<text class="tick" x="{x - 12:.2f}" y="{y1 + 18:.2f}">{value:g}</text>')
    parts.append(f'<text class="label" x="{(x0 + x1) / 2 - 24:.2f}" y="{y1 + 38:.2f}">lambda</text>')


def _draw_log_y_ticks(parts: list[str], x0: float, x1: float, tick_values: list[float], y_map) -> None:
    for value in tick_values:
        y = y_map(value)
        parts.append(f'<line class="grid" x1="{x0:.2f}" y1="{y:.2f}" x2="{x1:.2f}" y2="{y:.2f}" />')
        label = f"1e{int(round(math.log10(value)))}" if value < 1 else f"{value:g}"
        parts.append(f'<text class="tick" x="{x0 - 44:.2f}" y="{y + 4:.2f}">{label}</text>')


def _render_scaling_svg(uniform: UniformSphereToy, poly: N1PolytropeFamily) -> str:
    width, height = 1320, 980
    parts = [SVG_HEADER.format(width=width, height=height)]
    parts.append('<rect width="100%" height="100%" fill="#f8fafc" />')
    parts.append('<text class="title" x="50" y="48">Request 2: Internal-Structure Consistency</text>')
    parts.append('<text class="subtitle" x="50" y="74">Scalings for alpha = 1 - lambda when self-gravity is inserted literally into hydrostatic support.</text>')

    legend_y = 108
    parts.append(f'<line class="uniform" x1="52" y1="{legend_y:.2f}" x2="122" y2="{legend_y:.2f}" />')
    parts.append(f'<text class="legend" x="132" y="{legend_y + 4:.2f}">Uniform sphere toy (fixed M, fixed input R)</text>')
    parts.append(f'<line class="poly" x1="520" y1="{legend_y:.2f}" x2="590" y2="{legend_y:.2f}" />')
    parts.append(f'<text class="legend" x="600" y="{legend_y + 4:.2f}">n = 1 polytrope (fixed M, fixed K)</text>')

    panels = [
        (60, 150, 620, 470, "Radius ratio R(lambda) / R(0)", "radius"),
        (700, 150, 1260, 470, "Central density ratio rho_c(lambda) / rho_c(0)", "density"),
        (60, 560, 620, 880, "Central pressure ratio P_c(lambda) / P_c(0)", "pressure"),
        (700, 560, 1260, 880, "Binding-energy ratio |E_g(lambda)| / |E_g(0)|", "binding"),
    ]

    xs = uniform.lam
    for x0, y0, x1, y1, title, key in panels:
        _panel_frame(parts, x0, y0, x1, y1, title)
        _draw_x_ticks(parts, x0, y0, x1, y1)

        if key == "radius":
            tick_values = [1, 2, 5, 10, 20, 50]
            xy_map, y_map = _map_log_y(x0, x1, y0, y1, 1, 50)
            uniform_y = uniform.radius_ratio
            poly_y = poly.radius_ratio
            note = "Uniform-sphere radius is an imposed toy input, not an equilibrium prediction."
        elif key == "density":
            tick_values = [1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
            xy_map, y_map = _map_log_y(x0, x1, y0, y1, 1e-5, 1)
            uniform_y = uniform.density_ratio
            poly_y = poly.central_density_ratio
            note = "The n=1 family must dilute as lambda approaches 1 to keep M and K fixed."
        elif key == "pressure":
            tick_values = [1, 1e-2, 1e-4, 1e-6, 1e-8]
            xy_map, y_map = _map_log_y(x0, x1, y0, y1, 1e-9, 1)
            uniform_y = uniform.central_pressure_ratio
            poly_y = poly.central_pressure_ratio
            note = "Both families lose pressure support; the n=1 barotrope does so much faster."
        else:
            tick_values = [1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
            xy_map, y_map = _map_log_y(x0, x1, y0, y1, 1e-5, 1)
            uniform_y = uniform.binding_energy_ratio
            poly_y = poly.binding_energy_ratio
            note = "The binding energy collapses to zero, so no self-bound compact object survives."

        _draw_log_y_ticks(parts, x0, x1, tick_values, y_map)
        x_map = lambda x: xy_map(x, tick_values[0])[0]
        polyline_uniform = _polyline_path(xs, uniform_y, x_map, y_map)
        polyline_poly = _polyline_path(xs, poly_y, x_map, y_map)
        parts.append(f'<path class="uniform" d="{polyline_uniform}" />')
        parts.append(f'<path class="poly" d="{polyline_poly}" />')
        parts.append(f'<text class="note" x="{x0:.2f}" y="{y1 + 58:.2f}">{note}</text>')

    parts.append('<text class="note" x="60" y="940">Summary: literal self-gravity removal drives P_c, |E_g|, and the virial support integral to zero; for the n=1 barotrope the finite-mass radius diverges like (1-lambda)^(-1/2).</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def write_scaling_svg(path: str | Path, uniform: UniformSphereToy, poly: N1PolytropeFamily) -> Path:
    output_path = Path(path)
    output_path.write_text(_render_scaling_svg(uniform, poly))
    return output_path


def print_summary() -> None:
    uniform = uniform_sphere_symbolic()
    poly = n1_polytrope_symbolic()

    print("=== Uniform sphere toy ===")
    for key in ["rho", "Pc", "Eg", "Pi", "virial_residual"]:
        print(f"{key} = {sp.sstr(sp.simplify(uniform[key]))}")
    print()

    print("=== n=1 polytrope with fixed M and K ===")
    for key in ["R", "rho_c", "Pc", "Eg", "Pi", "virial_residual"]:
        print(f"{key} = {sp.sstr(sp.simplify(poly[key]))}")
    print()

    print("=== n=1 hydrostatic residual ===")
    print(sp.sstr(n1_hydrostatic_residual_symbolic()))
    print()

    print("=== lambda -> 1 asymptotics (alpha = 1 - lambda -> 0+) ===")
    print("uniform: Pc ~ alpha, |Eg| ~ alpha, Pi ~ alpha")
    print("n=1 fixed M,K: R ~ alpha^(-1/2), rho_c ~ alpha^(3/2), Pc ~ alpha^3, |Eg| ~ alpha^(3/2)")


def main(output_dir: str | Path = ".") -> None:
    output_root = Path(output_dir)
    lam = np.linspace(0.0, 0.999, 500)
    uniform = build_uniform_sphere_family(lam)
    poly = build_n1_polytrope_family(lam)

    table_path = write_scaling_table(output_root / "request2_internal_structure_scalings.tsv", uniform, poly)
    svg_path = write_scaling_svg(output_root / "request2_internal_structure_scalings.svg", uniform, poly)

    print_summary()
    print()
    print(f"wrote {table_path}")
    print(f"wrote {svg_path}")


if __name__ == "__main__":
    main()
