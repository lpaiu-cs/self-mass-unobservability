"""Proven: algebra checks. Imported: stored-table arithmetic; no new timing fit."""
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "symbolic"))
sys.path.insert(0, str(ROOT / "request10_external/scripts"))
from nonanalytic_jet_demo import nonanalytic_jet_summary
from sep_common import u95_of


def main():
    cases = {c.case_id: c for c in nonanalytic_jet_summary().cases}
    assert cases["smooth_flat_single_coordinate"].finite_taylor_jet_valid
    assert not cases["threshold_sqrt_activation"].finite_taylor_jet_valid
    y = sp.Symbol("y", positive=True)
    flat = sp.exp(-1 / y**2)
    for n in range(7):
        assert sp.limit(sp.diff(flat, y, n), y, 0, dir="+") == 0
        assert sp.limit(flat / y**n, y, 0, dir="+") == 0
    assert sp.limit(sp.sqrt(y) / y**5, y, 0, dir="+") == sp.oo

    z = sp.Symbol("z")
    tau = sp.Rational(2, 3)
    beta, cy = sp.Rational(-3, 2), sp.Rational(7, 4)
    transfer = cy + beta / (1 + tau * z)
    for k in range(1, 5):
        q = sp.prod(z*z + w*w for w in range(1, k+1))
        p = sp.cancel(cy + beta * (1 - q / q.subs(z, -1/tau)) / (1 + tau*z))
        assert sp.degree(p, z) == 2*k-1
        assert all(c.is_real for c in sp.Poly(p, z).all_coeffs())
        assert all(sp.simplify((p-transfer).subs(z, sign*sp.I*w)) == 0
                   for w in range(1, k+1) for sign in (-1, 1))
        coeffs = sp.symbols(f"a0:{2*k-1}", real=True)
        too_short = sum(a*z**n for n, a in enumerate(coeffs))
        residuals = [(too_short-transfer).subs(z, sp.I*w).expand(complex=True)
                     for w in range(1, k+1)]
        equations = [part for r in residuals for part in (sp.re(r), sp.im(r))]
        assert sp.linsolve(equations, coeffs) == sp.EmptySet
    p5 = -z**5/100 + z**4/100 - 3*z**3/20 + 3*z*z/20 - 16*z/25 + sp.Rational(16,25)
    assert all(sp.simplify((p5-1/(1+z)).subs(z, sign*sp.I*w)) == 0
               for w in (1,2,3) for sign in (-1,1))
    for n in range(5):
        p = cy + beta*sum((-tau*z)**j for j in range(n+1))
        assert sp.cancel(transfer-p-beta*(-tau*z)**(n+1)/(1+tau*z)) == 0
    t, w, relaxation = sp.symbols("t w relaxation", positive=True)
    alpha, f0 = sp.symbols("alpha f0", real=True)
    chi = alpha*f0*(sp.cos(w*t)+w*relaxation*sp.sin(w*t))/(1+w*w*relaxation**2)
    assert sp.simplify(relaxation*sp.diff(chi,t)+chi-alpha*f0*sp.cos(w*t)) == 0

    x, yy, zz = sp.symbols("x yy zz", real=True)
    r, gm = sp.symbols("r gm", positive=True)
    potential = -gm/sp.sqrt(x*x+yy*yy+zz*zz)
    tidal = sp.hessian(potential, (x, yy, zz)).subs({x:r, yy:0, zz:0})
    i2 = sp.simplify(sp.trace(tidal*tidal))
    assert i2 == 6*gm**2/r**6
    assert sp.diff(i2,r) == -36*gm**2/r**7
    assert sp.diff(-gm/r,r) == gm/r**2
    signal = sp.Matrix([0,1])
    assert (signal.T*(sp.eye(2)-sp.eye(2))*signal)[0] == 0
    assert (signal.T*(sp.eye(2)-sp.diag(1,0))*signal)[0] == 1

    folder = ROOT / "request10_external/sep_dynamic"
    data = json.loads((folder / "sep_phase_marg_10_8e.json").read_text())
    gate = json.loads((folder / "sep_gateG2wp.json").read_text())["anchor"]
    value = u95_of(gate["beta_hat_lin"], data["K_dyn"]*gate["sigma_fisher"])
    assert math.isclose(value, data["anchors"]["tau_2"]["u95pm_K10"], rel_tol=1e-12)
    manuscript = (ROOT / "paper/manuscript.md").read_text(encoding="utf-8")
    for key, row in data["anchors"].items():
        lag = key.removeprefix("tau_")
        table = next(line for line in manuscript.splitlines() if line.startswith(lag+" & "))
        numbers = re.findall(r"(\d+\.\d+)\\times10\^\{(-?\d+)\}", table)
        expected = [row[k] for k in ("u95pm_K10_fullrank", "u95pm_K10", "u95pm_fisher")]
        assert len(numbers) == 3
        for (mantissa, exponent), stored in zip(numbers, expected):
            assert math.isclose(float(mantissa)*10**int(exponent), stored, rel_tol=5e-4)
        ratio = float(table.split(" & ")[-1].rstrip("\\"))
        assert abs(ratio-expected[0]/expected[1]) < 0.0051
    assert data["detection_candidate"] is False
    # Verify the new descriptive table directly against the frozen per-condition counts.
    folder=ROOT/'outputs/research-completion'
    coverage=json.loads((folder/'coverage-audit.json').read_text())
    estimated=json.loads((folder/'estimated-covariance-audit.json').read_text())
    conditions=[('White noise','white',0.),('Omitted mean, norm 3','omitted_R3',None),
                ('Omitted mean, norm 10','omitted_R10',None),('Omitted mean, norm 30','omitted_R30',None),
                ('Extra Fourier RMS 0.25','extra_fourier_025',.25),('Extra Fourier RMS 1','extra_fourier_1',1.)]
    for label,scenario,amplitude in conditions:
        row=next(line for line in manuscript.splitlines() if line.startswith(label+' & '))
        values=row.rstrip('\\').split(' & ')[1:]
        for value,method in zip(values[:3],['truncated_diag','full_diag','full_matched_GLS']):
            measured=min(r['U_K1']['fraction'] for r in coverage['rows'] if r['scenario']==scenario and r['method']==method)
            assert abs(float(value)-measured)<5.1e-5,(label,method,measured)
        if amplitude is not None:
            measured=min(r['U_K1']['fraction'] for r in estimated['rows'] if r['true_a']==amplitude)
            assert abs(float(values[3])-measured)<5.1e-5
    for results in [coverage,estimated]:
        for row in results['rows']:
            cell=row['U_K1']; assert cell['n']==8192
            assert cell['fraction']==cell['hits']/cell['n']
            assert cell['lo95']<=cell['fraction']<=cell['hi95']
    assert coverage['gaussian_control_screen_pass'] and estimated['screen_pass']
    assert max(r['upper_fraction'] for r in estimated['amplitude_estimates'])==0
    matching=json.loads((folder/'physical-matching.json').read_text())
    assert len(matching['checks'])==14 and matching['common_origin_cannot_fix_closure']
    assert abs(matching['physical_closure_radians']-3.11837236)<1e-8
    comparison=json.loads((folder/'comparator-audit.json').read_text())
    assert len(comparison['rows'])==378 and len(comparison['positive_fast_spectrum_witnesses'])==54
    for row in comparison['rows']:
        if row['comparator']=='P5':
            assert row['exact_collapse'] and row['relative_information']<1e-18
        else:
            assert row['continuous_phase_information_lower_bound']>0 and row['unit_sigma_beta']>0
    for row in comparison['positive_fast_spectrum_witnesses']:
        assert row['slow_positive_beta_witness']>0 and row['minimum_distance_per_positive_beta']>0
    assert comparison['gap_is_physically_matched'] is False
    phases=json.loads((folder/'phase-state-audit.json').read_text())
    refined=json.loads((folder/'phase-refinement.json').read_text())
    assert len(phases['phase_rows'])==54 and len(refined['rows'])==18
    for row in refined['rows']:
        seeds=[r['grid_max_U'] for r in phases['phase_rows'] if r['a']==row['a'] and r['tau']==row['tau']]
        assert max(seeds)<=row['refined_U']*(1+1e-12)<=row['continuous_upper_bound']
        assert not row['global_maximum_certified']
        assert all(not r['move_cap_levels'] for r in row['runs'])
    assert max(r['refined_over_legacy'] for r in refined['rows'])<1.03722
    gap=phases['gaps']; assert gap['gaps']==565 and len(gap['rows'])==1130
    for a in [0.,1.]:
        rows=[r for r in gap['rows'] if r['a']==a]
        assert len({r['after_sorted_index'] for r in rows})==565
        assert all(r['gap_days']>1 for r in rows)
    weakest=min((r for r in gap['rows'] if r['a']==1),key=lambda r:min(r['delta_chi2_plus'],r['delta_chi2_minus']))
    assert math.isclose(weakest['gap_days'],223.3708406076462,rel_tol=1e-10)
    assert min(weakest['delta_chi2_plus'],weakest['delta_chi2_minus'])<.153
    assert all(r['transient_outside_known_drive_span']>.79 for r in phases['transients']['rows'])
    assert phases['numerical']['full_derivative_error_certificate'] is False
    import numpy as np
    physical=json.loads((folder/'corrected-physical-drive.json').read_text())
    assert math.isclose(sum(physical['drive']['normalized_amplitudes']),1.,abs_tol=1e-12)
    assert physical['drive']['mass_o_parameter_difference']<1e-12
    assert len(physical['fits'])==6 and physical['omitted_timing_error_bound'] is None
    calibrated=json.loads((folder/'simultaneous-calibration.json').read_text())
    validated=json.loads((folder/'simultaneous-validation.json').read_text())
    assert calibrated['seed']!=validated['seed'] and calibrated['threshold']==validated['threshold']
    assert validated['calibration_sha256']==hashlib.sha256((folder/'simultaneous-calibration.json').read_bytes()).hexdigest()
    for row in validated['rows']:
        assert row['inclusion']['hits']+row['null_false_positive']['hits']==8192
    for row in validated['data']['physical_lag_sections']:
        assert not row['empty'] and row['joint_region_abs_beta_upper']>=2*abs(row['beta'])
    live=json.loads((folder/'runtime12-analysis.json').read_text())
    assert len(live['transient'])==3 and all(r['convergence_5percent_pass'] for r in live['transient'])
    assert len(live['derivative_rows'])==28 and live['halfstep_basis']['rigorous_derivative_error_bound'] is None
    for path in (folder/'runtime12').glob('jac_*.npz'):
        z=np.load(path)
        assert np.allclose(z['dcol'],(z['plus']-z['minus'])/(2*z['h']))
    assert len(list((folder/'runtime12').glob('jac_*.npz')))==35
    assert all(r['eccentricity']>1 for r in live['original_displacement_rejection'])
    assert max(r['fraction'] for r in live['nonlinear_gap'])==.01
    assert len(json.loads((folder/'state-identifiability.json').read_text())['checks'])==14
    pairs=json.loads((folder/'gap-pair-audit.json').read_text())
    assert all(len(r['candidates'])==201 and any(r['coefficients']) for r in pairs['rows'])
    manifest = ROOT / "paper/revision-manifest.json"
    if manifest.exists():
        for path, digest in json.loads(manifest.read_text())["sha256"].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest, path
    print("PASS: analytic boundaries; historical tables; corrected drive; independent joint-region validation; live response and failed promotion gates; state identities; manifest.")


if __name__ == "__main__":
    main()
