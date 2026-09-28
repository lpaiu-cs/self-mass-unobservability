"""Independent saved-result gates and declared derivative-comparator projection."""
from pathlib import Path
import json
import time
import hashlib
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/direct-eos-gr33/def-orbital-charge-fem'


def projection(rows,degree,amplitudes,phi0):
    # Real coefficient fit, using actual harmonic amplitudes. Equal Euclidean
    # weights here are a declared signal norm, never an observational noise model.
    y=[];design=[]
    static=mp.mpc(*[str(v) for v in rows[0]['charge']])
    for n,row in enumerate(rows[1:],1):
        value=(mp.mpc(*[str(v) for v in row['charge']])-static)*amplitudes[n]/phi0
        basis=[amplitudes[n]/phi0*(-mp.j*n)**j for j in range(degree+1)]
        y.extend([value.real,value.imag]);design.extend([[v.real for v in basis],[v.imag for v in basis]])
    y=mp.matrix(y);A=mp.matrix(design)
    coeff=mp.lu_solve(A,y) if degree==5 else mp.qr_solve(A,y)[0]
    residual=y-A*coeff
    return dict(degree=degree,residual_norm=float(mp.norm(residual)),
        coefficients=[float(v) for v in coeff],residual=[float(v) for v in residual])


def main():
    start=time.monotonic();mp.mp.dps=70
    plan=json.loads((OUT/'plan.json').read_text());result=json.loads((OUT/'result.json').read_text());pilot=json.loads((OUT/'pilot.json').read_text())
    for path,digest in plan['bindings'].items():
        p=Path(path);p=p if p.is_absolute() else ROOT/p
        assert hashlib.sha256(p.read_bytes()).hexdigest()==digest,str(p)
    cases=result['cases'];gates=plan['gates'];assert set(cases)=={'p4','p2','outer','quadrature'}
    maxima={key:max(r[key] for rows in cases.values() for r in rows) for key in ['linear_residual','derivative_reaction_relative','radiative_current_error']}
    assert maxima['linear_residual']<gates['linear_residual']
    assert maxima['derivative_reaction_relative']<gates['endpoint_reaction']
    assert maxima['radiative_current_error']<gates['current']
    for row in result['rows']:
        for name,gate in [('p2','spatial_relative'),('outer','outer_relative'),('quadrature','quadrature_relative')]:
            assert row[name+'_relative']<gates[gate]
            if row['harmonic']:assert row[name+'_static_subtracted_relative']<gates[gate]
    benchmark=json.loads((ROOT/'outputs/direct-eos-gr33/def-resolved-scalar-pulse/companion-benchmark.json').read_text())
    b=benchmark;phi0=mp.mpf(str(b['parameters']['background_phi']));ecc=mp.mpf(str(b['parameters']['eccentricity']))
    scale=-mp.mpf(str(b['leading_drive']['alpha_companion']))*mp.mpf(str(b['inputs']['material_mass_geom_cm']))/(mp.mpf(str(b['parameters']['semimajor_axis_over_star_radius']))*mp.mpf(str(b['inputs']['radius_cm'])))
    amplitudes=[scale]+[2*scale*mp.besselj(n,n*ecc) for n in [1,2,3]]
    # Recompute the readouts from gains and inputs, not from saved amplitudes.
    for n in [1,2,3]:
        response=mp.mpc(*[str(v) for v in cases['p4'][n]['charge']]);static=mp.mpf(str(cases['p4'][0]['charge'][0]))
        direct=float(abs(amplitudes[n]*(response-static)/phi0));saved=result['rows'][n]['frozen_static_subtracted_over_phi0']
        assert abs(direct/saved-1)<1e-10
    projected={name:[projection(rows,d,amplitudes,phi0) for d in [0,1,2,3,4,5]] for name,rows in cases.items()}
    for rows in projected.values():assert rows[-1]['residual_norm']<1e-75
    residual4=projected['p4'][4]['residual'];controls={}
    for name in ['p2','outer','quadrature']:
        controls[name]=float(mp.norm(mp.matrix(projected[name][4]['residual'])-mp.matrix(residual4)))
    signal_sum=sum(r['delta_alpha_over_phi0'] for r in result['rows'][1:])
    difference_sum=sum(r['frozen_static_subtracted_over_phi0'] for r in result['rows'][1:])
    report=dict(classification='Counterexample candidate',all_registered_gates_passed=True,
        static_subtracted_gates_also_passed=True,maxima=maxima,projections=projected,
        projection_definition='Real polynomial in -i*n with freely shared coefficients; actual Kepler amplitudes; Euclidean norm of six real charge coefficients. No data, timing propagation or observational covariance.',
        degree4_projected_control_differences=controls,
        degree4_residual_robust_to_comparisons=projected['p4'][4]['residual_norm']>10*max(controls.values()),
        first_three_harmonics=dict(mechanical_charge_over_phi0_amplitude_sum=signal_sum,
            frozen_static_subtracted_amplitude_sum=difference_sum,
            leading_force_fraction_amplitude_sum=abs(b['leading_drive']['alpha_companion']*float(phi0))*difference_sum,
            scope='Triangle envelope of computed three-harmonic mechanical contrast only; not a certified bound for omitted harmonics, total charge, thermal response, mass normalization or timing residuals.'),
        budget_seconds=pilot['seconds']+result['seconds']+time.monotonic()-start,
        full_objective_complete=False,
        limits=['Mechanical free-minus-supported contrast, not full stellar charge.',
            'Zero frequency is the radiative ADM-fixed limit, not a re-equilibrated thermal sequence.',
            'Leading weak-body companion charge is declared, not microphysically matched.',
            'Background is not thermally stationary; no orbit-long stationary claim.',
            'No rigorous continuous/EOS/omitted-drive error enclosure or observational nuisance applied.'])
    (OUT/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='projections'},indent=2))
    print('DEGREES',[(x['degree'],x['residual_norm']) for x in projected['p4']])


if __name__=='__main__':main()
