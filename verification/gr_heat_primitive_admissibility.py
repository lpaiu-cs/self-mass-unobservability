"""Conserved-data velocity brackets and exact preservation of rest-split inputs."""
from fractions import Fraction as F
import json,sys
import numpy as np
from mpmath import iv
import gr_heat_primitive_global as prior
import gr_outer_product_pilot as outer

g=prior.g;previous=prior.previous;old=prior.old;OUT=g.OUT/'gr-heat-primitive-admissibility';sha=prior.sha;cusp=outer.cusp


def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def exact(value):
    a,b=value.as_integer_ratio();q=F(a,b)
    assert np.longdouble(q.numerator)/np.longdouble(q.denominator)==value
    return q


def prepare():
    previous.verify();prior.verify();assert not OUT.exists();OUT.mkdir();p=json.loads((previous.OUT/'plan.json').read_text())
    files=[g.ROOT/'verification/gr_heat_primitive_admissibility.py',previous.OUT/'plan.json',previous.OUT/'manifest.json',prior.OUT/'manifest.json',
        g.OUT/'initial-state-17-4.npz',prior.heat.OUT/'initial-rates.npz']+sorted(previous.OUT.glob('cell-*-case-*.npz'))
    assert len(list(previous.OUT.glob('cell-*-case-*.npz')))==27
    save('plan.json',dict(classification='Proven',checkpoint='27132476',bits=256,cells=p['cells'],cases=len(p['probes']),velocity_tolerance=p['root_velocity_tolerance'],
        bindings={p.relative_to(g.ROOT).as_posix():sha(p) for p in files},
        theorem='For specified Q, E>0 and any candidate primitive with P>=0, J-Q=v*(E+P) confines v to the interval with endpoints0 and(J-Q)/E. This requires no frozen EOS derivatives. If V=abs(J-Q)/E<1, rho=D*sqrt(1-v^2) lies in[D*sqrt(1-V^2),D]. Nonnegative pressure is an explicit branch-wide premise; it is not inferred from positive pressure at a few points.',
        preserved_input='Read the saved long-double D,E_star,J values by their exact integer ratios. Form E=E_star+D*C in exact rational arithmetic using exactly the prior binary64 C_X*c^2. Never cast the saved conserved arrays to binary64 or directly subtract rest energy.',
        numerical_control='Check all27 manufactured target velocities against the exact bracket; check the recovered approximate velocities using their original declared1e-12 velocity tolerance. Record any exact-bracket overshoot separately, preserving all earlier numerical verdicts.',
        boundary='Kinematic theorem and finite native-target controls. The heat variable must be supplied. No thermal branch existence/connectedness, native EOS interval error, nonuniform cell closure, atmosphere or finite GR trajectory is supplied.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(g.ROOT/rel)==digest,rel
    return p


def run():
    p=bindings();iv.prec=p['bits'];state=np.load(g.OUT/'initial-state-17-4.npz');heats=np.load(prior.heat.OUT/'initial-rates.npz');c=g.c.gr.C*100;rows=[]
    for i in p['cells']:
        C=F(float(state['CX'][i]*c*c));Q=F(float(heats['Q'][i]))
        for case in range(p['cases']):
            data=np.load(previous.OUT/f'cell-{i}-case-{case}.npz');D,Estar,J=map(exact,data['conserved']);E=Estar+D*C
            assert D>0 and E>0;end=(J-Q)/E;lo,hi=min(F(0),end),max(F(0),end);V=abs(end);assert V<1
            truth=F(float(data['true_primitive'][2]));recovered=F(float(data['recovered_primitive'][2]));overshoot=max(F(0),lo-recovered,recovered-hi)
            density_lower=cusp.I(D)*iv.sqrt(1-cusp.I(V)**2);fractional_width=1-iv.sqrt(1-cusp.I(V)**2)
            rows.append(dict(cell=i,case=case,exact_normal_density=str(D),exact_rest_subtracted_energy=str(Estar),exact_total_energy=str(E),exact_momentum=str(J),exact_heat=str(Q),
                velocity_bracket=list(map(str,[lo,hi])),maximum_absolute_velocity=str(V),manufactured_velocity=str(truth),recovered_velocity=str(recovered),
                manufactured_inside_exact_bracket=lo<=truth<=hi,recovered_exact_bracket_overshoot=str(overshoot),
                recovered_inside_original_tolerance=overshoot<=F(str(p['velocity_tolerance'])),
                density_lower_enclosure=cusp.interval_text(density_lower),relative_density_width_upper=str(cusp.high(fractional_width))))
    assert len(rows)==27
    save('result.json',dict(classification='Counterexample candidate',completed=True,passed=all(r['manufactured_inside_exact_bracket'] and r['recovered_inside_original_tolerance'] for r in rows),
        rows=rows,maximum_absolute_velocity=str(max(F(r['maximum_absolute_velocity']) for r in rows)),
        maximum_relative_density_width=str(max(F(r['relative_density_width_upper']) for r in rows)),
        largest_recovered_exact_bracket_overshoot=str(max(F(r['recovered_exact_bracket_overshoot']) for r in rows)),
        saved_longdouble_significand_bits=int(np.finfo(np.longdouble).nmant+1),saved_longdouble_storage_bits=int(np.dtype(np.longdouble).itemsize*8),
        actual_EOS_branch_certified=False,global_native_inverse_certified=False,finite_GR_evolution=False))
    save('manifest.json',dict(sha256={f.relative_to(g.ROOT).as_posix():sha(f) for f in OUT.iterdir() if f.is_file()}));verify()


def verify():
    p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['completed'] and len(r['rows'])==27
    for row in r['rows']:
        E,J,Q=map(F,[row['exact_total_energy'],row['exact_momentum'],row['exact_heat']]);end=(J-Q)/E
        assert row['velocity_bracket']==list(map(str,[min(F(0),end),max(F(0),end)])) and abs(end)==F(row['maximum_absolute_velocity'])<1
        lo,hi=map(F,row['velocity_bracket']);v=F(row['recovered_velocity']);over=max(F(0),lo-v,v-hi)
        assert over==F(row['recovered_exact_bracket_overshoot']) and (over<=F(str(p['velocity_tolerance'])))==row['recovered_inside_original_tolerance']
    assert r['passed']==all(x['manufactured_inside_exact_bracket'] and x['recovered_inside_original_tolerance'] for x in r['rows'])
    print('PASS exact conserved-data bracket audit; finite native-target membership:',r['passed'],flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
