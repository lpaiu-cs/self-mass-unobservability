"""Independent IVP control that visits every stored coefficient interval."""
import json, sys
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
import gr_scalar_residual_runner as proof

g=proof.g;regular=proof.original.regular;OUT=g.OUT/'gr-scalar-piecewise-control'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();proof.verify()
    paths=[g.ROOT/'verification/gr_scalar_piecewise_control.py',proof.OUT/'manifest.json',
        regular.OUT/'manifest.json',regular.original.OUT/'field-2e-12.npz']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='427ebf8',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        cases=[dict(tolerance=2e-12,centre_seed=1e-9),dict(tolerance=2e-13,centre_seed=1e-12)],
        method='Reuse DOP853 and the stored cubics but restart at every coefficient boundary, passing the previous field/flux onward. The first case retains the original tolerance and seed. The second lowers both. The same declared exterior normalization and R/M are used.',
        comparison='Check the computed coefficient against the already frozen rigorous residual enclosure; retain noninclusion. This is an independent finite numerical control, not a replacement for the certificate or a new physical model.',
        boundary='Seed approximation, floating coefficient evaluation and IVP error are not certified by these two runs. The independent exact-polynomial certificate includes the full centre interval and remains authoritative for the declared coefficients.'))


def run():
    mp.mp.dps=80;mp.iv.dps=80;plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    certified=json.loads((proof.OUT/'result.json').read_text());r=certified['exact_rationals']
    p={k:regular.Q(v) for k,v in json.loads((regular.OUT/'result.json').read_text())['exact_rationals'].items()}
    L=-mp.iv.ln(1-2*regular.iv(p['mu']))/(2*regular.iv(p['mu']))
    enclosure,_,_=proof.original.response(regular.Q(r['p_surface']),regular.Q(r['J_surface']),regular.Q(r['residual_bound']),
        p['kappa'],p['B2'],L,p['R_over_M'])
    data=dict(np.load(regular.original.OUT/'field-2e-12.npz'));x=data['x'];rows=[]
    for case in plan['cases']:
        tol=case['tolerance'];seed=case['centre_seed'];assert 0<seed<x[1]
        a0=data['p_coefficients'][-1,0];b0=data['b_coefficients'][-1,0]
        y=np.array([1+b0/a0*seed**2/6,b0*seed**3/3]);values=[y.copy()];calls=0
        for i,(left,right) in enumerate(zip(x,x[1:])):
            ac=data['p_coefficients'][:,i];bc=data['b_coefficients'][:,i]
            def rhs(t,y):
                z=t-left;a=((ac[0]*z+ac[1])*z+ac[2])*z+ac[3]
                b=((bc[0]*z+bc[1])*z+bc[2])*z+bc[3]
                return [y[1]/(t*t*a),t*t*b*y[0]]
            sol=solve_ivp(rhs,(seed if i==0 else left,right),y,method='DOP853',
                rtol=tol,atol=[tol*1e-2,tol*1e-6],max_step=right-left)
            assert sol.success;calls+=sol.nfev;y=sol.y[:,-1];values.append(y.copy())
        mu=float(p['mu']);normal=y[0]-y[1]*np.log1p(-2*mu)/(2*mu)
        value=float(p['R_over_M'])*y[1]/normal;v=regular.iv(regular.rational(value))
        inside=bool(enclosure.a<=v.a and v.b<=enclosure.b)
        rows.append(dict(**case,coefficient=value,inside_prior_certificate=inside,RHS_calls=calls,
            minimum_sampled_unit_central_field=float(np.min(np.array(values)[:,0]))))
        np.savez_compressed(OUT/f'field-{tol}.npz',x=np.r_[seed,x[1:]],values=np.array(values))
        print('PIECEWISE SCALAR',rows[-1],flush=True)
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        all_inside=all(r['inside_prior_certificate'] for r in rows),
        finite_difference=abs(rows[1]['coefficient']-rows[0]['coefficient']),
        old_global_IVP_coefficient=certified['prior_numerical_coefficient'],
        physical_EOS_certified=False,continuous_IVP_error_certified_by_this_control=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    proof.verify()
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS all-piece IVP control bindings; inspect certificate inclusion',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
