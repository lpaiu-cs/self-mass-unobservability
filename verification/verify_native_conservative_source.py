"""Independent saved-path replay, source conservation and controlled ablation."""
import argparse
import json
import resource
import signal
import time
import numpy as np
from pathlib import Path
import def_native_conservative_source as task

OUT=task.OUT;old=task.old


class GammaOnly(task.prior.Model):
    def __init__(self,degree):
        previous=task.prior.assemble;task.prior.assemble=task.assemble
        try:super().__init__(degree)
        finally:task.prior.assemble=previous
    evolve=task.Model.evolve


def plan():
    assert not (OUT/'control-plan.json').exists()
    result=json.loads((OUT/'result.json').read_text());assert result['passed']
    p4=json.loads((OUT/'p4-64.json').read_text());p2=json.loads((OUT/'p2-64.json').read_text())
    forecast=result['accounted_seconds']+1.5*(2*p4['seconds']+p2['seconds']+3*p4['setup_seconds']+12)
    assert forecast<240
    task.write(OUT/'control-plan.json',dict(classification='Counterexample candidate',
        decision='After the fixed candidate passed, test mass-rule dependence and isolate source reconstruction from Gamma1 correction. No mesh, source slopes, horizon or acceptance thresholds are tuned.',
        paths=['consistent mass C1 p4-64','Gamma1-only piecewise-constant source p4-64','Gamma1-only piecewise-constant source p2-64'],
        mass_rule_gates=dict(temperature=.02,velocity=.03,scalar=.03),
        ablation='Report both Gamma1-only paths with unchanged max-norm spatial gates. They are a causal contrast, not another physical result accepted without time controls.',
        conservation_gate=1e-9,total_budget_seconds=240,forecast_total_seconds=forecast,
        prior_accounted_seconds=result['accounted_seconds'],
        forecast_assumption='Actual p4/p2 costs plus three assemblies and12seconds, with50percent margin; the consistent-mass solve cost is still estimated.',
        bindings={str(p.relative_to(old.ROOT)):old.photons.digest(p) for p in [Path(__file__),Path(task.__file__),OUT/'result.json',OUT/'corrected-plan.json']}))


def compare(a,b):
    return {field:float(np.max(abs(a[field]-b[field]))/max(np.max(abs(b[field])),1e-100)) for field in ['temperature','velocity','scalar']}


def run():
    assert not (OUT/'control-result.json').exists()
    spec=json.loads((OUT/'control-plan.json').read_text())
    for name,h in spec['bindings'].items():assert old.photons.digest(old.ROOT/name)==h,name
    start=time.monotonic();signal.alarm(int(240-spec['prior_accounted_seconds']))
    resource.setrlimit(resource.RLIMIT_AS,(int(5e9),int(5e9)))
    horizon=json.loads((old.OUT/'plan.json').read_text())['horizon_seconds']
    m=task.Model(4,lumped=False);m.evolve(horizon,64,'consistent-p4-64',True)
    base=np.load(OUT/'p4-64.npz');other=np.load(OUT/'consistent-p4-64.npz')
    errors=compare(other,base);mass_pass=all(e<spec['mass_rule_gates'][f] for f,e in errors.items())
    # Audit the actual curved geometry by independent 8-point quadrature on
    # each original heat cell. This is not the endpoint identity reused.
    x,w=np.polynomial.legendre.leggauss(8)
    widths=np.diff(m.edges);r=(m.edges[:-1,None]+widths[:,None]*(x+1)/2).ravel();p=m.bg.sample(r)
    _,loss,J=m.sources(p);geo=old.G*m.bg.R**2/old.C**4
    volume_weight=4*np.pi*m.bg.R**3*p['N']/np.sqrt(1-2*p['m']/np.maximum(r,1e-100))*r*r
    local=(loss@base['E'])/geo*np.exp(-8*p['phi']**2)*volume_weight
    integral=widths/2*np.sum(local.reshape(m.n,8)*w,axis=1)
    expected=np.diff(np.r_[np.longdouble(0),base['E']])
    defect=float(max(abs(integral-expected))/max(max(abs(expected)),1e-100))
    assert defect<spec['conservation_gate'],defect
    # A common cumulative map supplies both the mass debit and heat lift.
    faces=m.face_map(m.edges)@base['E'];face_error=float(max(abs(faces-np.r_[0.,base['E']]))/max(abs(base['E'])))
    assert face_error<1e-12
    field=m.bg.sample(m.native);native_gamma=float(max(abs(field['gamma']/m.raw[:,4]-1)));assert native_gamma==0
    del m
    ablation=None
    if mass_pass:
        for degree in [4,2]:
            g=GammaOnly(degree);g.evolve(horizon,64,f'gamma-only-p{degree}-64',True);del g
        ablation=compare(np.load(OUT/'gamma-only-p2-64.npz'),np.load(OUT/'gamma-only-p4-64.npz'))
    data=dict(classification='Counterexample candidate',mass_rule_passed=mass_pass,mass_rule_relative=errors,
        gamma_only_spatial_relative=ablation,independent_curved_heat_integral_relative=defect,
        shared_face_inventory_relative=face_error,current_native_gamma_relative=native_gamma,
        seconds=time.monotonic()-start,accounted_total_seconds=spec['prior_accounted_seconds']+time.monotonic()-start,
        physical_source_profile_certified=False,moving_surface_solved=False,final_charge_solved=False,full_goal_complete=False)
    task.write(OUT/'control-result.json',data);signal.alarm(0);print('SOURCE CONTROLS',json.dumps(data),flush=True)


def audit():
    result=json.loads((OUT/'result.json').read_text());checks=json.loads((OUT/'control-result.json').read_text())
    source=task.controls();a=np.load(OUT/'p4-32.npz');b=np.load(OUT/'p4-64.npz');c=np.load(OUT/'p2-64.npz')
    spatial=compare(c,b);worst={}
    for field in spatial:
        t=float(np.max(abs(a[field]-b[field][::2]))/max(np.max(abs(b[field])),1e-100))
        assert abs(t-result['comparisons'][field]['time_relative'])<1e-14
        assert abs(spatial[field]-result['comparisons'][field]['space_relative'])<1e-14
        difference=abs(c[field]-b[field]);j,i=np.unravel_index(np.argmax(difference),difference.shape)
        worst[field]=dict(native_array_index=int(i),radius_fraction=float(b['radius'][i]),p4=float(b[field][j,i]),p2=float(c[field][j,i]))
    for label in ['pilot','p4-32','p4-64','p2-64','consistent-p4-64','gamma-only-p4-64','gamma-only-p2-64']:
        p=OUT/(label+'.npz')
        if not p.exists():continue
        d=np.load(p);report=json.loads((OUT/(label+'.json')).read_text())
        assert all(np.all(np.isfinite(d[k])) for k in d.files)
        assert all(np.max(abs(d[k][0]))==0 for k in spatial)
        assert report['max_linear_residual']<1e-9 and report['max_heat_identity']<1e-9
    data=dict(classification='Counterexample candidate',artifact_checks_passed=True,symbolic=source,worst=worst,
        registered_response_passed=result['passed'],mass_rule_passed=checks['mass_rule_passed'],
        original_piecewise_constant_failure_preserved=True,physical_interpolation_independence=False,
        moving_surface_solved=False,final_charge_solved=False,full_goal_complete=False)
    task.write(OUT/'audit.json',data);print(json.dumps(data),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['plan','run','audit']);globals()[parser.parse_args().action]()
