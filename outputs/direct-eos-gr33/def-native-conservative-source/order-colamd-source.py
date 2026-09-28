"""One budgeted consistent-mass order control, preserving the fixed C1 source."""
import argparse
import json
import signal
import resource
import time
import inspect
import textwrap
from pathlib import Path
import numpy as np
import def_native_conservative_source as task
from verify_native_conservative_source import compare

OUT=task.OUT
source=task.prior.source
anchor='cuts=np.unique(np.r_[self.cells,bg.x,self.edges,self.native]);gx,gw=np.polynomial.legendre.leggauss(6)'
assert source.count(anchor)==1
source=source.replace(anchor,anchor.replace('leggauss(6)','leggauss(8)'))
anchor='vx=(self.edges[:-1,None]+np.diff(self.edges)[:,None]*(gx+1)/2)'
assert source.count(anchor)==1
source=source.replace(anchor,'gx,gw=np.polynomial.legendre.leggauss(6)\n    '+anchor)
namespace=dict(task.assemble.__globals__);exec(compile(source,__file__,'exec'),namespace)
assemble=namespace['__init__']
stage_source=textwrap.dedent(inspect.getsource(task.old.Model.stage))
assert stage_source.count("permc_spec='NATURAL'")==2
stage_source=stage_source.replace("permc_spec='NATURAL'","permc_spec='COLAMD'",1)
stage_namespace=dict(task.old.Model.stage.__globals__)
exec(compile(stage_source,__file__,'exec'),stage_namespace)


class Model(task.Model):
    stage=stage_namespace['stage']
    def __init__(self):
        previous=task.assemble;task.assemble=assemble
        try:super().__init__(6,lumped=False)
        finally:task.assemble=previous
        # These quadrature evaluation matrices are assembly-only caches. The
        # timestep uses K/M/F, nativeV/nativeD and the compiled photon_test.
        # Releasing them changes no operator and keeps the declared5GB cap.
        del self.V,self.D


def plan():
    assert not (OUT/'order-plan.json').exists()
    control=json.loads((OUT/'control-result.json').read_text())
    assert not control['mass_rule_passed']
    task.write(OUT/'order-plan.json',dict(classification='Counterexample candidate',
        claim='Resolve whether the consistent-mass p4 response agrees with p6 for the identical conservative source, after the GLL/consistent-mass contrast failed.',
        rationale='Paired GLL p2/p4 agreement was insufficient. Use the consistent positive variational mass as reference; preserve the failed GLL control. This is one explicit budget reassessment, not a sequence of local patches or retuned sources.',
        changed='Degree6 and positive8-point mechanical quadrature. The p6 mass needs degree12 products, so the original6-point rule is not sufficient. Retain original6-point heat volumes and all source/profile/background/horizon data.',
        paths=['consistent-p6-8 pilot','consistent-p6-32','consistent-p6-64'],
        gates=dict(temperature=.02,velocity=.03,scalar=.03),
        comparisons='Time: consistent p6 32/64. Space: saved consistent p4-64 versus consistent p6-64. Source and mass-rule failures stay separate.',
        prior_accounted_seconds=control['accounted_total_seconds']+5,
        prior_audit_allowance_seconds=5,total_budget_seconds=240,
        stop='Use measured pilot to decide within remaining budget. No p8, more cells, changed source slopes or longer horizon follows this test.',
        scientific_boundary='Convergence of this fixed source approximation cannot certify the actual physical subcell heat profile or moving atmosphere.',
        bindings={str(p.relative_to(task.old.ROOT)):task.old.photons.digest(p) for p in [Path(__file__),Path(task.__file__),OUT/'control-result.json',OUT/'consistent-p4-64.npz']}))
    (OUT/'order-assembler.py').write_bytes(source.encode())


def run():
    assert not (OUT/'order-result.json').exists();spec=json.loads((OUT/'order-permuted-plan.json').read_text())
    for name,h in spec['bindings'].items():assert task.old.photons.digest(task.old.ROOT/name)==h,name
    previous=spec['prior_accounted_seconds'];signal.alarm(int(240-previous));start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(int(5e9),int(5e9)))
    horizon=json.loads((task.old.OUT/'plan.json').read_text())['horizon_seconds'];m=Model()
    pilot=m.evolve(horizon,8,'consistent-p6-pilot',True)
    forecast=previous+time.monotonic()-start+1.25*(pilot['seconds']*96/8+5)
    task.write(OUT/'order-pilot-budget.json',dict(classification='Counterexample candidate',forecast_total_seconds=forecast,
        setup_seconds=m.setup_seconds,pilot_evolution_seconds=pilot['seconds'],dofs=m.size,
        assumption='Actual p6 setup/pilot, scale96 remaining steps with25percent margin plus5seconds. No path beyond these two is allowed.'))
    assert forecast<240,'Stop: measured p6 controls do not fit remaining total budget'
    for n in [32,64]:m.evolve(horizon,n,f'consistent-p6-{n}',True)
    a=np.load(OUT/'consistent-p6-32.npz');b=np.load(OUT/'consistent-p6-64.npz');c=np.load(OUT/'consistent-p4-64.npz')
    time_error=compare(a,{f:b[f][::2] for f in ['temperature','velocity','scalar']})
    space_error=compare(c,b);passed=all(time_error[f]<g and space_error[f]<g for f,g in spec['gates'].items())
    data=dict(classification='Counterexample candidate',passed=bool(passed),time_relative=time_error,space_relative=space_error,
        seconds=time.monotonic()-start,accounted_total_seconds=previous+time.monotonic()-start,
        original_GLL_mass_control_passed=False,physical_source_profile_certified=False,moving_surface_solved=False,
        final_charge_solved=False,full_goal_complete=False)
    task.write(OUT/'order-result.json',data);signal.alarm(0);print('CONSISTENT ORDER',json.dumps(data),flush=True)


def correct():
    assert not (OUT/'order-corrected-plan.json').exists()
    spec=json.loads((OUT/'order-plan.json').read_text());key=str(Path(__file__).relative_to(task.old.ROOT))
    assert task.old.photons.digest(OUT/'order-registered-source.py')==spec['bindings'][key]
    spec['bindings'][key]=task.old.photons.digest(Path(__file__))
    spec['prior_accounted_seconds']+=10
    spec['resource_failure_allowance_seconds']=10
    spec['repair']='Release the assembly-only V/D quadrature caches before evolution. The initial p6 run reached the5GB virtual-address cap during first-stage scaling before completing the pilot. No numerical gate, source, degree, horizon or memory budget is changed.'
    task.write(OUT/'order-resource-failure.json',dict(classification='Counterexample candidate',passed=False,
        failure='ArrayMemoryError: allocation of11.6MiB failed under the5GB RLIMIT_AS during first-stage sparse column scaling.',
        evolution_completed=False,source=spec['repair']))
    task.write(OUT/'order-corrected-plan.json',spec)


def reorder():
    assert not (OUT/'order-permuted-plan.json').exists()
    spec=json.loads((OUT/'order-corrected-plan.json').read_text());key=str(Path(__file__).relative_to(task.old.ROOT))
    assert task.old.photons.digest(OUT/'order-cache-source.py')==spec['bindings'][key]
    spec['bindings'][key]=task.old.photons.digest(Path(__file__))
    spec['prior_accounted_seconds']+=10
    spec['second_resource_failure_and_diagnosis_allowance_seconds']=10
    spec['linear_ordering_repair']='Use SuperLU COLAMD for the coupled block. Keep the separate GR solve unchanged. The cache-only repair still failed during the first LU factorization at5GB; measured post-assembly VmRSS was0.452GB and VmSize0.699GB. This targets LU fill, not physical memory-budget expansion. All equation residual checks remain.'
    task.write(OUT/'order-second-resource-failure.json',dict(classification='Counterexample candidate',passed=False,
        failure='SuperLU reported not enough memory under unchanged5GB address-space cap before completing the first pilot stage.',
        evolution_completed=False,action=spec['linear_ordering_repair']))
    task.write(OUT/'order-permuted-plan.json',spec)
    (OUT/'order-stage.py').write_bytes(stage_source.encode())


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['plan','correct','reorder','run']);globals()[parser.parse_args().action]()
