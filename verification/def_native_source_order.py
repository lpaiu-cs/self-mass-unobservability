"""One budgeted consistent-mass order control, preserving the fixed C1 source."""
import argparse
import json
import signal
import resource
import time
import inspect
import textwrap
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.linalg.lapack import dgbtrf,dgbtrs
from scipy.sparse import coo_matrix
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


def band_factor(matrix):
    # The measured radial permutation has half-bandwidth50/52. LAPACK bounds
    # pivot fill to the band and avoids the multi-GB sparse-LU allocation.
    data=matrix.tocoo();lower=int(max(data.row-data.col));upper=int(max(data.col-data.row))
    assert lower<=64 and upper<=64,(lower,upper)
    ab=np.zeros((2*lower+upper+1,matrix.shape[0]),order='F')
    ab[lower+upper+data.row-data.col,data.col]=data.data
    lu,piv,info=dgbtrf(ab,lower,upper,overwrite_ab=True);assert info==0,info
    def solve(rhs):
        answer,info=dgbtrs(lu,lower,upper,np.asarray(rhs,float),piv)
        assert info==0,info
        return answer
    return SimpleNamespace(solve=solve)


def condense_factor(matrix,model):
    # C0 element bubbles couple to other bubbles only in their own element.
    # Eliminate those small blocks exactly, keeping every endpoint and heat
    # variable. This changes the linear solve, never the approximation space.
    group=np.full(model.size,-1,int);slot=np.full(model.size,-1,int)
    nodes=np.flatnonzero(np.arange(len(model.grid))%model.degree)
    for field in range(2):
        indices=model.indices[nodes,field];valid=indices>=0
        group[indices[valid]]=nodes[valid]//model.degree
        slot[indices[valid]]=2*(nodes[valid]%model.degree-1)+field
    grouped=np.r_[group,np.full(model.n,-1)][model.permutation]
    slotted=np.r_[slot,np.full(model.n,-1)][model.permutation]
    bubble=np.flatnonzero(grouped>=0);keep=np.flatnonzero(grouped<0)
    groups=grouped[bubble];slots=slotted[bubble];width=2*(model.degree-1);count=len(model.cells)-1
    bb=matrix[bubble,:][:,bubble].tocoo()
    assert np.all(groups[bb.row]==groups[bb.col]),'Bubble block must be element-local'
    blocks=np.broadcast_to(np.eye(width),(count,width,width)).copy()
    blocks[groups,slots,slots]=0
    blocks[groups[bb.row],slots[bb.row],slots[bb.col]]=bb.data
    inverse=np.linalg.inv(blocks)
    lookup=np.full((count,width),-1,int);lookup[groups,slots]=np.arange(len(bubble))
    columns=lookup[groups];rows=np.broadcast_to(np.arange(len(bubble))[:,None],columns.shape)
    values=inverse[groups[:,None],slots[:,None],np.arange(width)[None,:]];valid=columns>=0
    inv=coo_matrix((values[valid],(rows[valid],columns[valid])),shape=(len(bubble),len(bubble))).tocsc()
    inv.eliminate_zeros()
    kb=matrix[keep,:][:,bubble].tocsc();lift=(inv@matrix[bubble,:][:,keep]).tocsc()
    schur=(matrix[keep,:][:,keep]-kb@lift).tocsc();schur.eliminate_zeros()
    reduced=band_factor(schur)
    def solve(rhs):
        z=inv@rhs[bubble];x=reduced.solve(rhs[keep]-kb@z)
        answer=np.empty(len(rhs));answer[keep]=x;answer[bubble]=z-lift@x
        return answer
    return SimpleNamespace(solve=solve)


stage_source=textwrap.dedent(inspect.getsource(task.old.Model.stage))
assert stage_source.count("permc_spec='NATURAL'")==2
stage_source=stage_source.replace("factor=splu((aa@diags(1/col)).tocsc(),permc_spec='NATURAL')","factor=condense_factor((aa@diags(1/col)).tocsc(),self)")
stage_namespace=dict(task.old.Model.stage.__globals__)
stage_namespace['band_factor']=band_factor
stage_namespace['condense_factor']=condense_factor
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
    assert not (OUT/'order-result.json').exists();spec=json.loads((OUT/'order-condensed-plan.json').read_text())
    for name,h in spec['bindings'].items():assert task.old.photons.digest(task.old.ROOT/name)==h,name
    previous=spec['prior_accounted_seconds'];signal.alarm(int(240-previous));start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(int(5e9),int(5e9)))
    horizon=json.loads((task.old.OUT/'plan.json').read_text())['horizon_seconds'];m=Model()
    pilot=m.evolve(horizon,8,'condensed-p6-pilot',True)
    agreement=compare(np.load(OUT/'condensed-p6-pilot.npz'),np.load(OUT/'consistent-p6-pilot.npz'))
    assert max(agreement.values())<1e-7,agreement
    forecast=previous+time.monotonic()-start+1.25*(pilot['seconds']*96/8+5)
    task.write(OUT/'order-pilot-budget.json',dict(classification='Counterexample candidate',forecast_total_seconds=forecast,
        setup_seconds=m.setup_seconds,pilot_evolution_seconds=pilot['seconds'],dofs=m.size,
        full_banded_pilot_relative=agreement,
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


def banded():
    assert not (OUT/'order-banded-plan.json').exists()
    spec=json.loads((OUT/'order-permuted-plan.json').read_text());key=str(Path(__file__).relative_to(task.old.ROOT))
    assert task.old.photons.digest(OUT/'order-colamd-source.py')==spec['bindings'][key]
    spec['bindings'][key]=task.old.photons.digest(Path(__file__))
    spec['prior_accounted_seconds']+=10
    spec['third_resource_failure_and_band_diagnosis_allowance_seconds']=10
    measured=json.loads((OUT/'order-bandwidth.json').read_text());assert measured['lower']==50 and measured['upper']==52
    spec['banded_repair']='COLAMD also exhausted the same5GB cap at the first factorization. The measured151612-row matrix has lower/upper bandwidth50/52; use existing LAPACK dgbtrf/dgbtrs with bounded band fill. Equations, row/column scales and all extended-precision residual checks remain identical.'
    task.write(OUT/'order-third-resource-failure.json',dict(classification='Counterexample candidate',passed=False,
        failure='COLAMD first-stage sparse LU allocation failed under unchanged5GB cap.',
        evolution_completed=False,action=spec['banded_repair']))
    task.write(OUT/'order-banded-plan.json',spec)
    (OUT/'order-banded-stage.py').write_bytes(stage_source.encode())


def condensed():
    assert not (OUT/'order-condensed-plan.json').exists()
    spec=json.loads((OUT/'order-banded-plan.json').read_text());key=str(Path(__file__).relative_to(task.old.ROOT))
    assert task.old.photons.digest(OUT/'order-banded-source.py')==spec['bindings'][key]
    spec['bindings'][key]=task.old.photons.digest(Path(__file__))
    pilot=json.loads((OUT/'consistent-p6-pilot.json').read_text())
    spec['prior_accounted_seconds']+=pilot['seconds']+pilot['setup_seconds']+2
    budget=json.loads((OUT/'order-pilot-budget.json').read_text());assert budget['forecast_total_seconds']>=240
    task.write(OUT/'order-budget-stop.json',dict(classification='Counterexample candidate',production_started=False,
        measured_full_band_forecast_seconds=budget['forecast_total_seconds'],budget_seconds=240))
    spec['static_condensation']='Eliminate independent element-internal bubble blocks, retain all endpoint/heat unknowns, and reconstruct every bubble. Compare the identical8-step trajectory against the saved full banded pilot at1e-7 before budgeting production. Full original equation residuals remain checked after reconstruction.'
    spec['pilot_gate']=1e-7
    task.write(OUT/'order-condensed-plan.json',spec)
    (OUT/'order-condensed-stage.py').write_bytes(stage_source.encode())


def review_budget():
    assert not (OUT/'order-budget-review.json').exists()
    spec=json.loads((OUT/'order-condensed-plan.json').read_text())
    pilot=json.loads((OUT/'condensed-p6-pilot.json').read_text())
    measured=json.loads((OUT/'order-pilot-budget.json').read_text())
    assert max(measured['full_banded_pilot_relative'].values())<1e-7
    prior=spec['prior_accounted_seconds']+pilot['seconds']+pilot['setup_seconds']+2
    forecast=prior+1.25*(pilot['seconds']*96/8+pilot['setup_seconds']+5)
    assert forecast<300,forecast
    spec.update(total_budget_seconds=300,previous_budget_seconds=240,prior_accounted_seconds=prior,
        forecast_total_seconds=forecast,
        reason='The240s decision stopped both full-band and condensed production before launch. Reassess once: the remaining previously registered32/64 paths cost about90s from the measured7.27s eight-step pilot. They directly adjudicate the unresolved mass/order ambiguity. Cap the entire phase at300s, with no additional degree, mesh, source tuning or horizon extension.',
        cheaper_alternatives='Reuse background/EOS, all existing p4 paths and both completed p6 pilots. Static condensation preserves the full-band pilot to the recorded tolerance. Run no new pilot.')
    spec['bindings'][str(Path(__file__).relative_to(task.old.ROOT))]=task.old.photons.digest(Path(__file__))
    spec['bindings'][str((OUT/'condensed-p6-pilot.npz').relative_to(task.old.ROOT))]=task.old.photons.digest(OUT/'condensed-p6-pilot.npz')
    task.write(OUT/'order-budget-review.json',spec)


def finish():
    assert not (OUT/'order-result.json').exists();spec=json.loads((OUT/'order-budget-review.json').read_text())
    for name,h in spec['bindings'].items():assert task.old.photons.digest(task.old.ROOT/name)==h,name
    previous=spec['prior_accounted_seconds'];signal.alarm(int(300-previous));start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(int(5e9),int(5e9)))
    horizon=json.loads((task.old.OUT/'plan.json').read_text())['horizon_seconds'];m=Model()
    for n in [32,64]:m.evolve(horizon,n,f'consistent-p6-{n}',True)
    a=np.load(OUT/'consistent-p6-32.npz');b=np.load(OUT/'consistent-p6-64.npz');c=np.load(OUT/'consistent-p4-64.npz')
    time_error=compare(a,{f:b[f][::2] for f in ['temperature','velocity','scalar']})
    space_error=compare(c,b);passed=all(time_error[f]<g and space_error[f]<g for f,g in spec['gates'].items())
    data=dict(classification='Counterexample candidate',passed=bool(passed),time_relative=time_error,space_relative=space_error,
        seconds=time.monotonic()-start,accounted_total_seconds=previous+time.monotonic()-start,total_budget_seconds=300,
        original_GLL_mass_control_passed=False,physical_source_profile_certified=False,moving_surface_solved=False,
        final_charge_solved=False,full_goal_complete=False)
    task.write(OUT/'order-result.json',data);signal.alarm(0);print('CONSISTENT ORDER',json.dumps(data),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['plan','correct','reorder','banded','condensed','review_budget','run','finish']);globals()[parser.parse_args().action]()
