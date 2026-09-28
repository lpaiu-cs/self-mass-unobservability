"""Keep the closing SSP stage on its own source interval at a stored knot."""
from pathlib import Path
from types import FunctionType
import inspect,json,sys,time
import numpy as np
import sympy as sp
import verify_native_feedback_metric_return as base

OUT=base.OUT/'stage-aligned';GR=OUT/'gr';write=base.write;sha=base.sha
prior=base.prior;matter=base.matter
old_production=json.loads((base.OUT/'production.json').read_text())
old_pilot=json.loads((base.OUT/'pilot.json').read_text())
old_sources=json.loads((base.OUT/'sources.json').read_text())
production_cap=int(350-old_production['seconds']);pilot_cap=int(45-old_pilot['seconds']);source_cap=int(75-old_sources['seconds'])


field_source=prior.replace(inspect.getsource(matter.branch.base.Material.fields),"side='right'","side='left' if self.closing_stage else 'right'")
import textwrap
ns=dict(vars(matter.branch.base));exec(compile(textwrap.dedent(field_source),__file__,'exec'),ns)
run_source=textwrap.dedent(inspect.getsource(matter.branch.base.Material.run))
run_source=prior.replace(run_source,'r2,l2,cfl2=self.rhs(t+dt,trial)','r2,l2,cfl2=self.end_rhs(t+dt,trial)')
run_scope=dict(vars(matter.branch.base),OUT=OUT);exec(compile(run_source,__file__,'exec'),run_scope)


class Material(base.Material):
    closing_stage=False
    interval_fields=ns['fields']
    run=run_scope['run']

    def fields(self,t):
        # Arithmetic paths to the same declared knot may differ by a few ulps.
        k=int(np.argmin(abs(self.t-t)))
        if abs(self.t[k]-t)<=16*np.finfo(float).eps*self.t[-1]:t=self.t[k]
        return self.interval_fields(t)

    def end_rhs(self,t,z):
        self.closing_stage=True
        try:return self.rhs(t,z)
        finally:self.closing_stage=False


dispatch=FunctionType(matter.parallel.dispatch.__code__,dict(vars(matter.parallel),OUT=OUT,__file__=__file__))
worker=FunctionType(base.worker.__code__,dict(base.worker.__globals__,OUT=OUT,Material=Material))
pilot_source=inspect.getsource(matter.pilot).replace("previous.OUT/'material-audit.json'","base.OUT/'production.json'").replace('],45)','],pilot_cap)').replace('2*max(forecasts)<350','2*max(forecasts)<production_cap').replace('hard_seconds=350','hard_seconds=production_cap')
scope=dict(vars(matter),OUT=OUT,Material=Material,dispatch=dispatch,__file__=__file__,base=base,pilot_cap=pilot_cap,production_cap=production_cap)
exec(compile(pilot_source,__file__,'exec'),scope);pilot=scope['pilot']
production_source=prior.replace(inspect.getsource(matter.production),'for n,r in PATHS],350)','for n,r in PATHS],production_cap)')
exec(compile(production_source,__file__,'exec'),scope);production=scope['production']
source=prior.replace(base.source,'signal.alarm(75)','signal.alarm(source_cap)')
source_scope=dict(base.source_scope,OUT=OUT,GR=GR,Material=Material,source_cap=source_cap,material_path=lambda n,r:OUT/f'steps-{n}-reference-{r}.npz')
exec(compile(source,__file__,'exec'),source_scope)


class GRResponse(base.GRResponse):
    run=FunctionType(matter.previous.wave.base.Response.run.__code__,dict(vars(matter.previous.wave.base),OUT=GR))


fields_owner=FunctionType(base.fields_owner.__code__,dict(base.fields_owner.__globals__,OUT=OUT,GR=GR,GRResponse=GRResponse,__file__=__file__))
audit_owner=FunctionType(base.audit.__code__,dict(vars(base),OUT=OUT,GR=GR,__file__=__file__))


def prepare():
    assert not old_sources['passed'] and not OUT.exists();OUT.mkdir();GR.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        failure='Original material pressure time comparison3.17399percent exceeds2percent. Energy/H state comparisons pass. The unchanged primitive-pressure probe passes, so no arbitrary pressure rescaling or relaxed gate is allowed.',
        defect='At a stored source knot, the original SSP closing stage uses searchsorted right and consumes the next interval derivative. The starting stage needs the right interval; the closing stage needs the left interval. Metric-work and cumulative photon-transfer slopes are discontinuous there.',
        repair='Change only closing-stage interval selection and canonical-knot roundoff alignment. Keep the same shared fluxes,SSP2,CFL,macro clocks,paths and gates. Reuse all photon and GR histories. Preserve the failed material/source files.',
        budgets=dict(pilot_remaining_seconds=pilot_cap,production_remaining_seconds=production_cap,source_remaining_seconds=source_cap,GR_seconds=120,CPU_processes=3,threads_each=1,total_virtual_GiB=6),
        forecast='Use original completed Phase135 late-CFL counts and fresh concurrent prefixes. Require2x forecast within the original350s budget remainder. No automatic resolution or period increase.',
        stop='Stop at failed repaired comparison or remaining budget. Do not proceed to GR on a failed source result.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(base.__file__),Path(base.run.__file__),base.OUT/'production.json',base.OUT/'sources.json',base.OUT/'plan.json',base.run.RESPONSE/'result.json']}))
    a,b,h=sp.symbols('a b h');wrong=h*(a+b)/2;correct=h*a
    assert sp.simplify(wrong-correct-h*(b-a)/2)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='For a piecewise-constant derivative a before a knot and b after it, a step ending at the knot has integral h*a. Using the right-limit closing stage adds h*(b-a)/2. The left-limit closing stage removes that defect. This does not itself prove the full pressure comparison will pass.'))
    started=time.monotonic();matter.configure();m=Material(128)
    d=np.load(base.OUT/'steps-128-reference-128.npz');t=float(m.t[15]);i=int(np.argmin(abs(d['t']-t)));z=d['history_scaled'][i]
    right=m.rhs(t,z)[0];jr,_,_,rr=m.fields(t)
    left=m.end_rhs(t,z)[0];m.closing_stage=True;jl,_,_,rl=m.fields(t);m.closing_stage=False
    assert jr==15 and jl==14
    drive=lambda j:(m.transfer[j+1]-m.transfer[j])/(m.t[j+1]-m.t[j])
    expected=drive(jl)-drive(jr);point=m.point(15);rates=rl-rr
    expected[1]-=(rates[0]+rates[1])*point['Q'][1]
    expected[2]-=m.a*(point['Pr']*(rates[0]+rates[1])+2*point['Pg']*rates[0])
    norm=np.maximum.reduce([np.sum(abs(expected),axis=1),np.sum(abs(left),axis=1),np.sum(abs(right),axis=1),np.ones(4)*1e-300])
    err=float(np.max(np.sum(abs(left-right-expected),axis=1)/norm))
    write(OUT/'forcing-jump.json',dict(classification='Counterexample candidate',passed=err<1e-8,
        time=t,actual_jump_L1_physical_per_second=(np.sum(abs(left-right),axis=1)*matter.AMP).tolist(),
        predicted_jump_relative=err,seconds=time.monotonic()-started))
    assert err<1e-8


def sources():
    matter.configure();assert json.loads((OUT/'production.json').read_text())['passed']
    write(OUT/'source-plan.json',dict(classification='Counterexample candidate',budget_seconds=source_cap,
        claim='Re-evaluate actual pressure/trace and unchanged2percent source comparisons after the SSP endpoint repair; preserve original failed source result.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(base.__file__),OUT/'production.json',base.OUT/'sources.json']}))
    source_scope['sources']()


def fields():
    matter.configure();fields_owner()


def audit():
    audit_owner()
    d=json.loads((OUT/'audit.json').read_text());d.update(original_pressure_time_verdict=False,closing_stage_uses_left_source_interval=True,
        total_material_pilot_and_jump_seconds=old_pilot['seconds']+json.loads((OUT/'pilot.json').read_text())['seconds']+json.loads((OUT/'forcing-jump.json').read_text())['seconds'],
        total_material_production_seconds=old_production['seconds']+json.loads((OUT/'production.json').read_text())['seconds'],
        total_source_seconds=old_sources['seconds']+json.loads((OUT/'sources.json').read_text())['seconds'])
    assert d['total_material_production_seconds']<350 and d['total_source_seconds']<75 and d['total_material_pilot_and_jump_seconds']<45
    write(OUT/'audit.json',d)


if __name__=='__main__':
    if sys.argv[1]=='worker':worker(int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],None if sys.argv[5]=='None' else int(sys.argv[5]),None if sys.argv[6]=='None' else sys.argv[6])
    else:globals()[sys.argv[1]]()
