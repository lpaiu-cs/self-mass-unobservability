"""Actual-amplitude nonlinear material defect, without a small-state premise."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,signal,sys,time
import numpy as np
import def_native_collision_material_charge as old

OUT=old.OUT/'finite'/'zero-exact';GR=OUT/'gr';AMP=old.AMP;C=old.previous.prior.C
write=old.write;sha=old.sha;configure=old.configure;LD=np.longdouble


class Material(old.Material):
    def __init__(self,reference,steps=128):
        super().__init__(reference,steps);self.finite_resolution=0.;self.finite_calls=0
    # Evaluate the real donor at the actual amplitude; do not lock a linear branch.
    def raw(self,k,delta,field,eps,row=None):
        self.reference_flux=None
        return old.physical.branch.base.Material.raw(self,k,delta,field,eps,row)
    def rhs(self,t,z,probe=1.):
        j,f,field,rates=self.fields(t);assert not np.any(field) and not np.any(rates)
        if not np.any(z):
            # Identical arguments are the same evaluation, not two uncertain
            # independent fluxes. The initial material defect is exactly zero.
            drive=(self.transfer[j+1]-self.transfer[j])/(self.t[j+1]-self.t[j])
            dt=min(self.point(k)['dt'] for k,w in [(j,1-f),(j+1,f)] if w)
            return drive,np.sum(drive,axis=1),dt
        F=np.zeros((4,self.n+1));G=np.zeros((4,self.n));dt=np.inf;rounding=np.zeros(4)
        for k,w in [(j,1-f),(j+1,f)]:
            if w==0:continue
            p=self.point(k);a=self.raw(k,z,field,AMP)
            F+=w*np.asarray((a[0].astype(LD)-p['flux'].astype(LD))/LD(AMP),float)
            G[1]+=w*np.asarray((a[1].astype(LD)-p['gravity'].astype(LD))/LD(AMP),float);dt=min(dt,a[2])
            # A diagnostic arithmetic floor, not an EOS/recovery error theorem.
            rounding+=w*16*np.finfo(float).eps*np.sum(abs(a[0])+abs(p['flux']),axis=1)/AMP
            rounding[1]+=w*16*np.finfo(float).eps*np.sum(abs(a[1])+abs(p['gravity']))/AMP
            Q=p['Q']+AMP*z*p['active'][None]
            assert np.all(Q[0,p['active']]>0) and np.all(Q[3,p['active']]>0)
            self.finite_calls+=1
        G+=(self.transfer[j+1]-self.transfer[j])/(self.t[j+1]-self.t[j])
        rate=-np.diff(F,axis=1)+G
        self.finite_resolution=max(self.finite_resolution,float(np.max(rounding/np.maximum(np.sum(abs(rate),axis=1),1.))))
        return rate,F[:,0]-F[:,-1]+np.sum(G,axis=1),dt


# Retain the accepted SSP closing-stage convention, CFL and save format.
# Remove the infinitesimal-probe test and its small-state premise explicitly.
source=old.previous.previous.fixed.run_source
begin=source.index('        if k in [1,steps//2-1,steps-1]:');end=source.index('        balance=',begin)
source=source[:begin]+source[end:]
source=source.replace('max(balance)<1e-8 and relative<1e-6 and directional_error<.002','max(balance)<1e-8 and self.finite_resolution<.002')
scope=dict(old.physical.branch.base.__dict__,OUT=OUT);exec(compile(source,__file__,'exec'),scope);Material.run=scope['run']


def prepare():
    assert not OUT.exists();OUT.mkdir();GR.mkdir();gate=json.loads((old.OUT/'gate-result.json').read_text());assert not gate['linear_material_admissible']
    write(OUT.parent/'pilot-failure.json',dict(classification='Counterexample candidate',charged_pilot_seconds=60,
        cause='The gross-rate resolution estimate treated F(Q)-F(Q) at exactly zero defect as two independent uncertain evaluations. Its actual difference is identically zero. Original finite pilot retained.'))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        failure=gate,claim='Replace the inadmissible infinitesimal material approximation by actual-amplitude nonlinear native-interpolation flux/recovery differences, and return their actual finite pressure/trace to GR.',
        equation='delta_Qdot=F(Qbackground+delta_Q)-F(Qbackground)+saved_collision_drive; preserve actual donor, SSP left closing stage, CFL, finite primitives and shared-face conservation. No amplitude reduction. Background remains the saved17-knot trajectory and photons remain the accepted linear collision response.',
        distinction='The original1e-6 linear small-state verdict stays failed. The new finite-amplitude equation does not assume that bound. A16eps cancellation diagnostic is a numerical resolution check, not a rigorous EOS or recovery remainder bound. Do not call this full nonlinear radiation hydrodynamics.',
        budgets=dict(pilot_seconds=60,production_seconds=780,source_seconds=90,GR_seconds=90,CPU_threads=1,virtual_GiB=3),
        reallocation='Charge the first finite pilot60s, then use another60s from the original900s unstarted material production allowance; retain780s for the same two paths. Exactly zero defect returns the exact zero material difference, with unchanged collision drive. No mesh, clock or time-horizon expansion.',
        gates=dict(balance=1e-8,finite_flux_resolution=.002,time=.02,pressure_resolution=.002,quadrature=.002,independent_GR=1e-9),
        forecast='Reuse previous complete late-CFL call counts and fresh two-step finite prefixes; require twice summed forecast below840s. New finite RHS needs one actual state per knot instead of two amplified probes. Resume only the new finite prefixes.',
        stop='On nonlinear primitive, positivity, arithmetic resolution, time or budget failure preserve results and redesign; no automatic longer run or relaxed numerical gates.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(old.__file__),old.OUT/'gate-result.json',old.photons.OUT/'result.json',old.previous.OUT/'production.json']}))
    (OUT/'expanded-run.py').write_text(source)


def worker(n,label,limit=None,restart=None):
    configure();start=time.monotonic();m=Material(128,n);row=m.run(n,label,limit,restart)
    row.update(finite_amplitude=True,finite_flux_resolution=m.finite_resolution,finite_owner_calls=m.finite_calls,
        worker_wall_seconds=time.monotonic()-start,original_linear_gate_passed=row['maximum_true_relative_state']<1e-6,
        directional_check_not_applicable=True,maximum_owner_error=max(p['owner_error'] for p in m.cache.values()))
    row['passed']=row['passed'] and row['maximum_owner_error']<1e-8
    write(OUT/f'{label}.json',row);print(json.dumps(row),flush=True);assert row['passed'];return row


def pilot():
    start=time.monotonic();signal.alarm(60);rows=[worker(n,f'pilot-{n}',2) for n in [64,128]]
    prior=json.loads((old.previous.OUT/'production.json').read_text())['paths'];forecasts=[]
    for r in rows:
        p=next(v for v in prior if (v['steps'],v['reference'])==(r['steps'],128))
        forecasts.append(max(p['raw_owner_calls']*r['seconds']/r['raw_owner_calls']+r['worker_wall_seconds']-r['seconds']+15,p['worker_CPU_seconds']))
    result=dict(classification='Counterexample candidate',rows=rows,forecast_each_seconds=forecasts,upper_seconds=2*sum(forecasts),
        eligible=2*sum(forecasts)<780,seconds=time.monotonic()-start)
    write(OUT/'pilot.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def production():
    assert json.loads((OUT/'pilot.json').read_text())['eligible'];start=time.monotonic();signal.alarm(780)
    rows=[worker(n,f'steps-{n}-reference-128',restart=f'pilot-{n}') for n in [64,128]]
    write(OUT/'production.json',dict(classification='Counterexample candidate',passed=True,rows=rows,seconds=time.monotonic()-start));signal.alarm(0)


def pressure(m,k,z,field):
    assert not np.any(field);p=m.point(k);model=m.model;b=model.bulk;f=model.flow;nb=m.nb
    def value(delta,eps):
        row=m.raw(k,delta,field,eps);rr,v,lt,y=row[3]['primitive'];f.eos.y=y
        pg,u,*_=f.eos(rr,lt);pg=pg*f.eos.rho0*C*C
        H=(rr*(f.eos.cx+u)*f.eos.rho0*C*C+pg).astype(LD)
        gas=np.r_[b.eos.gas(row[3]['theta'],row[3]['eta'])[0],pg].astype(LD)*m.V
        radial=gas.copy();radial[:nb]+=2*model.kinetic();radial[nb:]+=H*(v.astype(LD)**2/(1-v.astype(LD)**2))*m.V[nb:]
        return np.array([gas,radial])
    before=value(np.zeros_like(z),0.);after=value(z,AMP);delta=(after-before)/LD(AMP)
    floor=16*np.finfo(float).eps*(abs(after)+abs(before))/AMP
    return np.asarray(before,float),np.asarray(delta,float),np.asarray(floor,float)


def sources():
    assert json.loads((OUT/'production.json').read_text())['passed'];configure()
    s=old.previous.previous.old.source_owner.source.replace('signal.alarm(75)','signal.alarm(90)')
    s=s.replace('[[64,128],[128,128],[128,64]]','[[64,128],[128,128]]').replace('background=compare(histories[2],histories[1]),','').replace(',stress_background=compare(allstress[2],allstress[1])','')
    scope=dict(old.previous.source_scope,OUT=OUT,GR=GR,Material=Material,stress=SimpleNamespace(pressure=pressure),
        photons=SimpleNamespace(path=old.photon_path),material_path=lambda n,r:OUT/f'steps-{n}-reference-{r}.npz')
    exec(compile(s,__file__,'exec'),scope);(OUT/'expanded-source.py').write_text(s);scope['sources']()
    assert json.loads((OUT/'sources.json').read_text())['passed']


charge=FunctionType(old.charge.__code__,dict(old.charge.__globals__,OUT=OUT,GR=GR,__file__=__file__))


if __name__=='__main__':
    signal.signal(signal.SIGALRM,old.repaired.forcing.history.flow.old.optical.timeout);action=sys.argv[1];started=time.monotonic()
    try:globals()[action]()
    except Exception as exc:write(OUT/f'{action}-failure.json',dict(error=repr(exc),seconds=time.monotonic()-started));raise
